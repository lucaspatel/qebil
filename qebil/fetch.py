from os import path, remove
import pandas as pd
import requests
from xmltodict import parse

from qebil.log import logger
from qebil.tools.fastq import (
    unpack_fastq_ftp,
    get_read_count,
    remove_index_read_file,
)
from qebil.tools.util import retrieve_ftp_file, retrieve_aspera_file


def fetch_ebi_info(accession):
    """Retrieve metadata from EBI

    Parameters
    ----------
    accession:string
        accession name from EBI/ENA

    Returns
    -------
    xml_dict: dict
       dictionary of accession information as provided by EBI

    """
    xml_dict = {}
    # could use this list valid_stems=["PRJEB", "PRJNA", "ERP",
    # "SRP", "SRR", "SRS"] to check format before requesting url
    url = (
        "http://www.ebi.ac.uk/ena/data/view/"
        + str(accession)
        + "&display=xml"
    )
    logger.info(url)
    try:
        response = requests.get(url)
        xml_dict = parse(response.content)
    except Exception:
        # TODO: add response exception, and separate out parse
        logger.error(
            "Could not obtain EBI information for "
            + str(accession)
            + " at URL: "
            + url
            + " . Please check connection and"
            + " accession id and try again."
        )

    return xml_dict


def fetch_ebi_metadata(study_accession, fields=[]):
    """Retrieves basic metadata information for a provided study

    Parameters
    ----------
    study_accession : pd.DataFrame
        string of study for retireval
    fields : list
        list of metadata columns to retrieve for samples
        If empty, will use default list

    Returns
    -------
    study_df: pd.DataFrame
       dataframe of study information as provided by the chosen repository
    """

    if len(fields) == 0:
        fields = [
            "sample_accession",
            "library_name",
            "secondary_sample_accession",
            "run_accession",
            "experiment_accession",
            "fastq_ftp",
            "fastq_aspera",
            "library_source",
            "instrument_platform",
            "submitted_format",
            "library_strategy",
            "library_layout",
            "tax_id",
            "scientific_name",
            "instrument_model",
            "library_selection",
            "center_name",
            "experiment_title",
            "study_title",
            "study_alias",
            "experiment_alias",
            "sample_alias",
            "sample_title",
            "fastq_md5",
            "study_accession",
        ]

    study_df = pd.DataFrame()

    host = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession="
    read_type = "&result=read_run&"
    fields = ",".join(fields)

    url = "".join([host, study_accession, read_type, "fields=", fields])
    logger.info("Querying study metadata at: " + url)

    try:
        study_df = pd.read_csv(url, sep="\t", dtype=str)
        study_df.dropna(axis=1, how="all", inplace=True)
    except Exception:
        logger.warning(
            "Issue retrieving study information for "
            + study_accession
            + " with url "
            + url
        )
    return study_df


def fetch_fastq_files(
    run_prefix, ftp_dict, output_dir, remove_index_file=False
):
    """Retrieves fastq files from EBI, checks their validity by
    comparing the remote and local checksums and then returns
    the number of reads in the downloaded file(s).

    Parameters
    -----------
    run_prefix: string
        prefix prepended to local files
    ftp_dict: dict
        dictionary of remote ftp files and checksums
    output_dir: string
        directory where to write the output file(s)

    Returns
    ---------
    raw_reads : int
        the number of reads in the downloaded file(s)
    """

    failed_list = []
    local_read_dict = {}

    for read in ftp_dict.keys():
        read_num = read[-1]
        local_fq_path = (
            output_dir
            + "/"
            + run_prefix
            + ".R"
            + str(read_num)
            + ".ebi.fastq.gz"
        )
        remote_fp = ftp_dict[read]["ftp"]
        remote_md5 = ftp_dict[read]["md5"]
        local_read_dict["read" + str(read_num)] = local_fq_path
        if run_prefix in failed_list:
            logger.warning(
                "Skipping download of " + remote_fp + " Paired file failed."
            )
        else:
            success = retrieve_ftp_file(remote_fp, local_fq_path, remote_md5)
            if not success:
                logger.warning(
                    "Download of "
                    + remote_fp
                    + " to "
                    + local_fq_path
                    + " failed."
                )
                failed_list.append(run_prefix)

    if len(failed_list) > 0:
        raw_reads = "error"
        logger.warning(
            "The fastq file(s) failed to download for:"
            + str(failed_list)
            + " removing partial"
            + " fastq files."
        )
        for r in local_read_dict.values():
            if path.isfile(r):
                logger.info("Removing " + r)
                remove(r)
            else:
                logger.warning(
                    "Could not remove" + r + "file does not exist."
                )
    else:
        if remove_index_file:
            logger.info("Removing index file")
            remove_index_read_file(run_prefix)

        if len(local_read_dict) == 1:
            raw_reads = get_read_count(local_read_dict["read1"])
            if raw_reads == "fqtools error":
                logger.warning(
                    "Check file validity failed for "
                    + str(local_read_dict["read1"])
                    + "Removing."
                )
                remove(local_read_dict["read1"])
                raw_reads = "error"
        elif len(local_read_dict) == 2:
            raw_reads = get_read_count(
                local_read_dict["read1"], local_read_dict["read2"]
            )
            if raw_reads == "fqtools error":
                logger.warning(
                    "Check file validity failed for "
                    + str(local_read_dict["read1"].replace("R1", "R*"))
                    + "Removing."
                )
                remove(local_read_dict["read1"])
                remove(local_read_dict["read2"])
                raw_reads = "error"
        else:
            logger.warning(
                "Read count not possible for"
                + str(len(local_read_dict))
                + " reads."
            )
            raw_reads = "error"

    return raw_reads


def fetch_fastq_files_aspera(
    run_prefix, ftp_dict, output_dir, remove_index_file=False
):
    """Retrieves fastq files from EBI, checks their validity by
    comparing the remote and local checksums and then returns
    the number of reads in the downloaded file(s).

    Parameters
    -----------
    run_prefix: string
        prefix prepended to local files
    ftp_dict: dict
        dictionary of remote ftp files and checksums
    output_dir: string
        directory where to write the output file(s)

    Returns
    ---------
    raw_reads : int
        the number of reads in the downloaded file(s)
    """

    failed_list = []
    local_read_dict = {}

    for read in ftp_dict.keys():
        read_num = read[-1]
        local_fq_path = (
            output_dir
            + "/"
            + run_prefix
            + ".R"
            + str(read_num)
            + ".ebi.fastq.gz"
        )
        remote_fp = ftp_dict[read]["ftp"]
        remote_md5 = ftp_dict[read]["md5"]
        local_read_dict["read" + str(read_num)] = local_fq_path
        if run_prefix in failed_list:
            logger.warning(
                "Skipping download of " + remote_fp + " Paired file failed."
            )
        else:
            success = retrieve_aspera_file(
                remote_fp, local_fq_path, remote_md5
            )
            if not success:
                logger.warning(
                    "Download of "
                    + remote_fp
                    + " to "
                    + local_fq_path
                    + " failed."
                )
                failed_list.append(run_prefix)

    if len(failed_list) > 0:
        raw_reads = "error"
        logger.warning(
            "The fastq file(s) failed to download for:"
            + str(failed_list)
            + " removing partial"
            + " fastq files."
        )
        for r in local_read_dict.values():
            if path.isfile(r):
                logger.info("Removing " + r)
                remove(r)
            else:
                logger.warning(
                    "Could not remove" + r + "file does not exist."
                )
    else:
        if remove_index_file:
            logger.info("Removing index file")
            remove_index_read_file(run_prefix)

        if len(local_read_dict) == 1:
            raw_reads = get_read_count(local_read_dict["read1"])
            if raw_reads == "fqtools error":
                logger.warning(
                    "Check file validity failed for "
                    + str(local_read_dict["read1"])
                    + "Removing."
                )
                remove(local_read_dict["read1"])
                raw_reads = "error"
        elif len(local_read_dict) == 2:
            raw_reads = get_read_count(
                local_read_dict["read1"], local_read_dict["read2"]
            )
            if raw_reads == "fqtools error":
                logger.warning(
                    "Check file validity failed for "
                    + str(local_read_dict["read1"].replace("R1", "R*"))
                    + "Removing."
                )
                remove(local_read_dict["read1"])
                remove(local_read_dict["read2"])
                raw_reads = "error"
        else:
            logger.warning(
                "Read count not possible for"
                + str(len(local_read_dict))
                + " reads."
            )
            raw_reads = "error"

    return raw_reads


def fetch_fastqs(study, output_dir, remove_index_file=False):
    """Helper method to consolidate calls to other methods for processing fastqs

    This method is potentially convoluted to try to handle the fact that runs
    need to be able to recover and some will want to minimize the storage
    footprint for heavily host contaminated files.

    This may be done more simply and elegantly so consider refactor?

    Parameters
    ----------
    study: qebil.Study
        Study object
    output_dir: string
        where to write the files
    remove_index_file: bool
        whether to resolve studies with three reads to remove the index file

    """
    md = study.metadata
    row_list = []

    for index, row in md.iterrows():
        run_prefix = ""
        ebi_dict = {}
        new_row = row
        try:
            run_prefix = row["run_prefix"]
        except KeyError:
            logger.warning("No run_prefix in metadata for: " + index)
        if len(run_prefix) > 0:
            logger.info("Unpacking: " + str(row["fastq_ftp"]))
            ebi_dict = unpack_fastq_ftp(
                str(row["fastq_ftp"]), str(row["fastq_md5"])
            )
            if len(ebi_dict) == 0:
                logger.warning(
                    "No fastq files to download found for\n" + run_prefix
                )
            else:
                row["qebil_raw_reads"] = fetch_fastq_files(
                    run_prefix, ebi_dict, output_dir, remove_index_file
                )
                if row["qebil_raw_reads"] == "error":
                    logger.warning("Issue retrieving files for " + run_prefix)
                else:
                    logger.info(
                        "Retrieved "
                        + run_prefix
                        + " with "
                        + str(row["qebil_raw_reads"])
                        + " reads."
                    )
        # even if it fails, we'll need to return the row
        row_list.append(new_row)

    # combine the rows back together
    md = pd.DataFrame()
    for r in row_list:
        md = md.append(r)

    return md


def fetch_fastqs_aspera(study, output_dir, remove_index_file=False):
    """Helper method to consolidate calls to other methods for processing fastqs

    This method is potentially convoluted to try to handle the fact that runs
    need to be able to recover and some will want to minimize the storage
    footprint for heavily host contaminated files.

    This may be done more simply and elegantly so consider refactor?

    Parameters
    ----------
    study: qebil.Study
        Study object
    output_dir: string
        where to write the files
    remove_index_file: bool
        whether to resolve studies with three reads to remove the index file

    """
    md = study.metadata
    row_list = []

    for index, row in md.iterrows():
        run_prefix = ""
        ebi_dict = {}
        new_row = row
        try:
            run_prefix = row["run_prefix"]
        except KeyError:
            logger.warning("No run_prefix in metadata for: " + index)
        if len(run_prefix) > 0:
            logger.info("Unpacking: " + str(row["fastq_aspera"]))
            ebi_dict = unpack_fastq_ftp(
                str(row["fastq_aspera"]), str(row["fastq_md5"])
            )
            if len(ebi_dict) == 0:
                logger.warning(
                    "No fastq files to download found for\n" + run_prefix
                )
            else:
                row["qebil_raw_reads"] = fetch_fastq_files_aspera(
                    run_prefix, ebi_dict, output_dir, remove_index_file
                )
                if row["qebil_raw_reads"] == "error":
                    logger.warning("Issue retrieving files for " + run_prefix)
                else:
                    logger.info(
                        "Retrieved "
                        + run_prefix
                        + " with "
                        + str(row["qebil_raw_reads"])
                        + " reads."
                    )
        # even if it fails, we'll need to return the row
        row_list.append(new_row)

    # combine the rows back together
    md = pd.DataFrame()
    for r in row_list:
        md = md.append(r)

    return md
