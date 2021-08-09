from os import path, remove
import pandas as pd
import requests
from xmltodict import parse
import urllib
from urllib.request import urlretrieve

from qebil.log import logger
from qebil.tools.fastq import (
    blast_for_type,
    get_read_count,
)

from qebil.tools.util import (
    get_checksum,
    unpack_fastq_ftp,
    remove_index_read_file,
)


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
    url = "https://www.ebi.ac.uk/ena/browser/api/xml/" + str(accession)
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

    host = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession="
    read_type = "&result=read_run&limit=0&offset=0&download=true&"
    # TODO: consider passing max-samples in here for limit to speed up process
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
    run_prefix, ftp_dict, output_dir, lib_layout, remove_index_file=False
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
    lib_layout: string
        layout of library; single or paired
    remove_index_file: bool
        whether to try to identify and remove any index file
        based on read length disparity

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
        local_read_dict["read" + str(read_num)] = {}
        local_read_dict["read" + str(read_num)]["fp"] = local_fq_path

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
            removed_index_dict = remove_index_read_file(
                local_read_dict, lib_layout
            )
            if len(removed_index_dict) < len(local_read_dict):
                for f in local_read_dict.keys():
                    local_read_dict = removed_index_dict

        if len(local_read_dict) == 1:
            raw_reads = get_read_count(local_read_dict["read1"]["fp"])
            if raw_reads == "fqtools error":
                logger.warning(
                    "Check file validity failed for "
                    + str(local_read_dict["read1"])
                    + "Removing."
                )
                remove(local_read_dict["read1"]["fp"])
                raw_reads = "error"
        elif len(local_read_dict) == 2:
            raw_reads = get_read_count(
                local_read_dict["read1"]["fp"], local_read_dict["read2"]["fp"]
            )
            if raw_reads == "fqtools error":
                logger.warning(
                    "Check file validity failed for "
                    + str(local_read_dict["read1"]["fp"].replace("R1", "R*"))
                    + "Removing."
                )
                remove(local_read_dict["read1"]["fp"])
                remove(local_read_dict["read2"]["fp"])
                raw_reads = "error"
        else:
            logger.warning(
                "Read count not possible for"
                + str(len(local_read_dict))
                + " reads."
            )
            raw_reads = "error"

    return raw_reads, local_read_dict


def fetch_fastqs(study, output_dir, remove_index_file=False, overwrite=False):
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
            if "local_fastq_fp" in row:
                logger.info(
                    "Using local filepaths and checksums."
                    + " Pass overwrite=True to force retrieval from ftp"
                )
                fp_string = str(row["local_fastq_fp"])
                md5_string = str(row["local_fastq_fp"])
            else:
                fp_string = str(row["fastq_ftp"])
                md5_string = str(row["fastq_md5"])

            logger.info("Unpacking: " + fp_string)
            ebi_dict = unpack_fastq_ftp(fp_string, md5_string)
            if len(ebi_dict) == 0:
                logger.warning(
                    "No fastq files to download found for\n" + run_prefix
                )
            else:
                layout = row["library_layout"]
                fetch_result = fetch_fastq_files(
                    run_prefix,
                    ebi_dict,
                    output_dir,
                    layout,
                    remove_index_file,
                )
                row["qebil_raw_reads"] = fetch_result[0]
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

                    corrected_fq_dict = fetch_result[1]
                    if len(corrected_fq_dict) < len(ebi_dict):
                        logger.warning(
                            "Number of read files changed, updating metadata."
                        )
                        fp_list = []
                        md5_list = []
                        for read in corrected_fq_dict.keys():
                            read_num = read[-1]
                            if read_num == "1":
                                fp_list = (
                                    corrected_fq_dict[read]["fp"] + fp_list
                                )
                                md5_list = (
                                    corrected_fq_dict[read]["md5"] + md5_list
                                )
                            elif read_num == "2":
                                fp_list += corrected_fq_dict[read]["fp"]
                                md5_list += corrected_fq_dict[read]["md5"]
                            else:
                                logger.error(
                                    "Should not be here."
                                    + " corrected_fq_dict is diff and >=3."
                                )

                        row["local_fastq_fp"] = ";".join(fp_list)
                        row["local_fastq_md5"] = ";".join(md5_list)
                    else:
                        row["local_fastq_fp"] = row["fastq_ftp"]
                        row["local_fastq_md5"] = row["fastq_md5"]

                    # adding loop to try to resolve ambiguous hits
                    # could be improved, but first stab at it
                    strategy = row["library_strategy"]
                    prep_type = row["qebil_prep_file"].split("_")[1]
                    if (
                        strategy in ["AMPLICON", "OTHER"]
                        and prep_type == "AMBIGUOUS"
                    ):
                        logger.warning(
                            "Attempting to resolve AMBIGUOUS"
                            + " assignment by blasting"
                        )
                        if (
                            len(ebi_dict) <= 2
                        ):  # for single/paired data, take R1. If 3, take R2
                            local_fq_path = (
                                output_dir
                                + "/"
                                + run_prefix
                                + ".R1.ebi.fastq.gz"
                            )
                        elif len(ebi_dict) == 3:
                            local_fq_path = (
                                output_dir
                                + "/"
                                + run_prefix
                                + ".R2.ebi.fastq.gz"
                            )
                        blast_type = blast_for_type(local_fq_path)
                        row["qebil_prep_file"] = row[
                            "qebil_prep_file"
                        ].replace(prep_type, blast_type)
        # even if it fails, we'll need to return the row
        row_list.append(new_row)

    # combine the rows back together
    md = pd.DataFrame()
    for r in row_list:
        md = md.append(r)

    return md


def retrieve_ftp_file(ftp_path, filepath, remote_checksum, overwrite=False):
    """Method to retrieve an ftp file and check accuracy
    of the download. If not overwriting, check the local copy for validity
    before downloading.

    Parameters
    ----------
    ftp_path: string:
        the ftp url to download from
    filepath: string
        the local path to save the file to
    remote_checksum: string
        hexadecimal md5 checksum to validate the downlod
    overwrite : bool
        whether to overwrite the local copy of the file

    Returns
    ---------
    checksum: str or bool
        returns either the string for the checksum or False
        if there is an issue with the download or integrity

    """
    # see if the file exists and make sure it's valid if so
    local = False
    if not overwrite:
        if path.isfile(filepath):
            checksum = get_checksum(filepath, remote_checksum)
            if checksum:
                local = True
                logger.info(
                    "Valid file found. Skipping download of file: "
                    + "ftp://"
                    + ftp_path
                    + " to "
                    + filepath
                )
                return checksum
            else:
                logger.warning(
                    "Local file found but invalid checksum."
                    + " Downloading again: "
                    + "ftp://"
                    + ftp_path
                    + " to "
                    + filepath
                )

    if not local:
        # add catch in case there is an issue with the connection
        try:
            urlretrieve("ftp://" + ftp_path, filepath)
            checksum = get_checksum(filepath, remote_checksum)
            return checksum
        except urllib.error.URLError:  # RequestException:
            logger.warning(
                "Issue with urlretrieve for "
                + "download of file:"
                + "ftp://"
                + ftp_path
                + " to "
                + filepath
            )
            # cleanup file if present
            if path.isfile(ftp_path):
                remove(ftp_path)
            return False
