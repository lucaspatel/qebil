from os import path, remove
import pandas as pd
import requests
from xmltodict import parse
import urllib
from urllib.request import urlretrieve
from time import sleep

from qebil.log import logger
from qebil.tools.fastq import (
    blast_for_type,
    get_read_count,
)

from qebil.tools.util import (
    get_checksum,
    unpack_fastq_ftp,
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
        logger.debug(f"URL is: {url}")
        xml_dict = parse(response.content)
        logger.debug(f"XML dict is: {xml_dict}")
    except Exception:
        # TODO: add response exception, and separate out parse
        logger.error(
            f"Could not obtain EBI information for {accession} at URL: {url}. Please check connection and accession id and try again."
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
            "fastq_bytes",
            "fastq_ftp",
            "fastq_md5",
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


def fetch_fastq_files(run_prefix, ftp_dict, output_dir, expected_reads=""):
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
    expected_reads: string
        number of reads in expected in the downloaded file

    Returns
    ---------
    raw_reads : int
        the number of reads in the downloaded file(s)
    """

    failed_list = []
    local_read_dict = {}
    qebil_notes = ""

    for read in ftp_dict.keys():
        skip = False
        read_num = read[-1]
        # to make backward compatible, check that run_prefix ends in .R
        if run_prefix[-2:] != ".R":
            run_prefix = run_prefix + ".R"

        if read_num == 0:
            local_fq_path = (
                output_dir
                + "/"
                + run_prefix.replace(".R", "_R")
                + str(read_num)
                + ".ebi.fastq.gz"
            )
        else:
            local_fq_path = (
                output_dir
                + "/"
                + run_prefix
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
            skip = True
        else:
            if path.isfile(local_fq_path):
                # local_md5 = get_checksum(local_fq_path)
                # local_read_dict["read" + str(read_num)]["md5"] = local_md5
                local_md5 = remote_md5

                if remote_md5 == local_md5:
                    logger.info(
                        "Valid file found. Skipping download of file: "
                        + "ftp://"
                        + remote_fp
                        + " to "
                        + local_fq_path
                    )
                    skip = True
                else:
                    logger.warning(
                        "Local file found but invalid checksum."
                        + " Downloading again: "
                        + "ftp://"
                        + remote_fp
                        + " to "
                        + local_fq_path
                    )
        if not skip:
            success = retrieve_ftp_file(remote_fp, local_fq_path, remote_md5)
            if success:
                local_read_dict["read" + str(read_num)]["md5"] = success
            else:
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
            read_fp = r["fp"]
            if path.isfile(read_fp):
                logger.info("Removing " + read_fp)
                remove(read_fp)
            else:
                logger.warning(
                    "Could not remove" + read_fp + "file does not exist."
                )
    else:
        if expected_reads.isnumeric():
            raw_reads = (
                expected_reads  # skips getting read count if already obtained
            )
        else:
            if "read2" in local_read_dict.keys():
                raw_reads = get_read_count(
                    local_read_dict["read1"]["fp"],
                    local_read_dict["read2"]["fp"],
                )
                if raw_reads == "fqtools error":
                    logger.warning(
                        "Check file validity failed for "
                        + str(
                            local_read_dict["read1"]["fp"].replace("R1", "R*")
                        )
                        + "Removing."
                    )
                    remove(local_read_dict["read1"]["fp"])
                    remove(local_read_dict["read2"]["fp"])
                    raw_reads = "error"
                    qebil_notes += "fastq file corrupted"
            else:
                raw_reads = get_read_count(local_read_dict["read1"]["fp"])
                if raw_reads == "fqtools error":
                    logger.warning(
                        "Check file validity failed for "
                        + str(local_read_dict["read1"])
                        + "Removing."
                    )
                    remove(local_read_dict["read1"]["fp"])
                    raw_reads = "error"
                    qebil_notes += "fastq file corrupted"

    return raw_reads, local_read_dict, qebil_notes


def fetch_fastqs(study, output_dir, overwrite=False):
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
    overwrite: boolean
        whether to overwrite the existing files

    """
    md = study.metadata
    row_list = []
    layout_dict = {"single": 1, "paired": 2}

    for index, row in md.iterrows():
        run_prefix = ""
        ebi_dict = {}
        new_row = row
        skip = False

        try:
            run_prefix = row["run_prefix"]
        except KeyError:
            logger.warning("No run_prefix in metadata for: " + index)
            skip = True

        try:
            notes = row["qebil_notes"]
        except KeyError:
            logger.info("No qebil_notes in metadata for: " + index)
            notes = ""
        if notes != "":
            if notes == "fastq file corrupted" and not overwrite:
                logger.warning(
                    "Previous download of files for "
                    + run_prefix
                    + " failed, skipping download."
                    + " Pass overwrite=True to override."
                )
                skip = True

        if len(run_prefix) > 0 and not skip:
            logger.debug("Row is: " + str(row))
            layout = str(row["library_layout"]).lower()
            if layout in layout_dict.keys():
                expected_num_read_files = layout_dict[layout]

            fp_string = str(row["fastq_ftp"])
            md5_string = str(row["fastq_md5"])
            bytes_string = str(row["fastq_bytes"])
            ebi_dict, error = unpack_fastq_ftp(
                fp_string, md5_string, bytes_string, expected_num_read_files
            )
            logger.debug("ebi_dict: " + str(ebi_dict))
            if len(ebi_dict) == 0:
                logger.warning(
                    "No fastq files to download found for\n" + run_prefix
                )
                skip = True
            elif error != "":
                logger.warning(error)
                if (
                    error == "More than 3 read files in ftp, skipping"
                    or "Fewer read files than expected in ftp."
                ):
                    skip = True

            if not skip:
                raw_reads = ""
                try:
                    raw_reads = row["qebil_raw_reads"]
                except KeyError:
                    logger.info(
                        "qebil_raw_reads not found for run prefix"
                        + " not expecting local file."
                    )
                read_count, local_dict, notes = fetch_fastq_files(
                    run_prefix, ebi_dict, output_dir, raw_reads
                )

                row["qebil_raw_reads"] = read_count
                row["qebil_notes"] = notes
                if row["qebil_raw_reads"] == "error":
                    logger.warning(
                        "Issue retrieving files for "
                        + run_prefix
                        + ": "
                        + notes
                    )
                else:
                    logger.info(
                        "Retrieved "
                        + run_prefix
                        + " with "
                        + str(row["qebil_raw_reads"])
                        + " reads."
                    )

                    # adding loop to try to resolve ambiguous hits
                    # could be improved, but first stab at it
                    strategy = row["library_strategy"]
                    prep_type = row["qebil_prep_file"].split("_")[1]
                    if (
                        strategy in ["AMPLICON", "OTHER"]
                        and prep_type == "AMBIGUOUS"
                    ):
                        blast_type = blast_for_type(local_dict["read1"]["fp"])
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


def retrieve_ftp_file(ftp_path, filepath, remote_checksum):
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

    Returns
    ---------
    checksum: str or bool
        returns either the string for the checksum or False
        if there is an issue with the download or integrity

    """
    # see if the file exists and make sure it's valid if so
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
