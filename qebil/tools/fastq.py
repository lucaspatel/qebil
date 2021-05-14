import json
from os import path, remove
from subprocess import Popen, PIPE
from shutil import move
import gzip

from qebil.log import logger


def check_valid_fastq(fastq_file, keep=True):
    """Method to check if a fastq file is valid using fqtools

    Parameters
    ----------
    fastq_file: string
        path to the fastq(.gz) file to check


    Returns
    ----------
    bool
        whether the file is valid

    Note: this method is no longer used routinely, being replace
    by get_read_count which checks validity and gets the count.

    """
    valid = False
    if path.isfile(fastq_file):
        fqtools_args = ["fqtools", "validate", fastq_file]
        fqtools_ps = Popen(fqtools_args)
        fqtools_ps.wait()

        if fqtools_ps.returncode == 0:
            valid = True
        else:
            if not keep:
                logger.warning(fastq_file + " is invalid. Removing")
                try:
                    remove(fastq_file)
                except FileNotFoundError:
                    logger.warning("Failed to remove: " + fastq_file)
            else:
                logger.warning(
                    fastq_file + " is invalid. Set keep = False to remove."
                )
    else:
        logger.warning("Could not find fastq file: " + fastq_file)

    return valid


def check_fastq_tail(fastq_file, keep=True):
    """Helper function to check that the end of a fastq file is valid.

     Parameters
    ----------
    fastq_file: string
        path to the fastq(.gz) file to check
    keep: bool
        whether to keep the file if it is corrupt

    Returns
    ----------
    valid: bool
        whether the file is valid

    Note: this method is no longer used routinely, being replace
    by get_read_count which checks validity and gets the count.

    """
    valid = False
    if path.isfile(fastq_file):
        with gzip.open(fastq_file, "rt") as f:
            lines = f.readlines()[-4]
        f.close()

        if len(lines) != 0:
            if lines[0][0] == "@" and len(lines[1]) == len(lines[3]):
                valid = True
        else:
            if not keep:
                logger.warning(fastq_file + " is invalid. Removing")
                try:
                    remove(fastq_file)
                except FileNotFoundError:
                    logger.warning("Failed to remove: " + fastq_file)
            else:
                logger.warning(
                    fastq_file + " is invalid. Set keep = False to remove."
                )
    else:
        logger.warning("Could not find fastq file: " + fastq_file)

    return valid


def get_read_count(forward_read, reverse_read=""):
    """Small helper function to get the number of reads in a file with fqtools
    Note: can likely speed up by wrapping count and validate together with read
    out of exitcode since only 0 when valid

    Parameters
    ----------
    forward_read: string
        path to the forward read
    reverse_read: string
        Optional: path to the reverse read

    Returns
    ----------
    total_reads: int
        the number of reads in the file

    """

    fqtools_args = [
        "fqtools",
        "count",
        forward_read,
    ]
    fqtools_ps = Popen(fqtools_args, stdout=PIPE)
    res = fqtools_ps.communicate()[0]
    res = res.decode().replace("\n", "")
    if fqtools_ps.returncode == 0:
        read_1_count = res
    else:
        return "fqtool error"

    if reverse_read != "":
        fqtools_args = ["fqtools", "count", reverse_read]
        fqtools_ps = Popen(fqtools_args, stdout=PIPE)
        res = fqtools_ps.communicate()[0]
        res = res.decode().replace("\n", "")

        if fqtools_ps.returncode == 0:
            read_2_count = res

            if read_2_count != read_1_count:
                return "fqtool error. R1 reads != R2 reads"
            else:
                logger.info("fqtools finished.")
                return read_1_count
        else:
            return "fqtool error"
    else:
        logger.info("fqtools finished.")
        return read_1_count


def get_read_count_fastp(forward_read, reverse_read=""):
    """Small helper function to get the number of reads in a file with fastp

    Note: now replaced by faster version for get_read_count using fqtools

    Parameters
    ----------
    forward_read: string
        path to the forward read
    reverse_read: string
        Optional: path to the reverse read

    Returns
    ----------
    total_reads: int
        the number of reads in the file


    """
    fastp_args = [
        "fastp",
        "-j",
        "temp.json",
        "-h",
        "",
        "-i",
        forward_read,
    ]

    if reverse_read != "":
        fastp_args += ["-I", reverse_read]

    # now generating a report of the current state and what would happen
    fastp_ps = Popen(fastp_args)
    fastp_ps.wait()

    if not path.isfile("temp.json"):
        logger.warning("Report generation failed.")
        return "fastp error"
    else:
        logger.info("Fastp finished.")
        with open("temp.json") as json_file:
            fastp_results = json.load(json_file)
        summary = fastp_results["summary"]
        total_reads = summary["before_filtering"]["total_reads"]
        remove("temp.json")

        return total_reads


def unpack_fastq_ftp(fastq_ftp, fastq_md5, sep=";"):
    """Unpacks the ftp and md5 field from EBI metadata

    Takes paired set of ebi-format fastq_ftp and fastq_md5
    string values and parses them into a dictionary to be used for
    downloading and using a checksum to validate the downloads

    Parameters
    ----------
    fastq_ftp: string
        string of semicolon separated ftp filepaths
    fastq_md5: string
        string of semicolon separated md5 checksums

    Returns
    ----------
    remote_dict: dict
        dict of paired ftp filepaths and md5 checksums
    """
    remote_dict = {}
    ftp_list = fastq_ftp.split(sep)
    logger.info(ftp_list)
    md5_list = fastq_md5.split(sep)

    if len(ftp_list) == 0:
        error_msg = (
            "No ftp files present. Check study details"
            + " as access may be restricted"
        )
        logger.warning(error_msg)
    else:
        read_counter = 0
        while read_counter < len(ftp_list):
            read_dict = {}
            read_dict["ftp"] = ftp_list[read_counter]
            read_dict["md5"] = md5_list[read_counter]
            read_counter += 1
            remote_dict["read_" + str(read_counter)] = read_dict

    return remote_dict


def remove_index_read_file(file_list):
    """Simple method to quickly remove and rename files when
    the user indicates that the first read (R1) is an index
    and the other read(s) are the sequences

    Parameters
    ----------
    file_list: list
        list of string filepaths to modify

    Returns
    ----------
    None

    """
    for f in file_list:
        if "_R1." in f:
            remove(f)
        elif "_R2." in f:
            move(f, f.replace("_R2.", "_R1."))
        elif "_R3." in f:
            move(f, f.replace("_R3.", "_R2."))
