from os import path, remove
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np

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
        return "fqtools error"

    if reverse_read != "":
        fqtools_args = ["fqtools", "count", reverse_read]
        fqtools_ps = Popen(fqtools_args, stdout=PIPE)
        res = fqtools_ps.communicate()[0]
        res = res.decode().replace("\n", "")

        if fqtools_ps.returncode == 0:
            read_2_count = res

            if read_2_count != read_1_count:
                logger.warning(
                    "fqtools error. R1 reads: "
                    + str(read_1_count)
                    + " != R2 reads:"
                    + str(read_2_count)
                )
                return "fqtools error"
            else:
                logger.info("fqtools finished.")
                return read_1_count
        else:
            return "fqtools error"
    else:
        logger.info("fqtools finished.")
        return read_1_count


def blast_for_type(fastq_file, db_dict={}):
    """Method to attempt to resolve prep type via blast

    This method is an attempt to handle the common case where
    the target_gene is not provided for AMPLICON or OTHER data.
    In this case, the top 100 reads in the fastq file will
    be extracted, converted to .fasta and blasted against the
    NCBI 16S, 18S, and ITS databases. The resulting evalues are
    then averaged and compared and the database with the lowest
    average evalue is presumed to be the target gene

    Parameters
    -----------
    fastq_file: string
        filepath to fastq file
    db_dict: dict
        dict of blast databases

    Returns
    ----------

    """
    if len(db_dict) == 0:  # assume we're on barnacle
        barnacle_fp = (
            "/panfs/panfs1.ucsd.edu/panscratch/qiita/qebil/databases/blast/"
        )
        db_dict = {
            "16S": barnacle_fp + "16S_ribosomal_RNA",
            "18S": barnacle_fp + "18S_fungal_sequences",
            "ITS_a": barnacle_fp + "ITS_eukaryote_sequences",
            "ITS_b": barnacle_fp + "ITS_RefSeq_Fungi",
        }

    # take the head of the fastq file to query
    fastq_head = get_fastq_head(fastq_file)
    fasta_file = fastq_to_fasta(fastq_head)

    # now query via blastn to find best match
    match = "AMBIGUOUS"
    min_eval = 1

    for db_type in db_dict.keys():
        db = db_dict[db_type]
        tmp_res = fastq_file + ".tmp_" + db_type + ".tsv"
        blastn_args = [
            "blastn",
            "-query",
            fasta_file,
            "-task",
            "blastn",
            "-db",
            db,
            "-outfmt",
            "6",
            "-max_hsps",
            "1",
            "-max_target_seqs",
            "5",
            "-num_threads",
            "4",
            "-out",
            tmp_res,
        ]
        blastn_string = " ".join(blastn_args)
        blastn_ps = Popen(blastn_args)
        blastn_ps.wait()

        if blastn_ps.returncode == 0:
            mean_eval = np.mean(
                list(pd.read_csv(tmp_res, sep="\t", header=None)[10])
            )
            logger.info(db_type + " " + str(mean_eval))
            if mean_eval < min_eval:
                match = db_type[:3]
                min_eval = mean_eval
        else:
            logger.warning(
                "Error running blastn with parameters: " + blastn_string
            )

        # clean up temp files
        remove(tmp_res)

    logger.info("blastn complete. Best match is :" + match)

    # clean up fasta and fastq_head
    remove(fastq_head)
    remove(fasta_file)

    return match


def get_read_length(fastq_file):
    """Helper method to get the average length of
    reads in fastq file"""

    length = "error"
    if path.isfile(fastq_file):
        fqtools_args = [
            "fqtools",
            "lengthtab",
            fastq_file,
        ]
        fqtools_ps = Popen(fqtools_args, stdout=PIPE)
        res = fqtools_ps.communicate()[0]
        res = res.decode().replace("\n", "")

        if fqtools_ps.returncode == 0:
            length = res.split("\t")[0]
        else:
            logger.warning(
                "failed to get count of "
                + fastq_file
                + ". Check format and filesize."
            )
    else:
        logger.warning("Could not find fastq file: " + fastq_file)

    return length


def get_fastq_head(fastq_file, head_fn="", head_size=100):
    """Helper method to get head of fastq file

    Parameters
    -----------


    Returns
    -----------


    """
    if head_fn == "":
        head_fn = fastq_file.replace(".fastq.gz", ".head.fastq.gz")

    if path.isfile(fastq_file):
        fqtools_args = [
            "fqtools",
            "head",
            "-n",
            str(head_size),
            "-o",
            head_fn.replace(".fastq.gz", ""),  # .fastq.gz auto-added
            fastq_file,
        ]
        fqtools_ps = Popen(fqtools_args)
        fqtools_ps.wait()

        if fqtools_ps.returncode != 0:
            logger.warning(
                "failed to get head of "
                + fastq_file
                + ". Check format and filesize."
            )
            head_fn = "error"
    else:
        logger.warning("Could not find fastq file: " + fastq_file)
        head_fn = "error"

    return head_fn


def fastq_to_fasta(fastq_file, fasta_filename=""):
    """Helper method to convert fastq file to fasta file"""
    if fasta_filename == "":
        fasta_filename = fastq_file.replace(".fastq.gz", ".fasta")

    if path.isfile(fastq_file):
        fqtools_args = [
            "fqtools",
            "-F",
            "F",
            "fasta",
            "-o",
            fasta_filename.replace(".fasta", ""),  # .fasta auto-added
            fastq_file,
        ]
        fqtools_ps = Popen(fqtools_args)
        fqtools_ps.wait()

        if fqtools_ps.returncode != 0:
            logger.warning(
                "failed to convert "
                + fastq_file
                + " to .fasta; check format and filesize."
            )
            fasta_filename = "error"
    else:
        logger.warning(
            "Could not find fastq file: "
            + fastq_file
            + " to convert to .fasta"
        )
        fasta_filename = "error"

    return fasta_filename
