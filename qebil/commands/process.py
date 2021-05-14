import click
from . import cli, add_options
from os import makedirs, path
import glob

from qebil.commands import _OUTPUT_OPTIONS, _STUDY_OPTIONS, _PROCESS_OPTIONS
from qebil.process import run_fastp
from qebil.tools.fastq import (
    check_valid_fastq,
    check_fastq_tail,
    get_read_count,
)
from qebil.tools.util import setup_output_dir

_PROCESS_OPTIONS = [
    click.option(
        "--fastq-file",
        multiple=True,
        default=[],
        type=click.Path(
            exists=True,
            file_okay=True,
            dir_okay=False,
            writable=False,
            readable=True,
            resolve_path=True,
            allow_dash=True,
        ),
        help=("fastq.gz file(s) to process."),
    ),
    click.option(
        "--input-dir",
        default="./",
        help=(
            "Directory of fastq.gz files to process. "
            + " Default is working directory unless"
            + " fastq-file is provided."
        ),
    ),
    click.option(
        "--cpus",
        type=int,
        default=4,
        help=(
            "Number of processors to use during host depletion."
            + " Default is 4."
        ),
    ),
    click.option(
        "--keep-files",
        is_flag=True,
        help=(
            "Whether or not to retain raw, intermediate, and/or invalid fastq files."
        ),
    ),
]


@cli.group()
def process():
    """Performs common functions on fastq files"""
    pass


@process.command(name="validate")
@add_options(_PROCESS_OPTIONS)
@add_options(_OUTPUT_OPTIONS)
@click.option(
    "--method",
    default="fqtools",
    type=click.Choice(
        [
            "fqtools",
            "tail",
        ],
        case_sensitive=False,
    ),
    help=("Tool to use for validating files, either fqtools or fastqc"),
)
def validate_fastq(
    fastq_file, input_dir, output_dir, prefix, quiet, cpus, keep_files, method
):
    """Validates files"""

    # output_dir directory
    if output_dir[-1] != "/":
        output_dir += "/"

    if not path.exists(output_dir):
        makedirs(output_dir)

    # prepare output and logs
    suffix = "_validate_fastq"
    setup_log(output_dir, prefix, suffix, quiet)

    from qebil.log import logger

    if len(fastq_file) > 0:
        for f in fastq_file:
            if path.isfile(f):
                file_list += f
    else:
        file_list = glob.glob(input_dir + "/*fastq.gz")

    # now validate them
    for f in file_list:
        valid = False
        if method == "fqtools":
            valid = check_valid_fastq(f, keep_files)
        elif method == "tail":
            valid = check_fastq_tail(f, keep_files)
        if valid:
            valid_files = open(prefix + "_valid_files.txt", "a")
            valid_files.write(fastq + "\n")
            valid_files.close()
        else:
            invalid_files = open(prefix + "_invalid_files.txt", "a")
            invalid_files.write(fastq + "\n")
            invalid_files.close()


@process.command(name="validate")
@add_options(_PROCESS_OPTIONS)
@add_options(_OUTPUT_OPTIONS)
def quality_filter(
    input_dir, output_dir, prefix, quiet, cpus, keep_files, min_length
):
    """Under development, likely will require changing on-the-fly assumptions"""
    # output_dir directory
    if output_dir[-1] != "/":
        output_dir += "/"

    if not path.exists(output_dir):
        makedirs(output_dir)

    # prepare output and logs
    suffix = "_validate_fastq"
    setup_log(output_dir, prefix, suffix, quiet)

    from qebil.log import logger

    if len(fastq_file) > 0:
        for f in fastq_file:
            if path.isfile(f):
                file_list += f
    else:
        file_list = glob.glob(input_dir + "/*fastq.gz")

    # now process them
    read_count_dict = {}
    for f in file_list:
        logger.info("Processing " + f)
        raw_read_count = get_read_count(f)


def host_deplete(
    input_dir,
    output_dir,
    prefix,
    quiet,
    cpus,
    keep_files,
    min_length,
    host_db,
):
    """
    Parameters
    ----------
    fastq_file_list: list
        list of raw files to quality filter and host deplete
    filter_db : string
        path to the databased to be used for filtering
    cpus : int
        number of threads to use for fastp and fastq
    output_dir : string
        the path to save the output files
    method : string
        method for host depletion
    model : string
       the model of instrument used for sequencing
    qc : boolean
        whether to run fastqc on output
    keep : boolean
        whether to keep the raw files after quality filtering,
        and fastp files after host filtering
    min_length : int
        minimum sequencing length for quality filtering; set to 45 to permit RNA-Seq data

    Returns
    -------
    None

    """
    pass
