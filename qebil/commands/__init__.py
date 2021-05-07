from importlib import import_module

import click

from qebil import __version__


# from https://stackoverflow.com/questions/40182157/
# shared-options-and-flags-between-commands
def add_options(options):
    def _add_options(func):
        for option in reversed(options):
            func = option(func)
        return func

    return _add_options


# these are the click options that subsequently get called and
# loaded by other command modules. Most of these terms are
# specified by the permitted values in EBI:
# https://ena-docs.readthedocs.io/en/latest/submit/reads/
# webin-cli.html?highlight=permitted%20values
_SUBSET_OPTIONS = [
    click.option(
        "--source",
        multiple=True,
        default=[],
        type=click.Choice(
            [
                "GENOMIC",
                "GENOMIC SINGLE CELL",
                "TRANSCRIPTOMIC",
                "TRANSCRIPTOMIC SINGLE CELL",
                "METAGENOMIC",
                "METATRANSCRIPTOMIC",
                "SYNTHETIC",
                "VIRAL RNA",
                "OTHER",
            ],
            case_sensitive=False,
        ),
        help=("Library sources for restricting search, not case"
              + " sensitive. Default is [] (no restriction)"),
    ),
    click.option(
        "--strategy",
        multiple=True,
        default=[
            "amplicon",
            "other",
            "wgs",
            "rna-seq",
            "wcs",
            "poolclone",
            "clone",
        ],
        type=click.Choice(
            [
                "POOLCLONE",
                "CLONE",
                "CLONEEND",
                "WGS",
                "WGA",
                "WCS",
                "WXS",
                "AMPLICON",
                "ChIP-Seq",
                "RNA-Seq",
                "MRE-Seq",
                "MeDIP-Seq",
                "MBD-Seq",
                "MNase-Seq",
                "DNase-Hypersensitivity",
                "Bisulfite-Seq",
                "EST",
                "FL-cDNA",
                "miRNA-Seq",
                "ncRNA-Seq",
                "FINISHING",
                "TS",
                "Tn-Seq",
                "VALIDATION",
                "FAIRE-seq",
                "SELEX",
                "RIP-Seq",
                "ChIA-PET",
                "RAD-Seq",
                "Other",
            ],
            case_sensitive=False,
        ),
        help=("Library strategy to restrict search; not case"
              + ' sensitive. Default is: ["amplicon", "other",'
              + ' "wgs", "rna-seq", "wcs", "poolclone", "clone"]'),
    ),
    click.option(
        "--platform",
        multiple=True,
        default=["illumina", "pacbio_smrt", "oxford_nanopore"],
        type=click.Choice(
            [
                "LS454",
                "Illumina",
                "Ion_Torrent",
                "PacBio_SMRT",
                "OXFORD_NANOPORE",
            ],
            case_sensitive=False,
        ),
        help=("Instrument platform to restrict search, not case"
              + 'sensitive. Default is: ["illumina","pacbio_smrt",'
              + ' "oxford_nanopore"]'),
    ),
    click.option(
        "--selection",
        multiple=True,
        default=[
            "random",
            "pcr",
            "random pcr",
            "rt-pcr",
            "cdna",
            "cdna_randompriming",
            "inverse rrna",
            "inverse rrna selection",
            "unspecified",
            "size fractionation",
            "repeat fractionation",
            "race",
            "other",
        ],
        type=click.Choice(
            [
                "RANDOM",
                "PCR",
                "RANDOM PCR",
                "RT-PCR",
                "HMPR",
                "MF",
                "repeat fractionation",
                "size fractionation",
                "MSLL",
                "cDNA",
                "cDNA_randomPriming",
                "cDNA_oligo_dT",
                "PolyA",
                "Oligo-dT",
                "Inverse rRNA",
                "Inverse rRNA selection",
                "ChIP",
                "ChIP-Seq",
                "MNase",
                "DNase",
                "Hybrid Selection",
                "Reduced Representation",
                "Restriction Digest",
                "5-methylcytidine antibody",
                "MBD2 protein methyl-CpG binding domain",
                "CAGE",
                "RACE",
                "MDA",
                "padlock probes capture method",
                "unspecified",
                "other",
            ],
            case_sensitive=False,
        ),
        help=("Library selection method to restrict search."
              + ' Default is: ["random", "pcr", "random pcr",'
              + ' "rt-pcr", "cdna", "cdna_randompriming",'
              + ' "inverse rrna", "inverse rrna selection",'
              + ' "unspecified", "size fractionation",'
              + ' "repeat fractionation", "race", "other"]'),
    ),
    click.option(
        "--scientific-name",
        multiple=True,
        default=[],
        help=("Scientific names to restrict search."
              + " Default is [] (no restriction)"),
    ),
    click.option(
        "--no-filter",
        is_flag=True,
        help=(
            "Ignore all defaults and do not filter samples with"
            + " any default selection criteria."
            + " Default is no flag (False)."
        ),
    ),
]

_SUBSAMPLE_OPTIONS = [
    click.option(
        "--max-samples",
        default="",
        help=("Max number of samples to grab from the study."
              + " Default is '"" (all samples)"),
    ),
    click.option(
        "--random-subsample",
        is_flag=True,
        help=(
            "When sampling, randomly select subset for processing."
            + "N.B. must supply a number with --max_samples."
            + " Default is no flag (False)."),
    ),
]

# these fields are specific to metadata manipulations using the
# normalize command, but can potentially be used elsewhere
_AUGMENT_OPTIONS = [
    click.option(
        "--add-metadata-file",
        multiple=True,
        type=click.Path(
            exists=True,
            file_okay=True,
            dir_okay=False,
            writable=False,
            readable=True,
            resolve_path=True,
            allow_dash=False,
        ),
        help=(
            "Supply additional metadata file for merging. "
            + " Uses sample_name unless otherwise specified with"
            + " --merge-column. Default is [] (none)."
        ),
    ),
    click.option(
        "--merge-column",
        default="auto",
        help=(
            "Column for merging supplemental metadata with"
            + " downloaded metadata. Default is 'auto' which"
            + " attempts to automatically determine the column"
            + " for merging."
        ),
    ),
    click.option(
        "--emp-protocol",
        is_flag=True,
        help=(
            "Update prep information with EMP protocol standards "
            + " for 16S rRNA sequencing. Default is no flag (False)."),
    ),
]

_OUTPUT_OPTIONS = [
    click.option(
        "--output-dir",
        default="./",
        help=(
            "Directory for output_dir files. "
            + " Default is working directory, './'"
        ),
    ),
    click.option(
        "--quiet/--verbose",
        default=True,
        help=(
            "Sets logging level for information. If --quiet (default)"
            + " logs will only be written to log file."),
    ),
]

# these terms are used when performing on-the-fly qc and host removal
_PROCESS_OPTIONS = [
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
        help=("Whether or not to retain raw and intermediate fastq files."
              + " Default is no flag (False)."),
    ),
]

_METADATA_OPTIONS = [
    click.option(
        "--metadata-file",
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
        help=("Metadata file(s) to use for processing specific samples."),
    ),
    click.option(
        "--prep-max",
        type=int,
        default=250,
        help=("Max number of samples per prep info file. Default is 250"),
    ),
    click.option(
        "--prefix",
        default="",
        help=("Prefix to prepend to results and log files."),
    ),
]

"""
# from redbiom repo on github
def _terribly_handle_brokenpipeerror():
    # based off http://stackoverflow.com/a/34299346
    import os
    import sys

    sys.stdout = os.fdopen(1, "w")
"""


@click.group()
@click.version_option(version=__version__)
@click.pass_context
def cli(ctx):
    pass
    # ctx.call_on_close(_terribly_handle_brokenpipeerror)


# modules to load for calling via cli
import_module("qebil.commands.search")
import_module("qebil.commands.fetch")
import_module("qebil.commands.metadata")
