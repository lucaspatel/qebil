import click
from . import cli, add_options
from os import path, makedirs

from qebil.log import setup_log
from qebil.tools.metadata import load_metadata
from qebil.normalize import qiimp_parser, apply_validation, normalize_lat_lon

from qebil.commands import _OUTPUT_OPTIONS, _METADATA_OPTIONS

this_dir, this_filename = path.split(__file__)
yaml_dir = path.join(this_dir, "..", "support_files", "validators")
default_mapping = path.join(
    this_dir, "..", "support_files", "default_mappings.yml"
)


@cli.group()
def metadata():
    """Normalize metadata"""
    pass


_normalize_options = [
    click.option(
        "--validator",
        help=("One or more yaml files in QIIMP format for validation."),
    ),
    click.option(
        "--qiita-standard",
        is_flag=True,
        help=("Whether to apply the Qiita standard normalization."),
    ),
]


@metadata.command(name="normalize")
@add_options(_normalize_options)
@add_options(_METADATA_OPTIONS)
@add_options(_OUTPUT_OPTIONS)
def normalize_metadata(
    output_dir,
    prefix,
    metadata_file,
    validator,
    prep_max,
    quiet,
    qiita_standard,
):
    """Applies the specified Qiimp-formatted .xlsx or .yml file
    to the provided metadata file(s) to normalize before uploading
    to Qiita

    Parameters
    ----------
    output_dir: string
        directory for writing out search results and logs
    prefix: string
        string to prepend to results and logs
    metadata-file: string
        tsv or csv file containing the samples that should
        be retrieved
    validator: string
        QIIMP-formatted .yml or .xlsx file to use for normalizing
        and validating the metadata
     prep_max: int
        the maximum number of files to add to any given prep
        info file
    quiet: bool
        whether to write out log messages to a file
    qiita_standard: bool
        whether to apply the Qiita-standard normalization for
        terms using default_mappings.yml

    \b
    Returns
    ----------
    None

    """
    # output_dir directory
    if output_dir[-1] != "/":
        output_dir += "/"

    if not path.exists(output_dir):
        makedirs(output_dir)

    suffix = "normalized.tsv"
    setup_log(output_dir, prefix, suffix, quiet)
    from qebil.log import logger

    # configure prefix
    if prefix != "" and prefix[-1] != "_":
        prefix += "_"

    # import metadata files if present for normalization
    for f in metadata_file:
        tmp_md = load_metadata(f)
        # set file name
        ext = f.split(".")[-1]
        source_dir = "/".join(f.split("/")[:-1])
        file_name = prefix + f.replace(source_dir, "").replace(ext, "")
        output_df, msg = apply_validation(tmp_md, qiimp_parser(validator))

        if qiita_standard:
            output_df, tmp_msg = apply_validation(
                output_df, qiimp_parser(default_mapping)
            )
            msg += tmp_msg
            output_df = normalize_lat_lon(output_df)

        output_df.to_csv(
            output_dir + file_name + suffix, sep="\t", index=False
        )

        if len(msg) > 0:
            logger.warning("Errors found during validation:")
            logger.warning(msg)
            if not quiet:
                valid_log_filename = (
                    output_dir + file_name + "validation_errors.txt"
                )
                errors = open(valid_log_filename, "w")
                errors.write(msg)
                errors.close()
                logger.warning(
                    "Validation errors written to " + valid_log_filename
                )
