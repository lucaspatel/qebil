import click
from . import cli, add_options
from os import makedirs, path

from qebil.commands import (
    _OUTPUT_OPTIONS,
    _METADATA_OPTIONS,
    _PROCESS_OPTIONS,
    _SUBSET_OPTIONS,
    _SUBSAMPLE_OPTIONS,
    _AUGMENT_OPTIONS,
)
from qebil.core import Study
from qebil.fetch import fetch_fastqs
from qebil.log import setup_log
from qebil.normalize import add_emp_info
from qebil.output import write_config_files, write_metadata_files
from qebil.process import deplete_on_the_fly

from qebil.tools.metadata import load_metadata, set_criteria, augment_metadata
from qebil.tools.util import load_project_file, parse_document, scrape_ebi_ids


def fetch_remote_studies(
    study_list, full_details=True, max_samples="", random_subsample=False
):
    """Primary handler for retrieving studies from EBI

    This method creates a dict of Study objects based on their
    EBI study IDs using the from_remote constructor of the Study
    class object

    Parameters
    -----------
    study_list: list
        list of EBI study/project IDs to retrieve
    full_details: bool
        whether to retrieve full study metadata or just the basics
    max_samples: int
        number of samples to subsample from the total
    random_subsample: bool
        whether to randomly take the subsample

    Returns
    -----------
    study_dict: dict
        dict of Study objects with study IDs as keys

    """
        
    study_dict = {}
    if max_samples != "":
        try:
            max_samples = int(max_samples)
        except ValueError:
            raise ValueError(
                "Max samples: " + str(max_samples) + " is not type int."
            )
    for p in study_list:
        study_dict[p] = Study.from_remote(
            p,
            full_details=full_details,
            max_samples=max_samples,
            random_subsample=random_subsample,
        )

    return study_dict


def check_existing_metadata(proj_list, output_dir="./", prefix=""):
    """Method to speed up reprocessing and interuptions with local files

    Parameters
    -----------
    proj_list: [string]
        list of EBI study/project IDs to retrieve
    output_dir: string
        directory where the metadata may be located
    prefix: string
        string to prepend files and keys with

    Returns
    -----------
    study_dict: dict
        dict of Study objects with study IDs as keys


    Returns
    -----------
    local_dict, msg: tuple of dict, string
        the study dict created from the local metadata files and
        a related message. If blank, no studies found

    """
    local_dict = {}
    msg = ""

    for p in proj_list:
        if prefix == "":
            ebi_prefix = p
        else:
            ebi_prefix = prefix + "_" + p
        suffix = ".EBI_metadata"
        ebi_md_file = output_dir + ebi_prefix + suffix + ".tsv"
        if path.isfile(ebi_md_file):
            local_dict[p] = ebi_md_file
            msg = "Loading study from local metadata."

    return local_dict, msg


@cli.group()
def fetch():
    """Fetch metadata from EBI/ENA for studies of
    interest matching criteria"""
    pass


_project_options = [
    click.option(
        "--ebi-id",
        multiple=True,
        default=[],
        help=("EBI/ENA project or study accession(s) to retrieve"),
    ),
    click.option(
        "--download-fastq",
        is_flag=True,
        help=("Whether to download the associated fastq files."),
    ),
    click.option(
        "--human-removal",
        is_flag=True,
        help=(
            "On-the-fly human read removal. See documentation for settings"
            + " and information for other organisms."
        ),
    ),
    click.option(
        "--qiita/--raw",
        default=True,
        help=(
            "Whether to format the output metadata as sample and prep info."
        ),
    ),
    click.option(
        "--overwrite",
        is_flag=True,
        help=("Do not check for existing metadata, overwrite results."),
    ),
    click.option(
        "--correct-index",
        is_flag=True,
        help=(
            "Fix the common issue of three reads by removing the index file."
        ),
    ),
]

_batch_options = [
    click.option(
        "--project-file",
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
        help=(
            "File with list of EBI/ENA projects or study accession(s)"
            + " to retrieve, one per line."
        ),
    ),
    click.option(
        "--publication",
        multiple=True,
        default=[],
        help=(
            "Publications (pdf, url, or plain text) containing list of"
            + " NCBI/SRA or EBI/ENA projects or study accession(s)"
            + " to retrieve."
        ),
    ),
]


@fetch.command(name="project")
@add_options(_project_options)
@add_options(_batch_options)
@add_options(_METADATA_OPTIONS)
@add_options(_OUTPUT_OPTIONS)
@add_options(_PROCESS_OPTIONS)
@add_options(_SUBSET_OPTIONS)
@add_options(_AUGMENT_OPTIONS)
@add_options(_SUBSAMPLE_OPTIONS)
def fetch_project(
    ebi_id,
    project_file,
    metadata_file,
    publication,
    output_dir,
    prefix,
    quiet,
    qiita,
    prep_max,
    download_fastq,
    human_removal,
    cpus,
    keep_files,
    source,
    strategy,
    platform,
    selection,
    scientific_name,
    max_samples,
    random_subsample,
    no_filter,
    add_metadata_file,
    merge_column,
    emp_protocol,
    overwrite,
    correct_index,
):
    """Retrieve metadata for samples from studies specified

    This method retrieves metadata (and data) for a study in EBI/ENA
    as requested either by --ebi-id or listed in the --project-file
    location, output from "qebil search ebi" or manually created.
    If a metadata file is supplied, only the samples in that file
    will be retreived rather than the whole study. Numerous options
    are provided to enable customization of the output, filtering for
    samples that match certain criteria, and augmentation with either
    additional metadata and/or with EMP protocol information.

    \b
    The steps are to:
    1) set up the output directory
    2) set the level of logging requested
    3) prepare a dictionary of filter criteria supplied by the
    user or using Qiita-compatible limits by default
    4) determine which studies and/or samples to download either
    from the list of --ebi-id entries, the --project-file entries,
    the list of --metadata-file entries, and/or the list of --publication
    entries.
    5) Download the study metadata, and details for automatic loading into
    Qiita
    6) Download the fastq files if requested

    The script is designed to be robust to restarts and performs multiple
    safety checks to ensure the files downloaded match those in the EBI/ENA
    repository.

    \b
    Several additional parameters to handle common situations are available:
    a) --max-samples and --random-subsample allow for users to probe a study
        without accessing all of the meatadata and data
    b) --add-metadata-file allows users to supply additional metadata to be
        merged with the metadata on EBI using --merge-column as the shared
        key. The term "auto" can be supplied to --merge-column to attempt to
        detect automatically
    c) --emp-protocol allows users to add standard EMP protocol metadata fields
        and valuesfor 16S V4 sequencing to preparation info files
    Untested/under development:
    d) --correct-index allows users to automatically resolve studies with three
        reads by renaming and removing the index file (presumed to be the first
        in the fastq_ftp list) but needs unittesting before going live.
    e) --human-removal performs on-the-fly (per sample) quality filtering with
        fastp and human read removal and quantification but needs unittesting
        before going live.

    \b
    Parameters
    ----------
    ebi-id: string
        string of EBI/ENA stud(y/ies) to retrieve
    project-file: string
        path to file with project IDs either one per line
        or unique IDs in the 'study_id' column as produced
        from a qebil search ebi result
    metadata-file:
        tsv or csv file containing the samples that should
        be retrieved
    publication: string
        url, pdf, or txt file containing project or study IDs
    output_dir: string
        directory for writing out search results and logs
    prefix: string
        string to prepend to results and logs
    quiet: bool
        whether to write out log messages to a file
    qiita/raw: bool
        whether to output the metadata into separate
        Qiita-compatible sample and preparation info files
    prep_max: int
        the maximum number of files to add to any given prep
        info file
    download_fastq: bool
        whether to download the fastq files from the study
    human_removal: bool
        whether to perform on-the-fly human read removal
    cpus: int
        number of threads available for processing, only
        used with human_removal
    source: (string1,string2,...stringN)
        tuple of library source(s) to filter for
    strategy: (string1,string2,...stringN)
        tuple of library strateg(y/ies) to filter for
    platform: (string1,string2,...stringN)
        tuple of sequencing platform(s) to filter for
    scientific_name: (string1,string2,...stringN)
        tuple of scientific name(s) to filter for
    no_filter: bool
        whether to omit the default Qiita-compatible filter
    max_samples: int
        max number of samples to populate with metadata and/or
        download fastq files from
    random_subsample: bool
        whether to randomly subsample the metadata
    add_metadata_file: string
        path to additional metadata to merge with the info
        retrieved from EBI/ENA
    merge_column: string
        column to use for merging the metadata. defaults to
        sample_name, or pass 'auto' to detect automatically
    emp_protocol: bool
        whether to add EMP protocol metadata fields for 16S V4
        sequencing to preparation info files
    overwrite: bool
        whether to overwrite existing files and metadata rather
        than use existing downloads. Note, sequencing files are
        automatically checked for md5checksum validity during
        download
    correct_index: bool
        whether to automatically resolve studies with three reads
        by renaming and removing the index file

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

    # prepare output and logs
    suffix = ".EBI_metadata"
    setup_log(output_dir, prefix, suffix, quiet)

    from qebil.log import logger

    # setup selection criteria
    if not no_filter:
        select_dict = set_criteria(
            strategy, platform, selection, source, scientific_name
        )
    else:
        select_dict = {}

    # intialize project_list if project ID(s) supplied
    project_list = list(ebi_id)
    metadata_list = list(metadata_file)

    # import list of project from project files if they exist
    if len(project_file) > 0:
        for pf in project_file:
            if path.isfile(pf):
                loaded_proj_ids = load_project_file(pf)
                project_list += loaded_proj_ids

    # read list of publications to scrape for project ids
    if len(publication) > 0:
        for pub in publication:
            parsed_doc = parse_document(pub)
            if len(parsed_doc) == 0:
                logger.warning("Could not parse: " + str(pub))
            else:
                found_ids = scrape_ebi_ids(parsed_doc)
                if len(found_ids) == 0:
                    logger.warning(
                        "Did not locate any project or study"
                        + " IDs in "
                        + str(pub)
                    )
                else:
                    logger.info(
                        "Found id(s): " + str(found_ids) + " in " + pub
                    )
                    project_list += found_ids

    # to speed up processing, use existing metadata if available
    if not overwrite:
        local_dict, msg = check_existing_metadata(
            project_list, output_dir, prefix
        )
        if len(msg) > 0:
            logger.info(msg)

        for local_proj in local_dict.keys():
            metadata_list.append(local_dict[local_proj])
            project_list.remove(local_proj)

    # get study information
    project_dict = fetch_remote_studies(
        project_list,
        full_details=qiita,
        max_samples=max_samples,
        random_subsample=random_subsample,
    )
    # now load the local files
    for f in metadata_list:
        if prefix == "":
            tmp_study = Study(load_metadata(f))
            if tmp_study.study_id != "not provided":
                tmp_study.populate_sample_names()
                project_dict[tmp_study.study_id] = tmp_study
            else:
                logger.error(
                    "study_accession not in metadata file, skipping" + f
                )

    updated_proj_dict = {}

    # write out files to speed up recovered processing
    write_metadata_files(
        project_dict, output_dir, prefix, suffix, False, prep_max
    )

    for proj_id in project_dict.keys():
        proj = project_dict[proj_id]
        proj.filter_samples(select_dict)

        if qiita:
            proj.populate_preps()
            if emp_protocol:
                proj.metadata = add_emp_info(proj.metadata)

        if human_removal:
            proj.metadata = deplete_on_the_fly(
                proj, cpus, output_dir, keep_files
            )
        elif download_fastq:
            proj.metadata = fetch_fastqs(proj, output_dir, correct_index)

        supp_md_list = list(add_metadata_file)

        if len(supp_md_list) > 0:
            proj.metadata = augment_metadata(
                proj.metadata, supp_md_list, merge_column, emp_protocol
            )

        updated_proj_dict[proj_id] = proj

    # write out metadata as a single table
    suffix = ".EBI_metadata"
    write_metadata_files(
        updated_proj_dict, output_dir, prefix, suffix, qiita, prep_max
    )

    if qiita:
        write_config_files(updated_proj_dict, output_dir, prefix)
