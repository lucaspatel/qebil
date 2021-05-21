import click
from os import path

from . import cli, add_options
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
from qebil.output import (
    write_config_file,
    write_metadata_files,
    write_qebil_info_files,
)
from qebil.process import deplete_on_the_fly
from qebil.tools.metadata import load_metadata, set_criteria, augment_metadata
from qebil.tools.util import (
    load_project_file,
    parse_document,
    scrape_ebi_ids,
    setup_output_dir,
    detect_qiita_study,
)


def fetch_remote_studies(
    study_list,
    full_details=True,
    max_samples="",
    random_subsample=False,
    output_dir="./",
    prefix="",
    overwrite=False,
):
    """Primary handler for retrieving studies from EBI

    This method creates a dict of Study objects based on their
    EBI study IDs using the from_remote constructor of the Study
    class object

    Parameters
    -----------
    study_list: list
        list of EBI study/project IDs to retrieve
    output_dir: string
        directory where the metadata may be located
    prefix: string
        string to prepend file with when searching for local metadata
    overwrite: bool
        whether to look for local metadata befre fetching
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
        local_study = False
        if not overwrite:
            # to speed up processing, use existing metadata if available
            local_study = check_existing_metadata(p, output_dir, prefix)

        if local_study:
            local_study.populate_sample_names()
            # adding step to skip metadata retrieval if already performed
            if "ebi_metadata_retrieved" not in local_study.metadata.columns:
                local_study.populate_details(full_details)
            study_dict[p] = local_study
        else:  # nothing local, so fetch
            study_dict[p] = Study.from_remote(
                p,
                full_details=full_details,
                max_samples=max_samples,
                random_subsample=random_subsample,
            )

    return study_dict


def check_existing_metadata(proj, output_dir="./", prefix=""):
    """Method to speed up reprocessing and interuptions with local files

    Parameters
    -----------
    proj: string
        EBI study/project ID to check
    output_dir: string
        directory where the metadata may be located
    prefix: string
        string to prepend file with when searching  for local metadata

    Returns
    -----------
    tmp_study: qebil.Study or False
        Study object if present


    """
    local_md = ""
    tmp_study = False

    if prefix == "":
        ebi_prefix = proj
    else:
        ebi_prefix = prefix + "_" + proj

    qiime_suffix = ".QIIME_mapping_file"
    qiime_md_file = output_dir + ebi_prefix + qiime_suffix + ".tsv"
    ebi_suffix = ".EBI_metadata"
    ebi_md_file = output_dir + ebi_prefix + ebi_suffix + ".tsv"

    if path.isfile(qiime_md_file):
        local_md = qiime_md_file
    elif path.isfile(ebi_md_file):
        local_md = ebi_md_file

    if local_md != "":
        tmp_study = Study(load_metadata(local_md), proj)

    return tmp_study


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

    # setup output directory
    output_dir = setup_output_dir(output_dir)

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

    # create empty project dict to populate
    project_dict = {}

    # get study information
    project_dict = fetch_remote_studies(
        project_list,
        full_details=qiita,
        max_samples=max_samples,
        random_subsample=random_subsample,
        output_dir=output_dir,
        prefix=prefix,
        overwrite=overwrite,
    )

    # now load provided metadata
    # TODO: move this out into its own command?
    for f in metadata_list:
        if prefix == "":
            md = load_metadata(f)
            if "study_accession" in self.metadata.columns:
                 # assume this is from EBI
                unique_study = self.metadata["study_accession"].unique()
                if len(unique_study) == 1:
                    tmp_study = Study(md,unique_study[0])
                    tmp_study.populate_sample_names()
                    self.populate_details(qiita)
                    logger.info("Loaded metadata file: " + f
                                + "\n Found EBI ID: " + tmp_study.ebi_id)
                    project_dict[tmp_study.study_id] = tmp_study
                elif len(unique_study) == 0:
                    logger.warning("No study accession in metadata")
                else:
                    raise ValueError(
                        "Metadata file should contain only one "
                        + " unique study_accession value."
                        + " Found: "
                        + str(unique_study)
                    )
            else:
                logger.error(
                    "study_accession not in metadata file, skipping" + f
                )                

    # write out files to speed up recovery processing
    write_metadata_files(
        project_dict, output_dir, prefix, suffix, False, prep_max
    )

    for proj_id in project_dict.keys():
        proj = project_dict[proj_id]
        proj.filter_samples(select_dict)

        if prefix == "":
            proj_prefix = proj_id
        else:
            proj_prefix = prefix + "_" + proj_id

        if len(proj.metadata) == 0:
            logger.error(
                "No metadata retreived for EBI ID: "
                + proj_id
                + "Check connection and study webpage: "
                + +"https://www.ebi.ac.uk/ena/browser/view/"
                + proj_id
            )
        else:
            qiita_id = detect_qiita_study(proj.metadata)
            if qiita_id:
                logger.error("Study already present in Qiita:" + qiita_id)
            else:
                if qiita:
                    write_config_file(proj.details, output_dir + proj_prefix)
                    proj.populate_preps()
                    if emp_protocol:
                        proj.metadata = add_emp_info(proj.metadata)

                    write_qebil_info_files(
                        proj, output_dir, proj_prefix, max_prep=prep_max
                    )

                    supp_md_list = list(add_metadata_file)

                    if len(supp_md_list) > 0:
                        proj.metadata = augment_metadata(
                            proj.metadata,
                            supp_md_list,
                            merge_column,
                            emp_protocol,
                        )
                    suffix = ""

                if human_removal:
                    proj.metadata = deplete_on_the_fly(
                        proj, cpus, output_dir, keep_files
                    )
                elif download_fastq:
                    proj.metadata = fetch_fastqs(
                        proj, output_dir, correct_index
                    )

                # write out updated metadata
                # bit of a kludge until if/when this is refactored
                # but goal is to write out each loop instead of at end to avoid
                # issues if a job dies early
                write_metadata_files(
                    {proj_id: proj},
                    output_dir,
                    prefix,
                    suffix,
                    qiita,
                    prep_max,
                )
