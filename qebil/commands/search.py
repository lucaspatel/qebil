import click
import pandas as pd
import requests
from xmltodict import parse

from . import cli, add_options
from qebil.core import Study
from qebil.tools.metadata import set_criteria
from qebil.commands import _OUTPUT_OPTIONS, _SUBSET_OPTIONS
from qebil.log import get_timestamp, setup_log
from qebil.tools.util import setup_output_dir


def keyword_search(search, filter_dict={}, summarize=[]):
    """fetch the study that has the search string in the
        study title generate an txt file that has all
        these studies search in these studies

    Parameters
    ----------
    search: string
        string to query for matching EBI/ENA studies
    filter_dict: dict
        dictionary of filtering criteria to apply
    summarize: list of strings
        metadata fields to be used for summarizing

    Returns
    ----------
    results_df: pd.DataFrame
        summarized table of results from EBI search
    """

    md_fields = list(filter_dict.keys())
    for s in summarize:
        if s not in md_fields:
            md_fields.append(s)

    search_results = []
    from qebil.log import logger

    # download the search result in xml file
    url = (
        "https://www.ebi.ac.uk/ena/browser/api/xml/textsearch?"
        + "domain=sra-study&query="
        + search
    )
    logger.info("Searching EBI: " + url)
    search_response = requests.get(url)
    logger.info("Response code: " + str(search_response.status_code) + "\n Content length: " + str(len(search_response.content)))
    if len(search_response.content) == 0:
        raise Exception("No studies found for search: " + search)
    else:
        search_xml_dict = parse(search_response.content)
        study_list = []
        if type(search_xml_dict["STUDY_SET"]["STUDY"]) is list:
            for study in search_xml_dict["STUDY_SET"]["STUDY"]:
                study_id = study["IDENTIFIERS"]["PRIMARY_ID"]
                study_list.append(study_id)
        else:  # handling issues with parser for single study result
            study_id = search_xml_dict["STUDY_SET"]["STUDY"]["IDENTIFIERS"][
                "PRIMARY_ID"
            ]
            study_list.append(study_id)

        if len(summarize) == 0:
            # if the user specifies no fields to summarize by then just get IDs
            return pd.DataFrame(study_list, columns=["study_id"])
        else:
            # create studies
            for study_id in study_list:
                ebi_study = Study.from_remote(study_id, md_fields)
                ebi_study.filter_samples(filter_dict)
                sample_count = len(ebi_study.metadata)
                logger.info(
                    str(sample_count)
                    + " compatible samples found in study ID "
                    + study_id
                )
                if sample_count > 0:
                    study_summ = ebi_study.summarize(summarize)
                    search_results.append(study_summ)

        # generate the output file
        if len(search_results) == 0:
            error_string = (
                "No compatible studies found with the following criteria:\n"
                + "Search term: "
                + search
                + "\n"
            )
            for filt in filter_dict.keys():
                error_string += filt + ": " + str(filter_dict[filt])
            logger.warning(error_string)
            return pd.DataFrame()
        else:
            return pd.concat(search_results)


@cli.group()
def search():
    """Search EBI/ENA for studies of interest matching criteria"""
    pass


@search.command(name="ebi")
@click.option(
    "--query",
    multiple=True,
    default=[],
    help=("EBI search string to identify studies to download."),
)
@add_options(_OUTPUT_OPTIONS)
@click.option(
    "--prefix",
    default="",
    help=(
        "Prefix to prepend to results and log files. Defaults to timestamp."
    ),
)
@add_options(_SUBSET_OPTIONS)
@click.option(
    "--summarize",
    multiple=True,
    default=[
        "study_title",
        "scientific_name",
        "library_source",
        "library_strategy",
        "library_selection",
        "instrument_platform",
    ],
    type=click.Choice(
        [
            "center_name",
            "host_scientific_name",
            "instrument_model",
            "instrument_platform",
            "library_layout",
            "library_selection",
            "library_source",
            "library_strategy",
            "scientific_name",
            "study_title",
            "submitted_format",
            "tax_id",
        ],
        case_sensitive=False,
    ),
    help=("List of metadata to group summary by, in order."),
)
def search_ebi(
    query,
    output_dir,
    prefix,
    quiet,
    source,
    strategy,
    platform,
    selection,
    scientific_name,
    summarize,
    no_filter,
):
    """Search EBI for studies that match the criteria specified

    This method performs a search in EBI for the provided terms,
    applies the selection criteria as a filter, and then aggregates
    and summarizing the results found.

    \b
    The steps are to:
    1) set up the output directory
    2) set the level of logging requested
    3) prepare a dictionary of filter criteria supplied by the
    user or using Qiita-compatible limits by default
    4) concatenate the query strings into one term
    5) perform the search and aggregate the results
    6) write out the results

    \b
    Parameters
    ----------
    query: string
        string to query for matching EBI/ENA studies
    output_dir: string
        directory for writing out search results and logs
    prefix: string
        string to prepend to results and logs
    quiet: bool
        whether to write out log messages to a file
    source: (string1,string2,...stringN)
        tuple of library source(s) to filter for
    strategy: (string1,string2,...stringN)
        tuple of library strateg(y/ies) to filter for
    platform: (string1,string2,...stringN)
        tuple of sequencing platform(s) to filter for
    scientific_name: (string1,string2,...stringN)
        tuple of scientific name(s) to filter for
    summarize: (string1,string2,...stringN)
        tuple of metadata fields to groupby when summarizing
    no_filter: bool
        whether to omit the default Qiita-compatible filter

    \b
    Returns
    ----------
    None

    """
    # setup output directory
    output_dir = setup_output_dir(output_dir)

    if prefix == "":
        prefix = get_timestamp()

    suffix = "_EBI_query"
    setup_log(output_dir, prefix, suffix, quiet)

    # setup selection criteria
    if not no_filter:
        select_dict = set_criteria(
            strategy, platform, selection, source, scientific_name
        )
    else:
        select_dict = {}

    summary = [str(s).lower() for s in summarize]
    query = "%20".join([q.replace(" ","%20") for q in list(query)]).strip('"')
    search_result = keyword_search(query, select_dict, summary)

    search_result.to_csv(
        output_dir + prefix + "_EBI_search_results.tsv", sep="\t", index=False
    )
