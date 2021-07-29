import pandas as pd
from os import path
import yaml

from qebil.log import logger
from qebil.normalize import add_emp_info


_THIS_DIR, _THIS_FILENAME = path.split(__file__)


_QEBIL_KNOWN_SAMPLE_TYPES = path.join(
    _THIS_DIR, "..", "support_files", "known_sample_types.yaml"
)


_QIITA_RESTRICTED_TERMS = path.join(
    _THIS_DIR, "..", "support_files", "reserved_words.yaml"
)


_SAMPLE_NAME_ALIASES = [
    "sample id",
    "sample_id",
    "sampleid",
    "sample-id",
    "sample name",
    "#SampleID",
    "#sampleid",
]


NULL_DICT = {
    "n/a": "not applicable",
    "na": "not applicable",
    "": "not provided",
    "null": "not provided",
    "nd": "no data",
}


def load_metadata(filename):
    """
    Lightly modified from Qiita test_util.py to allow for
    csv or tsv input

    Parameters
    ----------
    filename : string
        the path to load the metadata

    Returns
    -------
    load_df : pd.DataFrame
        dataframe of metadata that was successfully loaded

    """

    load_df = pd.DataFrame()

    if path.isfile(filename):
        tmp_df = pd.read_csv(
            filename,
            sep="\t",
            dtype=str,
            encoding="utf-8",
            infer_datetime_format=False,
            keep_default_na=False,
            index_col=False,
            comment="\t",
            low_memory=False,
            header=0,
        )

        if "sample_name" not in list(tmp_df.columns):
            # handle QIIME formatted names
            found_sn = False
            for s in _SAMPLE_NAME_ALIASES:
                if not found_sn:
                    if s in tmp_df.columns:
                        found_sn = True
                        tmp_df = tmp_df.rename(columns={s: "sample_name"})
                        logger.warning("'sample_name' mapped from: " + s)

        if "sample_name" not in list(tmp_df.columns):
            # try again as csv
            tmp_df = pd.read_csv(
                filename,
                sep=",",
                dtype=str,
                encoding="utf-8",
                infer_datetime_format=False,
                keep_default_na=False,
                index_col=False,
                comment="\t",
                low_memory=False,
                header=0,
            )

            if "sample_name" not in tmp_df.columns:
                # handle QIIME formatted names
                found_sn = False
                for s in _SAMPLE_NAME_ALIASES:
                    if not found_sn:
                        if s in tmp_df.columns:
                            found_sn = True
                            tmp_df = tmp_df.rename(columns={s: "sample_name"})
                            logger.warning("'sample_name' mapped from: " + s)

        if len(tmp_df) > 0 and "sample_name" in list(tmp_df.columns):
            # remove newlines and tabs from fields
            tmp_df.replace(
                to_replace="[\t\n\r\x0b\x0c]+",
                value="",
                regex=True,
                inplace=True,
            )
            # removing columns with empty values
            tmp_df.dropna(axis="columns", how="all", inplace=True)

            # set sample_name as index
            load_df = tmp_df.set_index("sample_name")

            logger.info(
                "Loaded " + filename + " with columns " + str(load_df.columns)
            )
        else:
            logger.warning(
                "Could not find sample_name in "
                + filename
                + ", check format and try again."
                + "valid sample_name equivalents: sample_name,sample name,"
                + "sample id,sampleid,sample-id,#SampleID,#sampleid"
            )
    else:
        logger.warning(
            "Could not load file "
            + filename
            + ", check path and format and try again."
        )

    return load_df


def format_prep_type(row, sample_name, seq_methods=[]):
    """Maps EBI library strategies to Qiita prep types

    Method to map  EBI library strategies to Qiita prep types. For amplicon
    and 'other' types, also looks for target_gene to set prep type. Any
    ambigious prep type is labeled as such and this dictionary can be updated
    as new prep types are supported. Also attempts to account for common typos
    and format changes for target_gene.

    Parameters
    ----------
    expt_df : pd.DataFrame
        The dataframe of experiment data from EBI to be normalized
    row: index value
        The specific line of the dataframe to look up

    Returns
    -------
    library_strat_to_qebil_dict[ebi_strat] OR
    amplicon_dict[tg] OR 'AMBIGUOUS': string
        Qiita-normalized prep type term
    """

    amplicon_list = ["AMPLICON", "OTHER"]
    amplicon_dict = {
        "16s rna": "16S",
        "16s rrna": "16S",
        "16s": "16S",
        "16s rdna": "16S",
        "16s rna": "16S",
        "16s dna": "16S",
        "16s ": "16S",
        "its": "ITS",
        "its1": "ITS",
        "its2": "ITS",
        "its ": "ITS",
        "18s rrna": "18S",
        "18s rdna": "18S",
        "18s": "18S",
        "18s rna": "18S",
        "18s dna": "18S",
        "18s ": "18S",
    }
    library_strat_to_qebil_dict = {
        "POOLCLONE": "AMBIGUOUS",
        "CLONE": "Genome_Isolate",
        "CLONEEND": "AMBIGUOUS",
        "WGS": "Metagenomic",
        "WGA": "Metagenomic",
        "WCS": "Genome_Isolate",
        "WXS": "Metagenomic",
        "ChIP-Seq": "Metagenomic",
        "RNA-Seq": "Metatranscriptomic",
        "MRE-Seq": "AMBIGUOUS",
        "MeDIP-Seq": "AMBIGUOUS",
        "MBD-Seq": "AMBIGUOUS",
        "MNase-Seq": "AMBIGUOUS",
        "DNase-Hypersensitivity": "AMBIGUOUS",
        "Bisulfite-Seq": "Metagenomic",
        "EST": "AMBIGUOUS",
        "FL-cDNA": "AMBIGUOUS",
        "miRNA-Seq": "Metatranscriptomic",
        "ncRNA-Seq": "Metatranscriptomic",
        "FINISHING": "AMBIGUOUS",
        "CTS": "AMBIGUOUS",
        "Tn-Seq": "AMBIGUOUS",
        "VALIDATION": "AMBIGUOUS",
        "FAIRE-seq": "AMBIGUOUS",
        "SELEX": "AMBIGUOUS",
        "RIP-Seq": "AMBIGUOUS",
        "ChIA-PET": "AMBIGUOUS",
        "RAD-Seq": "AMBIGUOUS",
    }
    ebi_strat = row["library_strategy"]
    logger.info(
        "EBI strategy is:" + str(ebi_strat) + " for " + str(sample_name)
    )

    prep_type = "AMBIGUOUS"

    if ebi_strat not in amplicon_list:
        if ebi_strat in library_strat_to_qebil_dict.keys():
            prep_type = library_strat_to_qebil_dict[ebi_strat]
        else:
            logger.info(
                ebi_strat
                + " not found in EBI-Qiita library strategy mapping. "
                + "Setting to 'AMBIGUOUS'."
            )
    else:
        # since this is amplicon data, there should be a target gene,
        # if not check the seq_methods, then return AMBIGUOUS
        try:
            tg = row["target_gene"].lower()
            if tg in amplicon_dict.keys():
                prep_type = amplicon_dict[tg]
        except Exception:
            logger.warning(
                "target_gene not found for "
                + str(sample_name)
                + ". Checking study information."
            )
            if len(seq_methods) == 1:
                if seq_methods[0] in amplicon_dict.keys():
                    prep_type = amplicon_dict[seq_methods[0]]
                    logger.info(
                        "Found "
                        + str(seq_methods[0])
                        + " as only target gene in study"
                        + " information. Using as prep_type"
                        + " but user should check and specify"
                        + " target subfragment."
                    )
                else:
                    logger.warning(
                        "seq_methods "
                        + str(seq_methods)
                        + " found in study information."
                        + " Setting to 'AMBIGUOUS'."
                    )

    return prep_type


def set_criteria(
    strategies=[], platforms=[], selections=[], sources=[], names=[]
):
    """Parses the input parameters to create a dictionary of selection criteria to
    return for downstream use.

    Parameters
    ----------
    strategies: list
        list of library strategies
    platforms: list
        list of instrument platforms
    selections: list
        list of library selection strategies
    sources: list
        list of library sources
    names: list
        list of scientific names

    Returns
    -------
    selection_dict : dict
        dict of list of criteria with keys that match normalized
        metadata columns
    """
    selection_dict = {}

    # try:
    #    valid =

    # try:
    if len(strategies) > 0:
        selection_dict["library_strategy"] = [
            str(s).lower() for s in strategies
        ]
    if len(platforms) > 0:
        selection_dict["instrument_platform"] = [
            str(p).lower() for p in platforms
        ]
    if len(selections) > 0:
        selection_dict["library_selection"] = [
            str(s).lower() for s in selections
        ]
    if len(sources) > 0:
        selection_dict["library_source"] = [str(s).lower() for s in sources]
    if len(names) > 0:
        selection_dict["scientific_name"] = [str(n).lower() for n in names]

    return selection_dict


def check_sample_type(input_df, validation_dict={}):
    """Helper method to try to infer sample_type when missing from metadata

    Parameters
    -----------
    input_df: pd.DataFrame
        pandas DataFrame to check for sample type

    Returns:
    ---------
    output_df: pd.DataFrame
        pandas DataFrame with updated sample_type
    """

    # parse the validation dict to see what keys are available
    # this doesn't guarantee a match since this will be across
    # hosts, but it gives it the best shot

    available_sample_types = []

    for k1 in validation_dict.keys():
        k1_tmp = validation_dict[k1]
        for k2 in k1_tmp.keys():
            k2_tmp = k1_tmp[k2]
            for k3 in k2_tmp.keys():
                available_sample_types.append(k3)

    st_col_dict = {}

    for col in input_df.columns:
        unique_entries = [x.lower() for x in list(input_df[col].unique())]
        u_count = 0
        for u in unique_entries:
            if u in available_sample_types:
                u_count += 1

        if u_count > 0:
            st_col_dict[col] = u_count

    if len(st_col_dict) == 0:
        logger.warning(
            "Did not detect sample_type column. Setting to"
            + " 'unspecified sample type'"
        )
        input_df["qebil_sample_type"] = "unspecified sample type"
    else:
        max_matches = 0
        match_col = []
        for s in st_col_dict.keys():
            if st_col_dict[s] > max_matches:
                max_matches = st_col_dict[s]
                match_col = s
        logger.info(
            "Automatically selected column: " + match_col + "to copy."
        )
        input_df["qebil_sample_type"] = input_df[match_col].apply(
            lambda x: x.lower()
        )

    return input_df


def detect_merger_column(base_md, supp_md):
    """Helper method to compare the values in the columns of two dataframes
    to determine the column to use for merging

    Parameters
    -----------
    base_md: pd.DataFrame
        pandas DataFrame to be merged into
    supp_md: pd.DataFrame
        pandas DataFrame to be merged with

    Returns:
    ---------
    supp_col: str
        string of column name to use for merging
    """
    auto_col_list = []
    num_base_samples = len(base_md)
    num_supp_samples = len(supp_md)

    for col in supp_md.columns:
        num_supp_unique = supp_md[col].nunique()
        if num_supp_unique == num_supp_samples:
            if col in base_md.columns:
                base_unique = base_md[col].nunique()
                if base_unique == num_base_samples:
                    auto_col_list.append(col)
    found_col = False

    for col in auto_col_list:
        if not found_col:
            base_ids = list(base_md[col].unique())
            supp_ids = list(supp_md[col].unique())
            if all(x in base_ids for x in supp_ids):
                supp_col = col
                logger.info(
                    "Automatically selected column: "
                    + supp_col
                    + " for merging supplemental metadata"
                )
                found_col = True
    if not found_col:
        logger.warning(
            "Autodetection of merger column failed.\n"
            + " Skipping merger. You may specify the"
            + " column with --merge-column."
        )
        supp_col = ""
    return supp_col


def merge_metadata(base_md, supp_md, supp_col="sample_name"):
    """Helper method to merge two dataframes together

    Parameters
    -----------
    base_md: pd.DataFrame
        pandas DataFrame to be merged into
    supp_md: pd.DataFrame
        pandas DataFrame to be merged with
    supp_col: str
        string of column name to use for merging,
        passing 'auto' attempts to detect the column
        to use

    Returns:
    ---------
    merged_df: pd.DataFrame
        the resulting pandas DataFrame after attempting to merge
    """
    merged_df = base_md  # will return base_md if merging fails
    num_supp_samples = len(supp_md)

    if num_supp_samples > 0:
        logger.info("Adding supplemental metadata")
        supp_md.columns = [clean_column_name(col) for col in supp_md.columns]
        if supp_col == "auto":
            supp_col = detect_merger_column(base_md, supp_md)

        if supp_col in base_md.columns:
            # update columns with values from supplement
            base_md.update(supp_md)

            # remove shared columns having already updated them
            supp_df = supp_md[
                supp_md.columns[~supp_md.columns.isin(base_md.columns)]
            ]

            # kludge, adding back the column needed to merge on
            supp_df[supp_col] = supp_md[supp_col]

            merged_df = base_md.merge(supp_df, how="left", on=supp_col)

        else:
            logger.warning(
                "Omitting merge of supplemental info. Column '"
                + supp_col
                + "' not found in base metadata."
            )
    return merged_df


def qebil_format(md, null_term="not provided"):
    """Reformats the metadata to be ready for loading into Qiita

    First each column is cleaned up using the clean_column_name
    method, then blanks are populated with the provided null term.

    Parameters
    -----------
    md: pd.DataFrame
        pandas DataFrame to format

    Returns:
    ---------
    md: pd.DataFrame
        pandas DataFrame with Qiita-compliant columns and no blanks
    """
    # column cleanup
    clean_columns_dict = {}
    for col in md.columns:
        clean_columns_dict[col] = clean_column_name(col)

    md = md.rename(columns=clean_columns_dict)

    # null value cleanup
    md = md.fillna(null_term)

    # resolve duplicate columns
    # https://stackoverflow.com/questions/24685012/pandas-
    # dataframe-renaming-multiple-identically-named-columns
    cols = pd.Series(md.columns)
    dupe_cols = md.columns[md.columns.duplicated(keep=False)]
    if len(dupe_cols) > 0:
        logger.warning(
            "Duplicate columns found:"
            + str(dupe_cols)
            + " Renaming with suffixed numbers."
        )
    for dup in dupe_cols:
        cols[md.columns.get_loc(dup)] = [
            dup + "_" + str(d_idx) if d_idx != 0 else dup
            for d_idx in range(md.columns.get_loc(dup).sum())
        ]
    md.columns = cols

    # clean up blanks and null values
    for col in md.columns:
        md[col] = md[col].apply(lambda x: clean_nulls(x))

    return md


def augment_metadata(base_md, add_md_list=[], merge_col="", emp=False):
    """Combines a list of metadata files with a base metadata Dataframe

    This method attemps to load and merge a list of metadata files with
    a base metadata object, merging on the specified column, or attempting
    to automatically detect the merger column for each file if merge_col ="".
    Optionally also adds metadata fields according to the EMP standard 16S V4
    preparation information if requested.

    Parameters
    -----------
    base_md: pd.DataFrame
        pandas DataFrame to be merged into
    add_md_list: list
        list of csv, tsv, or txt files to load and merge
    merge_col: str
        string of column name to use for merging,
        passing "" attempts to detect the column
        to use automatically
    emp: bool
        whether or not to add the standard EMP 16S V4 prep information

    Returns:
    ---------
    merged_df: pd.DataFrame
        the resulting pandas DataFrame after attempting to merge
    """
    supp_df_dict = {}
    supp_col = ""
    output_df = base_md
    if len(add_md_list) > 0:
        if len(merge_col) <= 1:
            logger.warning(
                "No column name found for merging supplemental"
                + " metadata, defaulting to sample_name."
            )
            supp_col = "sample_name"
        else:
            supp_col = merge_col

        for f in add_md_list:
            supp_df = load_metadata(f)

            if len(supp_df) == 0:
                logger.warning(
                    "Loaded df was empty, check path and format "
                    + " and try again. Skipping merge."
                )
            elif supp_col not in supp_df.columns:
                logger.warning(
                    "Loaded datframe from file: "
                    + f
                    + " was missing merge column: "
                    + supp_col
                    + "Check names and try again. Skipping merge."
                )
            else:
                supp_df_dict[f] = supp_df

    # merge supplemental metadata
    for md in supp_df_dict.keys():
        supp_md = supp_df_dict[md]
        output_df = merge_metadata(base_md, supp_md, supp_col, md)

    # add EMP information to prep info
    if emp:
        output_df = add_emp_info(output_df)

    return output_df


def scrub_special_chars(input_string, custom_dict={}, sub="_"):
    """Removes special characters from a string

    This method scrubs special characters from any input string and
    replaces them with underscores for punctuation and spelled-out
    strings for other characters. Users may supply a custom dictionary
    which will be applied prior to the base rules to ensure that user
    preference for renaming is respected, e.g. "/" may be converted to
    '_or_' by the user to avoid conversion to '_per' by the default
    dictionary. Returned strings are also compatible with Qiita column
    requirements.

    Parameters
    ----------
    input_string : string
        the string to be scrubbed

    Returns
    -------
    clean_string: string
        the scrubbed string
    """
    # priority replacements for common biological terms
    priority_dict = {
        "+/+": sub + "plus_plus" + sub,
        "+/-": sub + "plus_minus" + sub,
        "-/-": sub + "minus_minus" + sub,
        ">=": sub + "greater_than_or_equal_to" + sub,
        "<=": sub + "greater_than_or_equal_to" + sub,
        "^-": sub + "to_power_negative" + sub,
        "[]": sub + "concentration" + sub,
    }
    replace_dict = {
        " ": sub,
        "-": sub,
        "(": sub + "leftparen" + sub,
        ")": sub + "rightparen" + sub,
        "/": sub + "per" + sub,
        "|": sub + "bar" + sub,
        "~": sub + "approximately" + sub,
        "`": sub,
        "@": sub + "at" + sub,
        "#": sub + "number" + sub,
        "$": sub + "dollar" + sub,
        "%": sub + "percent" + sub,
        "^": sub + "to_power" + sub,
        "&": sub + "and" + sub,
        "*": sub + "star" + sub,
        "+": sub + "plus" + sub,
        "=": sub + "equals" + sub,
        "\\": sub + "per" + sub,
        "{": sub + "leftbracket" + sub,
        "}": sub + "rightbracket" + sub,
        "[": sub + "leftbracket" + sub,
        "]": sub + "rightbracket" + sub,
        "?": sub + "question" + sub,
        "<": sub + "less_than" + sub,
        ">": sub + "greater_than" + sub,
        ",": sub + "comma" + sub,
        ".": sub + "dot" + sub,
        "_": sub,
        "__": sub,
    }
    # trim whitespace around string
    clean_string = str(input_string).strip()

    try:
        for k in custom_dict.keys():
            clean_string = clean_string.replace(k, custom_dict[k])
        for k in priority_dict.keys():
            if k != sub:
                clean_string = clean_string.replace(k, priority_dict[k])
        for k in replace_dict.keys():
            if k != sub:
                clean_string = clean_string.replace(k, replace_dict[k])
        # remove leading and trailing substitute chars
        clean_string = clean_string.strip(sub)

        # collapse multiple sub strings into one
        clean_string = sub.join(clean_string.split(sub))

    except AttributeError:
        logger.error(
            "Issue with input_string:"
            + str(input_string)
            + " of type:"
            + str(type(input_string))
            + "or keys in custom_dict object provided."
        )

    if clean_string != str(input_string):
        logger.info(
            "String "
            + input_string
            + " contained invalid "
            + " characters. Scrubbed string to:"
            + clean_string
        )

    return clean_string


def check_qebil_restricted_column(col):
    """
    Reformats a list of column ids to ensure that
    no restricted terms are used

    Parameters
    ----------
    col: string
        column name to be checked

    Returns
    -------
    rename_dict : dict
        dict to be used for renaming problem columns

    """

    # _qebil_restricted_terms file is copied from qiimp github
    with open(_QIITA_RESTRICTED_TERMS) as file:
        restricted_qebil_terms = yaml.load(file, Loader=yaml.BaseLoader)

    restricted_term_lower = [term.lower() for term in restricted_qebil_terms]

    col = col.lower()
    if col in restricted_term_lower:
        logger.warning(
            col
            + " in list of SQL/Qiita restricted terms,"
            + " prepending with 'user_' but keeping values."
        )
        col = "user_" + col
    # DEBUG: logger.info("Returning column " + col)
    return col


def enforce_start_characters(col):
    """
    Reformats a list of column ids to ensure that
    no numbers or special characters are present
    at the start

    Parameters
    ----------
    col: string
        column name to be checked

    Returns
    -------
    clean_col : string
        cleaned column name

    """
    clean_col = col
    test_char = col[0]
    if test_char == "_":
        clean_col = "qebil" + clean_col
    elif test_char.isnumeric():
        clean_col = "qebil_" + clean_col

    if clean_col != col:
        # TODO: pythonify strings?
        logger.info(
            "Column "
            + str(col)
            + " contained invalid"
            + " start character: "
            + str(test_char)
            + " Scrubbed"
            + " column name to new column: "
            + str(clean_col)
        )

    return clean_col


def clean_nulls(potential_null, supp_dict={}):
    """Simple function to replace common null value terms
    If the value is in the null_dict, it will be replaced,
    otherwise the original value is returned

     Parameters
    ----------
    potential_null : string
        value to try to swap for a valid null term

    Returns
    -------
    potential_null : string
        valid null

    """
    resolved_null = potential_null

    for s in supp_dict.keys():
        NULL_DICT[s] = supp_dict[s]

    if type(potential_null) in {int, float}:
        return resolved_null
    else:
        test_null = str(resolved_null).lower()

        if test_null in NULL_DICT:
            resolved_null = NULL_DICT[test_null]

        # DEBUG: logger.info("potential null is " + str(potential_null))
        # DEBUG: logger.info("resolved null is " + str(resolved_null))

        if resolved_null != potential_null:
            logger.warning(
                "Cleaned null value "
                + potential_null
                + " to"
                + " standardardized null value: "
                + resolved_null
            )

        return resolved_null


def clean_column_name(col):
    """Helper function to call other column cleanup methods"""
    return check_qebil_restricted_column(
        enforce_start_characters(scrub_special_chars(col).lower())
    )


def subset_metadata(md, max_samples, rand_sample=False):
    """Subsamples a dataframe for the specified number of samples

    Simple helper function to return a subset of the metadata
    a random subset based on the number requested

    Parameters
    ----------
    metadata: pd.DataFrame
        pandas Dataframe of meatadata
    max_samples: int
        The number of samples to take
    rand_sample: boolean
        Whether to perform a random subsample

    Returns
    -------
    sample_df
        subsampled dataframe

    """
    if isinstance(max_samples, int):
        if isinstance(rand_sample, bool):
            if rand_sample:
                md = md.sample(max_samples)
            else:
                md = md.head(max_samples)
        else:
            raise ValueError(
                "Expected boolean, received: " + str(rand_sample)
            )
    else:
        raise ValueError("Expected int, received: " + str(rand_sample))

    return md
