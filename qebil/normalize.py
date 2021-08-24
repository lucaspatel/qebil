import collections
import pandas as pd
from pathlib import Path
import numpy as np
import subprocess
import yaml

from qebil.log import logger

_VALID_PREP_TYPES = [
    "16S",
    "18S",
    "ITS",
    "Genome_Isolate",
    "Metagenomic",
    "Metatranscriptomic",
]


def qiimp_parser(filename):
    """Parses Qiimp-format yaml or xlsx files to obtain a dictionary of
    rules for metadata normalization

    This method will accept an input xlsx from Qiimp (qiita.ucsd.edu/qiimp)
    with a yaml formatted set of rules in cell A1 of the Sheet metadata_schema,
    or a yaml/yml file with the same format which may be created by extracting
    the rules from such a Qiimp template into a new yaml/yml file. Throws
    errors if unexpected .xlsx formatted file or corrupted yaml is found.

    N.B. this method can be used to parse any yaml/yml file into a Python dict.

    Parameters
    ----------
    filename : string
        The filename of the Qiimp .xlsx or yml to be parsed

    Returns
    -------
    parse_yml
        returns the yaml file parsed into a dict object
    """
    ext = Path(filename).suffix
    parsed_yml = {}
    if ext == ".xlsx":
        try:
            temp_yml = pd.read_excel(
                filename, sheet_name="metadata_schema", header=None
            )  # assume QIIMP-style excel
            parsed_yml = yaml.load(temp_yml.at[0, 0], Loader=yaml.FullLoader)
        except Exception:
            logger.warning(
                "Invalid .xlsx file. Please ensure file is from QIIMP or "
                + "contains a compliant yaml "
                + "in cell A1 of the sheet labelled 'metadata_schema'."
            )
    elif ext == ".yml" or ext == ".yaml":
        try:
            with open(filename) as file:
                parsed_yml = yaml.load(file, Loader=yaml.FullLoader)
        except Exception:
            logger.warning(
                "Could not load yaml from " + filename + "Contents:\n"
            )
            logger.warning(subprocess.run(["cat", filename]))
    else:
        logger.warning("Invalid file extension for yaml parsing: " + str(ext))

    valid_format = True

    if len(parsed_yml) == 0:
        valid_format = False
        logger.warning(
            " Skipping validator load for "
            + str(filename)
            + " The file contains no yaml data."
            + " Please check contents and try again."
        )
    else:
        for k in parsed_yml.keys():
            try:
                key_list = parsed_yml[k].keys()
                logger.info(
                    "Keys loaded for "
                    + str(k)
                    + ": "
                    + str(key_list)
                    + " in file"
                    + filename
                )
            except Exception:
                valid_format = False
                logger.warning(
                    "Invalid format. No keys for entry: "
                    + str(k)
                    + " Skipping validator load for "
                    + str(filename)
                    + " Please check contents and try again."
                )

    if not valid_format:
        parsed_yml = {}

    return parsed_yml


def apply_validation(prevalid_df, validator_yaml):
    msg = ""
    for k in validator_yaml.keys():
        mapped = False
        if k not in prevalid_df.columns:
            if k.split("_")[0] != "qebil":
                qebil_k = "qebil_" + k
            else:
                qebil_k = k
            msg = (
                msg
                + k
                + " not found in metadata. Attempting to generate "
                + qebil_k
                + " from available sources.\n"
            )

            if "sources" in validator_yaml[k].keys():
                source_list = validator_yaml[k]["sources"]
                for src in source_list.keys():
                    if src in prevalid_df.columns and not mapped:
                        msg = msg + "Mapping from " + src + " to " + qebil_k
                        map_dict = validator_yaml[k]["sources"][src]
                        if not isinstance(map_dict, collections.Mapping):
                            msg = msg + " by direct copy \n"
                            prevalid_df[qebil_k] = prevalid_df[src]
                            if len(prevalid_df[qebil_k].unique()) == 1:
                                if prevalid_df[qebil_k][0] != "not provided":
                                    mapped = True
                                else:
                                    msg = (
                                        msg
                                        + " Mapping applied but all "
                                        + " values are still 'not "
                                        + " provided'. Trying additional"
                                        + " source fields if available."
                                    )
                            else:
                                mapped = True
                        elif "mapping" in map_dict.keys():
                            mapping_dict = map_dict["mapping"]
                            prevalid_df[qebil_k] = prevalid_df[src].apply(
                                lambda x: mapping_dict[x]
                                if x in mapping_dict.keys()
                                else x
                            )
                            msg = (
                                msg
                                + " using mapping: "
                                + str(mapping_dict)
                                + "\n"
                            )
                            if len(prevalid_df[qebil_k].unique()) == 1:
                                if prevalid_df[qebil_k][0] != "not provided":
                                    mapped = True
                                else:
                                    msg = (
                                        msg
                                        + " Mapping applied but all "
                                        + " values are still 'not "
                                        + " provided'. Trying additional"
                                        + " source fields if available."
                                    )
                            else:
                                mapped = True
                        else:
                            msg = (
                                msg
                                + " but dictionary missing keyword: "
                                + "mapping. "
                                + "skipping attempt to map with source: "
                                + src
                            )

            else:  # if there isn't a source then revert to default
                if not mapped:
                    try:
                        prevalid_df[qebil_k] = validator_yaml[k]["default"]
                        msg = (
                            msg
                            + "Setting "
                            + qebil_k
                            + " to "
                            + validator_yaml[k]["default"]
                            + "\n"
                        )
                    except Exception:
                        prevalid_df[qebil_k] = "not provided"
                        msg = (
                            msg
                            + k
                            + " has no default in yaml template."
                            + " Encoding as 'not provided'\n"
                        )
        else:
            # construct rules
            uniq = prevalid_df[k].unique()
            allowed_list = []
            min_value = ""
            max_value = ""
            min_value_excl = ""
            max_value_excl = ""
            if "anyof" in validator_yaml[k].keys():
                anyof_list = validator_yaml[k]["anyof"]
                for r in anyof_list:
                    if r["type"] == "string":
                        for a in r["allowed"]:
                            allowed_list.append(a)
                    elif r["type"] == "number":
                        if "min" in r.keys():
                            min_value = r["min"]
                        if "max" in r.keys():
                            max_value = r["max"]
                        if "min_exclusive" in r.keys():
                            min_value = r["min_exclusive"]
                        if "max_exclusive" in r.keys():
                            max_value = r["max_exclusive"]
            elif "type" in validator_yaml[k].keys():
                if validator_yaml[k]["type"] == "string":
                    if "allowed" in validator_yaml[k].keys():
                        allowed_list = validator_yaml[k]["allowed"]
                if (
                    validator_yaml[k]["type"] == "number"
                    or validator_yaml[k]["type"] == "integer"
                ):
                    if "min" in validator_yaml[k].keys():
                        min_value = validator_yaml[k]["min"]
                    if "max" in validator_yaml[k].keys():
                        max_value = validator_yaml[k]["max"]
                    if "min_exclusive" in validator_yaml[k].keys():
                        min_value_excl = validator_yaml[k]["min"]
                    if "max_exclusive" in validator_yaml[k].keys():
                        max_value_excl = validator_yaml[k]["max"]

            # alert user of issues
            for u in uniq.astype(str):
                if not u.isnumeric():
                    if u not in allowed_list and len(allowed_list) > 0:
                        msg = (
                            msg
                            + "Warning "
                            + u
                            + " found in column "
                            + k
                            + " but not allowed per Qiimp template."
                            + "valid values: "
                            + str(allowed_list)
                            + "\n"
                        )
                else:
                    if u not in allowed_list:  # assume it's actually a number
                        if min_value != "" and u < min_value:
                            msg = (
                                msg
                                + "Warning "
                                + u
                                + " found in column "
                                + k
                                + " but less than min value per yaml: "
                                + str(min_value)
                                + "\n"
                            )
                        if max_value != "" and u > max_value:
                            msg = (
                                msg
                                + "Warning "
                                + u
                                + " found in column "
                                + k
                                + " but more than max value per yaml: "
                                + str(max_value)
                                + "\n"
                            )
                        if min_value_excl != "" and u <= min_value_excl:
                            msg = (
                                msg
                                + "Warning "
                                + u
                                + " found in column "
                                + k
                                + " but less than min value per yaml: "
                                + str(min_value_excl)
                                + "\n"
                            )
                        if max_value_excl != "" and u >= max_value_excl:
                            msg = (
                                msg
                                + "Warning "
                                + u
                                + " found in column "
                                + k
                                + " but not allowed per yaml: "
                                + str(max_value_excl)
                                + "\n"
                            )
    return prevalid_df, msg


def normalize_lat_lon(valid_df):
    if "lat_lon" in valid_df.columns:
        valid_df["qebil_latitude"] = valid_df["lat_lon"].apply(
            lambda x: split_lat_lon(x, "lat")
        )
        valid_df["qebil_longitude"] = valid_df["lat_lon"].apply(
            lambda x: split_lat_lon(x, "long")
        )

    return valid_df


def add_emp_info(input_df):
    """
    Automatically adds/replaces information with the standard EMP protocol
    information for 16S rRNA gene sequencing.

    Parameters
    ----------
    input_df: pd.DataFrame
        dataframe to update with EMP protocol standard information

    Returns
    ----------
    output_df: pd.DataFrame
        dataframe with EMP protocol standard information
    """
    output_df = input_df
    output_df["qebil_prep_file"] = np.where(
        output_df["library_strategy"] == "AMPLICON",
        output_df["qebil_prep_file"].str.replace("AMBIGUOUS", "16S"),
        output_df["qebil_prep_file"],
    )
    output_df["target_gene"] = np.where(
        output_df["library_strategy"] == "AMPLICON",
        "16S rRNA",
        "not applicable",
    )
    output_df["target_subfragment"] = np.where(
        output_df["library_strategy"] == "AMPLICON",
        "V4",
        "not applicable",
    )
    if "library_construction_protocol" in output_df.columns:
        output_df["library_construction_protocol"] = np.where(
            output_df["library_strategy"] == "AMPLICON",
            "Illumina EMP protocol 515fbc, "
            + "806r amplification of 16SrRNA V4",
            output_df["library_construction_protocol"],
        )
    else:
        output_df["library_construction_protocol"] = np.where(
            output_df["library_strategy"] == "AMPLICON",
            "Illumina EMP protocol 515fbc, "
            + "806r amplification of 16SrRNA V4",
            "not applicable",
        )

    output_df["pcr_primers"] = np.where(
        output_df["library_strategy"] == "AMPLICON",
        "FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT",
        "not applicable",
    )

    output_df["primer"] = np.where(
        output_df["library_strategy"] == "AMPLICON",
        "GTGTGCCAGCMGCCGCGGTAA",
        "not applicable",
    )
    if "sequencing_meth" in output_df.columns:
        output_df["sequencing_meth"] = np.where(
            output_df["library_strategy"] == "AMPLICON",
            "Sequencing by synthesis",
            output_df["sequencing_meth"],
        )
    else:
        output_df["sequencing_meth"] = np.where(
            output_df["library_strategy"] == "AMPLICON",
            "Sequencing by synthesis",
            "not applicable",
        )

    return input_df


def split_lat_lon(lat_lon_string, coordinate=""):
    """Helper method to parse lat_lon into Qiita-friendly lat and long"""
    if "N" in lat_lon_string:
        lat = lat_lon_string.split("N")[0].strip()
        long = lat_lon_string.split("N")[-1].strip()
        if "W" in long:
            long = "-" + long.split("W")[0].strip()
        elif "E" in long:
            long = long.split("E")[0].strip()
        else:
            logger.info("E nor W found in " + +long + " returning unsplit.")
            lat = lat_lon_string
            long = lat_lon_string
    elif "S" in lat_lon_string:
        lat = "-" + lat_lon_string.split("S")[0].strip()
        long = lat_lon_string.split("S")[-1].strip()
        if "W" in long:
            long = "-" + long.split("W")[0].strip()
        elif "E" in long:
            long = long.split("E")[0].strip()
        else:
            logger.info("E nor W found in " + +long + " returning unsplit.")
            lat = lat_lon_string
            long = lat_lon_string
    else:
        logger.info(
            "N nor S found in:" + lat_lon_string + " returning original."
        )
        lat = lat_lon_string
        long = lat_lon_string
    if coordinate == "lat":
        return lat
    elif coordinate == "long":
        return long
    else:
        return lat_lon_string


def update_preps(input_df, prep_type):
    """
    Automatically replaces the prep type for all samples with the one
    provided

    Parameters
    ----------
    input_df: pd.DataFrame
        dataframe to update with EMP protocol standard information
    prep_type: string
        overwrite prep type to one of the valid Qiita forms:
            "16S", "18S", "ITS", "Genome_Isolate",
            "Metagenomic","Metatranscriptomic",

    Returns
    ----------
    output_df: pd.DataFrame
        dataframe with EMP protocol standard information
    """
    output_df = input_df.copy()
    replace_preps = _VALID_PREP_TYPES + ["AMBIGUOUS"]
    if prep_type not in _VALID_PREP_TYPES:
        logger.warning(
            "forced prep type: "
            + prep_type
            + " not in valid prep type list: "
            + str(_VALID_PREP_TYPES)
        )
    elif "qebil_prep_file" not in output_df.columns:
        logger.warning(
            "No qebil_prep_file column found, skipping prep type update"
        )
    else:
        for r in replace_preps:
            output_df["qebil_prep_file"] = output_df[
                "qebil_prep_file"
            ].str.replace(r, prep_type)

    return output_df
