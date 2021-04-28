import glob

from qebil.log import logger


def write_config_files(proj_dict, output_dir, prefix):
    """Helper function for handling Study objects to write config files

    Takes a dictionary of Study objects and sends the details to the
    write_config_file function to create the actual config and title files
    needed for automatic import into Qiita

    Parameters
    ----------
    proj_dict : dict
        dict of qebil.Study objects populated with details to write out
    output_dir: string
        path where the files should be written
    prefix: string
        prefix to prepend to output info files

    Returns
    -------
    None

    """

    for proj in proj_dict.keys():
        if prefix == "":
            file_prefix = proj
        else:
            # if there is a user prefix, prepend to the project ID
            file_prefix = prefix + "_" + proj

        file_prefix = output_dir + file_prefix

        # gets the xml_dict in the Study object
        study_details = proj_dict[proj].details

        if len(study_details) > 0:
            # writes out the files needed for Qiita study generation
            print(file_prefix)
            write_config_file(study_details, file_prefix)
        else:
            logger.warning("No study details found for " + str(proj))


def write_config_file(xml_dict, prefix=""):
    """Writes the files needed for Qiita loading

    This method writes out the study_config and study_title
    files needed for automatic loading in Qiita. The format
    is very specific for downstream processing.



    Parameters
    ----------
    xml_dict: dict
        dictionary of study information for parsing
    prefix: string
        prefix to add to the output files, includes path

    Returns
    -------
    None

    """

    # prepend file prefix if not empty
    if prefix != "":
        prefix += "_"

    # This block of code is specified by Qiita test configuration
    config_string = (
        "[required]\n"
        + "timeseries_type_id = 1\n"
        + "metadata_complete = True\n"
        + "mixs_compliant = True"
    )

    # for now at least, setting the PI to default to Qiita-Help
    config_string = (
        config_string
        + "\nprincipal_investigator = Qiita-EBI Import,"
        + " qiita.help@gmail.com, See study details"
    )

    config_string = config_string + "\nreprocess = False"
    parse_dict = xml_dict["STUDY_SET"]["STUDY"]
    null_val = "XXEBIXX"

    title = null_val
    alias = "\nstudy_alias = " + null_val
    abstract = "\nstudy_abstract = " + null_val
    description = "\nstudy_description = " + null_val
    study_id = null_val
    proj_id = null_val

    if "IDENTIFIERS" in parse_dict.keys():
        if "PRIMARY_ID" in parse_dict["IDENTIFIERS"].keys():
            study_id = parse_dict["IDENTIFIERS"]["PRIMARY_ID"]
        if "SECONDARY_ID" in parse_dict["IDENTIFIERS"].keys():
            proj_id = parse_dict["IDENTIFIERS"]["SECONDARY_ID"]
        if study_id == null_val and proj_id != null_val:
            study_id = proj_id
        if proj_id == null_val and study_id != null_val:
            proj_id = study_id

    # setting alias to project ID and ignoring EBI alias unless
    # no description found
    alias = alias.replace(null_val, proj_id)
    description = description.replace(null_val, study_id) + "; "

    # adding study ID to description

    desc_dict = {}
    if "DESCRIPTOR" in parse_dict.keys():
        desc_dict = parse_dict["DESCRIPTOR"]
    else:
        logger.warning(
            "No DESCRIPTOR values found. Using " + null_val + " for values."
        )

    if len(desc_dict) > 0:
        if "STUDY_ABSTRACT" in desc_dict.keys():
            abstract = abstract.replace(null_val, desc_dict["STUDY_ABSTRACT"])
        elif "ABSTRACT" in desc_dict.keys():
            abstract = abstract.replace(null_val, desc_dict["ABSTRACT"])
        else:
            logger.warning(
                "No abstract found, using " + null_val + " for abstract"
            )

        if "STUDY_DESCRIPTION" in desc_dict.keys():
            description = description.replace(
                "XXEBIXX", desc_dict["STUDY_DESCRIPTION"]
            )
        elif "DESCRIPTION" in desc_dict.keys():
            description = description.replace(
                null_val, desc_dict["DESCRIPTION"]
            )
        else:
            if "@alias" in parse_dict.keys():
                description = description.replace(
                    null_val, parse_dict["@alias"]
                )
                logger.warning(
                    "No description found, using EBI alias '"
                    + description
                    + "' for description"
                )
            else:
                logger.warning(
                    "No description found, using "
                    + null_val
                    + " for description"
                )

        if "STUDY_TITLE" in desc_dict.keys():
            title = title.replace(null_val, desc_dict["STUDY_TITLE"])
        elif "TITLE" in desc_dict.keys():
            title = title.replace(null_val, desc_dict["TITLE"])
        else:
            logger.warning("No title found, using " + null_val + " for title")

    config_string = (
        config_string
        + alias
        + description
        + abstract
        + "\nefo_ids = 1\n[optional]"
    )

    study_config_file = prefix + "study_config.txt"
    study_title_file = prefix + "study_title.txt"

    # Write out files
    c_file = open(study_config_file, "w")
    c_file.write(config_string)
    c_file.close()

    t_file = open(study_title_file, "w")
    t_file.write(title)
    t_file.close()


def write_metadata_files(
    proj_dict,
    output_dir="./",
    prefix="",
    suffix="",
    output_qiita=True,
    prep_max=250,
):
    """Helper function for writing out metadata

    Takes a dictionary of Study objects, extracts the metadata and
    writes out a basic tsv using pandas. If output_qiita is True also
    wrtie out separate sample and prep information files with the
    write_qiita_info_file method


    Parameters
    ----------
    proj_dict : dict
        dict of qebil.Study objects populated with details to write out
    output_dir: string
        path where the files should be written
    prefix: string
        prefix to prepend to output info files
    suffix: string
        suffix to append to output info files
    max_prep: int
        max number of samples to write into any prep info file

    Returns
    -------
    None

    """
    for proj_id in proj_dict.keys():
        proj = proj_dict[proj_id]
        # prepend file prefix if not empty

        if prefix != "":
            file_prefix = prefix + "_" + proj.study_id
        else:
            file_prefix = proj.study_id

        # extract the Study metadata
        md = proj.metadata

        # check that the metadata has at least one sample
        if len(md) > 0:
            # writes out the metadata as a single file
            ebi_prefix = output_dir + file_prefix
            md.to_csv(ebi_prefix + suffix + ".tsv", sep="\t", index=True)

            if output_qiita:
                # use helper to split into info files
                write_qiita_info_files(
                    proj, output_dir, file_prefix, suffix, prep_max
                )
        else:
            logger.warning("No metadata to write for study: " + proj.study_id)


def write_qiita_info_files(
    study, output_dir, file_prefix, file_suffix="", max_prep=250
):
    """Writes out the prep and sample information files

    This is the method that writes out the final metadata file(s).
    It relies on the column name "prep_file", established by calling
    populate_preps() for the Study to parse the combined sample and prep
    information metadata into separate files. It then checks to see if
    the files specified by run_prefix are present in the output_dir, and
    creates appropriately suffixed files for the following conditions:

    "MISSING": run_prefixes not located in the output_dir
    "TOOMANYREADS": run_prefixes matching more than two fastq.gz files

    This separates files that Qiita will accept as valid for matching the
    prep information file from those that it will reject. Typical causes of
    "TOOMANYREADS" are EBI studies where the index file(s) have been included,
    or where fastp and/or minimap2 have been run into the same directory without
    cleaning up the old files. The one circumstance currently not examined is the
    case where single-end data has been uploaded with an index file, which is
    interpreted as being paire-end. A check could be added, though likely elsewhere
    to detect and fix this, along with a fix for sorting out index files included
    with paired-end data, but this is not in scope for the first launch.

    Parameters
    ----------
    study : qebil.Study
        Study object with metadata to write out
    output_dir: string
        path where the files should be written
    file_prefix: string
        prefix to prepend to output info files
    file_suffix: string
        suffix to append to output info files
    max_prep: int
        max number of samples to write into any prep info file

    Returns
    -------
    None

    """
    prefix = output_dir + "/" + file_prefix
    sample_info_filename = prefix + "_sample_info" + file_suffix + ".tsv"
    prep_file_suffix = file_suffix + ".tsv"

    final_df = study.metadata
    prep_info_columns = study.prep_columns

    if max_prep > len(final_df):
        max_prep = len(final_df)

    # write sample_info
    sample_info_cols = final_df.columns[
        ~final_df.columns.isin(prep_info_columns)
    ]
    output_sample_df = final_df[sample_info_cols]

    # check for duplicates here. N.B. Need to retain previously to enable
    # download of all runs without faffing around with prep files
    output_sample_df = output_sample_df[~output_sample_df.index.duplicated()]

    output_sample_df.to_csv(
        sample_info_filename, sep="\t", index=True, index_label="sample_name"
    )

    # now handle preps
    prep_info_columns = [
        "sample_name"
    ] + prep_info_columns  # add to list for writing out prep files

    if "qiita_prep_file" not in final_df.columns:
        logger.warning(
            "No qiita_prep_file column found, skipping prep"
            + " info file(s) writing."
        )
    else:
        amplicon_min_prep_list = [
            "target_gene",
            "target_subfragment",
            "primer",
        ]
        amplicon_type_preps = ["16S", "ITS", "18S"]

        for prep_file in final_df["qiita_prep_file"].unique():
            # adding way to check for min essential target gene
            # information where needed
            prep_df = final_df[final_df["qiita_prep_file"] == prep_file]
            if (
                prep_file.split("_")[0] in amplicon_type_preps
            ):  # check to see if the prep is amplicon-style,
                # specified by list above
                for min_prep in amplicon_min_prep_list:
                    # if amplicon-style, enforce presence or
                    # null values for minimum prep info information
                    # will throw warning, but okay with current pandas
                    if min_prep not in prep_df.columns:
                        prep_df[min_prep] = "XXEBIXX"

            # now write out the prep info files
            prep_df = prep_df[
                prep_df.columns[prep_df.columns.isin(prep_info_columns)]
            ]
            prep_df = prep_df.dropna(axis=1, how="all")
            prep_df_list = [
                prep_df[i : i + max_prep]
                for i in range(0, prep_df.shape[0], max_prep)
            ]
            prep_count = 0
            for prep in prep_df_list:
                # adding loop to write out missing files
                file_status_dict = {
                    "VALID": [],
                    "MISSING": [],
                    "TOOMANYREADS": [],
                }
                for f in prep["run_prefix"]:
                    f_list = glob.glob(output_dir + f + "*fastq.gz")
                    if len(f_list) == 0:
                        file_status_dict["MISSING"].append(f)
                    elif len(f_list) > 2:
                        file_status_dict["TOOMANYREADS"].append(f)
                    else:
                        file_status_dict["VALID"].append(f)

                for k in file_status_dict.keys():
                    subset_df = prep[
                        prep["run_prefix"].isin(file_status_dict[k])
                    ]
                    if len(subset_df) > 0:
                        if k != "VALID":  # only write for errors
                            subset_df.to_csv(
                                prefix
                                + "_prep_info_"
                                + prep_file
                                + "_part"
                                + str(prep_count)
                                + "."
                                + k
                                + prep_file_suffix,
                                sep="\t",
                                index=True,
                                index_label="sample_name",
                            )
                        else:
                            subset_df.to_csv(
                                prefix
                                + "_prep_info_"
                                + prep_file
                                + "_part"
                                + str(prep_count)
                                + prep_file_suffix,
                                sep="\t",
                                index=True,
                                index_label="sample_name",
                            )
                prep_count += 1
