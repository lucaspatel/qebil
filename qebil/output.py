import glob
from os import path, remove
from os.path import join

from qebil.log import logger
from qebil.tools.util import get_ebi_ids, parse_details

from functools import partial


def update_qebil_status(output_dir, prefix, msg="", overwrite=False):
    """Method to track progress of qebil

    Parameters
    ----------
    output_dir: string
        path where the files should be written
    prefix: string
        prefix to prepend to output info files
    status: string
        type of status file to write
    msg: string
        text to include in the status file

    """
    # TODO: refactor out_folder as Study property
    out_folder = partial(join, output_dir)
    status_file = out_folder(prefix + ".qebil_status")

    # Write out status file with message if provided
    if overwrite:
        write_file(status_file, msg, "w")
    else:
        write_file(status_file, msg, "a")


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

        config_prefix = output_dir + file_prefix

        # gets the xml_dict in the Study object
        study_details = proj_dict[proj].details

        if len(study_details) > 0:
            # writes out the files needed for Qiita study generation
            write_config_file(study_details, config_prefix)
        else:
            msg = "No study details found for " + str(proj)
            logger.warning(msg)
            update_qebil_status(output_dir, file_prefix, msg)


def write_config_file(xml_dict, prefix="", null_val="XXEBIXX"):
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

    study_config_file = prefix + "study_config.txt"
    study_title_file = prefix + "study_title.txt"

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

    # TODO: refactor study details to be parsed dict from get-go
    # TODO: refactor further to store and retrieve ids as part of desc_dict
    desc_dict = parse_details(xml_dict, null_val)
    study_id, proj_id = get_ebi_ids(xml_dict)

    abstract = "\nstudy_abstract = " + desc_dict["abstract"]

    # setting alias to project ID for Qiita tracking
    alias = "\nstudy_alias = " + str(proj_id) + "; "

    # starting description with study ID
    description = "\nstudy_description = " + str(study_id) + "; "

    # then adding details
    description += desc_dict["description"]

    config_string = (
        config_string
        + alias
        + description
        + abstract
        + "\nefo_ids = 1\n[optional]"
    ).replace(
        "%", "%%"
    )  # need to avoid % sign

    title = desc_dict["title"].replace("%", "%%")  # need to avoid % sign

    # Write out files
    write_file(study_config_file, config_string)
    write_file(study_title_file, title)


def write_metadata_files(
    proj_dict,
    output_dir="./",
    prefix="",
    suffix="",
    output_qiita=True,
    prep_max=250,
    fastq_prefix=".ebi",
):
    """Helper function for writing out metadata

    Takes a dictionary of Study objects, extracts the metadata and
    writes out a basic tsv using pandas. If output_qiita is True also
    wrtie out separate sample and prep information files with the
    write_qebil_info_file method


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
    output_qiita: bool
        whether to write out Qiita-related files
    prep_max: int
        max number of samples to write into any prep info file
    fastq_prefix: string
        the prefix before the .fastq.gz extension to use for checking files

    Returns
    -------
    None

    """
    for proj_id in proj_dict.keys():
        proj = proj_dict[proj_id]
        # prepend file prefix if not empty

        if prefix != "":
            file_prefix = prefix + "_" + proj_id
        else:
            file_prefix = proj_id

        # extract the Study metadata
        md = proj.metadata

        # check that the metadata has at least one sample
        if len(md) > 0:
            if output_qiita:
                # use helper to split into info files
                write_qebil_info_files(
                    proj, output_dir, file_prefix, suffix, prep_max
                )
                suffix = ".QIIME_mapping_file"

            # writes out the metadata as a single file
            proj_prefix = output_dir + file_prefix
            md.to_csv(
                proj_prefix + suffix + ".tsv",
                sep="\t",
                index=True,
                index_label="sample_name",
            )
        else:
            msg = "No metadata to write for: " + proj_id
            logger.warning(msg)
            update_qebil_status(output_dir, file_prefix, msg)


def write_qebil_info_files(
    study,
    output_dir,
    file_prefix,
    file_suffix="",
    max_prep=250,
    update_status=True,
    fastq_prefix=".ebi",
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
    or where fastp and/or minimap2 have been run into the same directory
    without cleaning up the old files. The one circumstance currently not
    examined is the case where single-end data has been uploaded with an index
    file, which is interpreted as being paire-end. A check could be added,
    though likely elsewhere to detect and fix this, along with a fix for
    sorting out index files included with paired-end data, but this is not in
    scope for the first launch.

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
    fastq_prefix: string
        the prefix before the .fastq.gz extension to use for checking files

    Returns
    -------
    None

    """
    prefix = output_dir + "/" + file_prefix
    sample_info_filename = prefix + "_sample_info" + file_suffix + ".tsv"
    prep_file_suffix = file_suffix + ".tsv"
    valid_prep = False  # flag to enable writing of .qebil_status file
    blasted_samples = []
    analytical_notes = ""

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

    if "qebil_prep_file" not in final_df.columns:
        logger.warning(
            "No qebil_prep_file column found, skipping prep"
            + " info file(s) writing."
        )
    else:
        amplicon_min_prep_list = [
            "target_gene",
            "target_subfragment",
            "primer",
        ]
        amplicon_type_preps = ["16S", "ITS", "18S"]

        for prep_file in final_df["qebil_prep_file"].unique():
            # adding way to check for min essential target gene
            # information where needed
            prep_df = final_df[final_df["qebil_prep_file"] == prep_file]
            layout = prep_file.split("_")[0]
            if (
                prep_file.split("_")[1] in amplicon_type_preps
            ):  # check to see if the prep is amplicon-style,
                # specified by list above
                for min_prep in amplicon_min_prep_list:
                    # if amplicon-style, enforce presence or
                    # null values for minimum prep info information
                    # will throw warning, but okay with current pandas
                    if min_prep not in prep_df.columns:
                        prep_df[min_prep] = "XXEBIXX"
                        if min_prep == "target_gene":
                            # now that we're blasting to determine type,
                            # add note to this effect since we can only
                            # be in this loop if the type a specific
                            # amplicon prep despite not having a target_gene
                            blasted_samples += list(prep_df["run_prefix"])

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
                valid_fp = (
                    prefix
                    + "_prep_info_"
                    + prep_file
                    + "_part"
                    + str(prep_count)
                    + prep_file_suffix
                )
                missing_fp = (
                    prefix
                    + "_prep_info_"
                    + prep_file
                    + "_part"
                    + str(prep_count)
                    + ".MISSING"
                    + prep_file_suffix
                )
                toomany_fp = (
                    prefix
                    + "_prep_info_"
                    + prep_file
                    + "_part"
                    + str(prep_count)
                    + ".TOOMANYREADS"
                    + prep_file_suffix
                )

                file_status_dict = {
                    "VALID": {"fp": valid_fp, "files": []},
                    "MISSING": {"fp": missing_fp, "files": []},
                    "TOOMANYREADS": {"fp": toomany_fp, "files": []},
                }
                for f in prep["run_prefix"]:
                    f_list = glob.glob(
                        output_dir + str(f) + "*" + fastq_prefix + ".fastq.gz"
                    )
                    if len(f_list) == 0:
                        file_status_dict["MISSING"]["files"].append(str(f))
                        logger.warning("fastq file(s) missing for " + str(f))
                    elif len(f_list) == 1:
                        if layout == "SINGLE":
                            file_status_dict["VALID"]["files"].append(str(f))
                        elif layout == "PAIRED":
                            logger.warning("fastq file missing for " + str(f))
                            file_status_dict["MISSING"]["files"].append(
                                str(f)
                            )
                            # for f in f_list:
                            # remove(f)
                        else:
                            logger.error("Layout is unexpected: " + layout)
                    elif len(f_list) == 2:
                        if layout == "SINGLE":
                            file_status_dict["TOOMANYREADS"]["files"].append(
                                str(f)
                            )
                            logger.warning(
                                "Too many reads("
                                + len(f_list)
                                + ") for "
                                + str(f)
                            )
                            # for f in f_list:
                            # remove(f)
                        elif layout == "PAIRED":
                            file_status_dict["VALID"]["files"].append(str(f))
                        else:
                            logger.error("Layout is unexpected: " + layout)
                    else:
                        file_status_dict["TOOMANYREADS"]["files"].append(
                            str(f)
                        )
                        logger.warning(
                            "Too many reads("
                            + str(len(f_list))
                            + ") for "
                            + str(f)
                            + "with layout "
                            + layout
                            + " Try running again with --correct-index"
                        )

                # see if there are valid files first
                valid_fp = file_status_dict["VALID"]["fp"]
                if path.isfile(valid_fp):
                    remove(file_status_dict["VALID"]["fp"])
                valid_df = prep[
                    prep["run_prefix"].isin(
                        file_status_dict["VALID"]["files"]
                    )
                ]
                if len(valid_df) > 0:
                    # update as complete
                    if not valid_prep and update_status:
                        update_qebil_status(
                            output_dir,
                            file_prefix,
                            msg="complete\n",
                            overwrite=True,
                        )
                        valid_prep = True

                    valid_df.to_csv(
                        valid_fp,
                        sep="\t",
                        index=True,
                        index_label="sample_name",
                    )

                for k in file_status_dict.keys():
                    if k != "VALID":
                        # add loop to cleanup previous files to avoid confusion
                        output_fp = file_status_dict[k]["fp"]
                        if path.isfile(output_fp):
                            remove(file_status_dict[k]["fp"])

                        subset_df = prep[
                            prep["run_prefix"].isin(
                                file_status_dict[k]["files"]
                            )
                        ]

                        if len(subset_df) > 0:
                            if update_status:
                                update_qebil_status(
                                    output_dir,
                                    file_prefix,
                                    "with issue: "
                                    + k
                                    + " files for "
                                    + prep_file
                                    + "\n",
                                )
                            subset_df.to_csv(
                                output_fp,
                                sep="\t",
                                index=True,
                                index_label="sample_name",
                            )
                prep_count += 1

                # add analytic notes if needed
                if len(blasted_samples) > 0:
                    analytical_notes += (
                        "The files with the following "
                        + " run_prefix were blasted to determine type"
                        + " due to missing target_gene information:"
                        + str(blasted_samples)
                    )

                # if there are analytical notes write them
                if len(analytical_notes) > 0:
                    write_file(
                        output_dir + file_prefix + "_analytical_notes.txt",
                        analytical_notes,
                        "a",
                    )


def write_file(filename, contents, mode="w"):
    """Helper method to write out text files"""

    if mode not in ["w", "a"]:
        logger.error("Write mode" + mode + " is invalid.")
    else:
        out_file = open(filename, mode)
        out_file.write(contents)
        out_file.close()
