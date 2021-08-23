import pandas as pd
from qebil.tools.metadata import (
    format_prep_type,
    qebil_format,
    clean_column_name,
    scrub_special_chars,
    subset_metadata,
)
from qebil.log import logger
from qebil.fetch import fetch_ebi_info, fetch_ebi_metadata
from qebil.tools.util import get_ebi_ids, parse_details


_QEBIL_PREP_INFO_COLUMNS = [
    "run_prefix",
    "study_accession",
    "experiment_accession",
    "platform",
    "instrument_model",
    "library_strategy",
    "library_source",
    "library_layout",
    "library_selection",
    "ena_checklist",
    "ena_spot_count",
    "ena_base_count",
    "ena_first_public",
    "ena_last_update",
    "instrument_platform",
    "submitted_format",
    "sequencing_method",
    "target_gene",
    "target_subfragment",
    "primer",
    "run_accession",
    "qebil_ebi_import",
    "experiment_alias",
    "experiment_title",
    "experiment_title_specific",
    "library_name",
    "qebil_prep_file",
    'experiment_design_description',
]


_READ_COLUMNS = [
    "qebil_raw_reads",
    "qebil_quality_filtered_reads",
    "qebil_non_host_reads",
    "qebil_frac_reads_passing_filter",
    "qebil_frac_non_host_reads",
]

_QEBIL_COLUMNS = [
    "local_fastq_fp",
    "local_fastq_md5",
    "fastq_ftp",
    "fastq_md5",
    'fastq_bytes',
    "qebil_notes",
]


class Study:
    def __init__(self, md=pd.DataFrame(), proj_id="not provided"):
        """Basic constructor

        Creates a basic Study object with an empty DataFrame
        unless provided and sets this to be metadata. If a
        DataFrame is provided, it will look for a the study
        ID in the study_accession column, set this property,
        and retrieve the associated EBI details. This allows
        for previously downloaded metadata to be used without
        too much overhead by retrieving the otherwise missing
        details needed to write out the Qiita support files.

        Parameters
        ----------
        md : pd.DataFrame
            dataframe of samples

        Returns
        -------
        None

        """
        self.metadata = md
        self.ebi_id = proj_id
        self.prep_columns = []

    @property
    def metadata(self):
        """Gets the metadata being stored as a pandas DataFrame"""
        return self._metadata

    @metadata.setter
    def metadata(self, md):
        """Gets the metadata being stored as a pandas DataFrame"""
        if isinstance(md, pd.DataFrame):
            self._metadata = md
        else:
            raise Exception(
                "Invalid data type, expected pd.DataFrame"
                + " received "
                + str(md)
                + " of type"
                + str(type(md))
            )

    @property
    def details(self):
        """Gets the details being stored as a dict"""
        return self._details

    @details.setter
    def details(self, det):
        """Sets the details as a dict, and rejects if other type"""
        if isinstance(det, dict):
            self._details = det
        else:
            raise ValueError("Expected dict, received " + str(type(det)))

    @property
    def ebi_id(self):
        """Gets the EBI ID used to create the study"""
        return str(self._ebi_id)

    @ebi_id.setter
    def ebi_id(self, value):
        """Sets the EBI ID if a string"""
        if not isinstance(value, str):
            raise ValueError("Did not receive string, instead: " + str(value))
        else:
            self._ebi_id = value

    @property
    def study_id(self):
        """Gets the Study ID"""
        return str(self._study_id)

    @study_id.setter
    def study_id(self, value):
        """Sets the Study ID if a string"""
        if not isinstance(value, str):
            raise ValueError("Did not receive string, instead: " + str(value))
        else:
            self._study_id = value

    @property
    def proj_id(self):
        """Gets the Project ID"""
        return str(self._proj_id)

    @proj_id.setter
    def proj_id(self, value):
        """Sets the Project ID if a string"""
        if not isinstance(value, str):
            raise ValueError("Did not receive string, instead: " + str(value))
        else:
            self._proj_id = value

    @property
    def prep_columns(self):
        """Gets the columns to be used to create prep information files"""
        return self._prep_columns

    @prep_columns.setter
    def prep_columns(self, value):
        """Sets the columns to be used to create prep information files"""
        if not isinstance(value, list):
            raise ValueError("Did not receive list, instead: " + value)
        else:
            self._prep_columns = value

    @property
    def qebil_columns(self):
        """Gets the columns to be excluded from Qiita information files"""
        return self._qebil_columns

    @qebil_columns.setter
    def qebil_columns(self, value):
        """Gets the columns to be excluded from Qiita information files"""
        if not isinstance(value, list):
            raise ValueError("Did not receive list, instead: " + value)
        else:
            self._qebil_columns = value

    @property
    def cpus(self):
        """Gets the number of cpus available"""
        return self._cpus

    @cpus.setter
    def cpus(self, value):
        """Sets the number of cpus available if int, rejects if other type"""
        if isinstance(value, int):
            self._cpus = value
        else:
            raise ValueError("Expected int, received " + str(type(value)))

    @property
    def out_dir(self):
        """Gets the output directory"""
        return str(self._out_dir)

    @out_dir.setter
    def out_dir(self, value):
        """Sets the output directory if a string"""
        if not isinstance(value, str):
            raise ValueError("Did not receive string, instead: " + str(value))
        else:
            self._out_dir = value

    @classmethod
    def from_remote(
        cls,
        ebi_id,
        fields=[],
        full_details=False,
        max_samples="",
        random_subsample=False,
    ):
        """Allows creation of study from metadata in EBI/ENA repository

        Parameters
        ----------
        ebi_id: string
            ID of the EBI/ENA project number, starts with ERP, SRP, PRJ, etc.
        fields: list
            the list of metadata columns to retrieve, if empty get all
        full_details: bool
            whether to get the sample and run metadata for each sample
        max_samples: int
            max number of samples to retrieve, if empty get all
        random_subsample: bool
            whether to randomly subsample. Used with max_samples

        Returns
        -------
        tmp_study: Study
            Study object containing metadata and details

        """
        if not isinstance(ebi_id, str):
            raise ValueError(
                "Did not receive string, instead: " + str(ebi_id)
            )
        else:
            md = fetch_ebi_metadata(ebi_id, fields)
            tmp_study = cls(md, ebi_id)
            if len(md) > 0:
                # subset the study if the user requests before the
                # time-consuming per sample/run metadata retrieval
                if isinstance(max_samples, int):
                    tmp_study.metadata = subset_metadata(
                        md, max_samples, random_subsample
                    )

                tmp_study.populate_sample_names()
                tmp_study.populate_details(full_details)
            else:
                logger.warning("No study metadata retrieved for " + ebi_id)

            return tmp_study

    def populate_details(self, full_details=False):
        """Retrieves the details of a study in EBI/ENA

        Core method to retrieve basic details of an EBI/ENA study,
        and then detailed information for each sample and run
        if full_details is set to True and update these for the
        provided Study object.

        Parameters
        ----------
        full_details: bool
            whether to retrieve detailed sample/run info

        Returns
        -------
        None


        """
        ebi_accession = self.ebi_id
        ebi_xml_dict = {}

        if ebi_accession != "not provided":
            ebi_xml_dict = fetch_ebi_info(ebi_accession)
            study_accession, project_accession = get_ebi_ids(ebi_xml_dict)
        else:
            study_accession = False

        retrieved_details = False

        if not study_accession:
            logger.error("No matching study ID found for " + ebi_accession)
            self.study_id = "not provided"
        else:  # consider simplifying loop to display info by default
            if study_accession == ebi_accession:
                logger.info(
                    ebi_accession + " is study ID," + " retrieved details."
                )
                retrieved_details = True
                self.study_id = str(study_accession)
                self.proj_id = str(project_accession)
            elif project_accession == ebi_accession:
                try:
                    # consider refactor to sort specific Errors
                    # could define custom error type for Qebil
                    ebi_xml_dict = fetch_ebi_info(study_accession)
                    if len(ebi_xml_dict) > 0:
                        logger.info(
                            ebi_accession
                            + " is project ID."
                            + " Retrieved secondary ID: "
                            + study_accession
                        )
                        retrieved_details = True
                        self.study_id = str(study_accession)
                        self.proj_id = str(project_accession)
                except Exception:
                    logger.error(
                        "No matching study ID found for project "
                        + ebi_accession
                    )
            elif project_accession:
                logger.error(
                    "No matching study ID found for " + str(ebi_accession)
                )
                self.study_id = "not provided"
                self.proj_id = project_accession
            else:
                logger.error("No EBI entry found for " + str(ebi_accession))

        if retrieved_details:
            logger.info(
                "Study info retrieved for"
                + ebi_accession
                + " Study ID: "
                + self.study_id
                + " Project ID: "
                + self.proj_id
            )
            if full_details:
                self.populate_sample_info()
                self.populate_expt_info()
                self.metadata["ebi_metadata_retrieved"] = "TRUE"
        elif full_details:
            logger.warning(
                "Skipping sample and run info retrieval."
                + " No study details retrieved for "
                + ebi_accession
            )
        else:
            logger.warning("No study info retrieved for" + ebi_accession)

        logger.info(ebi_xml_dict)
        self.details = ebi_xml_dict

    def populate_sample_info(self):
        """Retrieve information for each sample including experiment/prep
        information

        This method retrieves a list of samples with the requested study
        accession from EBI or NCBI including all provided sample and prep
        metadata in the EBI ENA repository.

        """

        identifier = "secondary_sample_accession"

        for index, row in self.metadata.iterrows():
            # TODO: could be refactored for speed, though need to retrieve
            #  per sample may be a bottleneck here. Marcus provided this url:
            # https://engineering.upside.com/a-beginners-guide-to-optimizing-
            # pandas-code-for-speed-c09ef2c6a4d6?gi=fc949808a74c
            # one option may be to extract the identifiers, create a dict of
            # the retrieved metadata, convert to a DataFrame and merge with
            # the Study.metadata
            sample_accession = row[identifier]
            sample_xml_dict = fetch_ebi_info(sample_accession)

            if "SAMPLE_SET" not in sample_xml_dict.keys():
                logger.warning(
                    "No metadata found for sample named: "
                    + sample_accession
                    + " omitting."
                )
            else:
                # add sample specific information
                sample_xml_dict = sample_xml_dict["SAMPLE_SET"]["SAMPLE"]
                self.metadata.at[
                    index, "sample_title_specific"
                ] = sample_xml_dict["TITLE"]

                sn = sample_xml_dict["SAMPLE_NAME"]
                for s in sn.keys():
                    col = clean_column_name(s)
                    self.metadata.at[index, col] = sn[s]

                sa = sample_xml_dict["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]
                for s in sa:
                    col = s["TAG"]
                    col = clean_column_name(col)
                    try:
                        self.metadata.at[index, col] = s["VALUE"]
                    except Exception:
                        logger.warning(
                            "No value found for sample attribute: "
                            + col
                            + "for sample "
                            + sample_accession
                            + '. Setting to "not provided".'
                        )
                        self.metadata.at[index, col] = "not provided"

    def populate_expt_info(self):
        """Get additional prep info from EBI

        For each sample, retrieve the associated metadata, and
        clean up the column names and null values.

        Parameters
        ----------
        None

        Returns
        -------
        None


        """
        identifier = "secondary_sample_accession"

        for index, row in self.metadata.iterrows():
            sample_accession = row[identifier]
            expt_accession = row["experiment_accession"]
            expt_xml_dict = fetch_ebi_info(expt_accession)

            if "EXPERIMENT_SET" not in expt_xml_dict.keys():
                logger.warning(
                    "No experimental metadata found for sample: "
                    + sample_accession
                    + " with experiment_accession "
                    + expt_accession
                    + ", omitting."
                )
            else:
                expt_xml_dict = expt_xml_dict["EXPERIMENT_SET"]["EXPERIMENT"]
                self.metadata.at[
                    index, "experiment_title_specific"
                ] = expt_xml_dict["TITLE"]
                self.metadata.at[
                    index, "experiment_design_description"
                ] = expt_xml_dict["DESIGN"]["DESIGN_DESCRIPTION"]

                try:
                    ea = expt_xml_dict["EXPERIMENT_ATTRIBUTES"][
                        "EXPERIMENT_ATTRIBUTE"
                    ]
                    for e in ea:
                        col = e["TAG"]
                        col = clean_column_name(col)
                        self.prep_columns = self.prep_columns + [col]
                        try:
                            self.metadata.at[index, col] = e["VALUE"]
                        except Exception:
                            logger.warning(
                                "No value found for experiment "
                                + " attribute: "
                                + col
                                + '. Setting to "not provided".'
                            )
                            self.metadata.at[index, col] = "not provided"
                except Exception:
                    logger.warning(
                        "No experiment attributes found for "
                        + expt_accession
                        + " corresponding to sample "
                        + sample_accession
                    )

    def filter_samples(self, filter_dict={}):
        """Allows users to reduce the metadata to match criteria

        For each key in the filter_dict, the matching Study metadata column
        will be filtered to only retain samples with values that match.

        Parameters
        ----------
        filter_dict : dict
            dict of metadata columns as keys and acceptable terms as values

        Returns
        -------
        None


        """
        md = self.metadata
        if isinstance(filter_dict, dict):
            for crit in filter_dict.keys():
                # make sure the column is there
                if crit not in self.metadata.columns:
                    logger.warning(crit + " not in metadata. Skipping.")
                else:
                    if isinstance(filter_dict[crit], list):
                        md = md[md[crit].str.lower().isin(filter_dict[crit])]
                    else:
                        raise ValueError(
                            "Expected list, received: "
                            + str(filter_dict[crit])
                        )
        else:
            raise ValueError("Expected dict, received: " + str(filter_dict))
        self.metadata = md

    def summarize(self, groupby_cols):
        """Groupby summary of the Study

        Creates a groupby summary of the Study using the
        columns provided


        Parameters
        ----------
        groupby_cols: list
            list of columns to summarize

        Returns
        -------
        summary_df: pd.DataFrame
            a pd.DataFrame with study_id and sample counts

        """
        summary_gb = self.metadata.groupby(groupby_cols)[
            "sample_accession"
        ].count()
        summary_df = pd.DataFrame({"samples": summary_gb}).reset_index()
        summary_df["study_id"] = self.study_id
        output_columns = ["study_id", "samples"] + groupby_cols

        return summary_df[output_columns]

    def populate_sample_names(
        self, identifier="sample_accession", run_accession="run_accession"
    ):
        """Resolves common sample name user error for EBI metadata

        We need to catch common issue where identifier is identical
        for unique samples. for now assume that in this case,
        library_name will be unique, and if it isn't combine
        sample and run names and then update the Study metadata.

        Parameters
        ----------
        identifier : string
            column name that defines the sample name used by EBI/ENA
        run_accession: string
            column name that defines the run accession used by EBI/ENA

        """
        md = self.metadata

        # if the data for some reason has a run_prefix already, replace it
        # TODO: determine if safe to drop now since Qiita check elsewhere
        # if "run_prefix" in md.columns:
        #    logger.warning(
        #        "This may be a Qiita study." + " Check library names."
        #    )
        #    md = md.rename({"run_prefix": "user_run_prefix"}, axis=1)

        lib_strats = md["library_strategy"].unique()
        lib_dfs_to_combine = []

        for lib in lib_strats:
            lib_subset = md[md["library_strategy"] == lib]
            num_sample = len(lib_subset)
            unique_samples = lib_subset[identifier].nunique()
            if num_sample == unique_samples:
                lib_subset["sample_name"] = lib_subset[identifier]
                # prepend sample name for easier alignment
                lib_subset["run_prefix"] = (
                    lib_subset["sample_name"]
                    + "."
                    + lib_subset[run_accession]
                )
            elif "library_name" in md.keys():
                # users like to put their helpful info here
                unique_lib = lib_subset["library_name"].nunique()
                if num_sample == unique_lib:
                    lib_subset["sample_name"] = lib_subset[
                        "library_name"
                    ].apply(lambda x: scrub_special_chars(str(x), sub="."))
                    # prepend sample name for easier alignment
                    lib_subset["run_prefix"] = (
                        lib_subset["sample_name"]
                        + "."
                        + lib_subset[run_accession]
                    )
                else:
                    # fall back to sample + run id
                    lib_subset["sample_name"] = (
                        lib_subset[identifier]
                        + "."
                        + lib_subset[run_accession]
                    )
                    lib_subset["run_prefix"] = lib_subset["sample_name"]
            else:
                # fall back to sample + run id
                lib_subset["sample_name"] = (
                    lib_subset[identifier] + "." + lib_subset[run_accession]
                )
                # in this case sample_name and run_prefix are the same
                lib_subset["run_prefix"] = lib_subset["sample_name"]

            lib_dfs_to_combine.append(lib_subset)

        # now recombine subsets into one
        if len(lib_dfs_to_combine) == 0:
            logger.warning(
                "No DataFrames to combine. The list of library "
                + "strategies detected were: "
                + str(lib_strats)
                + " for project: "
                + str(self.ebi_id)
            )
        else:
            md = pd.concat(lib_dfs_to_combine)
            md = md.set_index("sample_name")

        self.metadata = md

    def populate_preps(self, overwrite=False):
        """Populates the metadata with the information needed for Qiita

        Stages the study with qebil_prep_file entries in the metadata to
        enable the metadata to be split into prep info files based on this
        column for automatic loading in Qiita. Also gets the Qiita-recognized
        term for the library strategy employed, corrects column labels, and
        updates the list of prep info columns to ensure that sample and prep
        metadata are separated.

        Parameters
        ----------
        identifier : string
            column name that defines the sample name used by EBI/ENA
        run_accession: string
            column name that defines the run accession used by EBI/ENA

        Returns
        -------
        None
        """

        md = qebil_format(self.metadata)
        sample_count_dict = {}
        row_list = []

        if len(md) == 0:
            logger.error("Failed to create preps, no samples in metadata")
        elif "qebil_prep_file" in md.columns:
            logger.warning(
                "qebil_prep_file found in metadata, "
                + "skipping population of preps. "
                + " Pass overwrite=True to force generation."
            )
        else:
            try:
                md["platform"] = md["instrument_platform"]
            except KeyError:
                logger.error(
                    "'instrument_platform' not found in metadata:"
                    + str(md.columns)
                )

            # adding catch for single files, to detect genome isolate
            layout = md.iloc[0]["library_layout"]
            source = md.iloc[0]["library_source"]
            name = md.iloc[0]["scientific_name"]
            strategy = md.iloc[0]["library_strategy"]

            if (
                len(md) == 1
                and layout.lower() not in ["amplicon"]
                and (
                    "meta" not in source.lower()
                    or "metagenom" not in name.lower()
                )
            ):
                logger.info(
                    "Only one sample in study "
                    + self.proj_id
                    + " with library_strategy: "
                    + strategy
                    + " and scientific_name: "
                    + name
                    + " . Setting prep type to Genome_Isolate"
                )
                md["qebil_prep_file"] = layout + "_Genome_Isolate_0"
            else:
                # frequently the prep_type is returning as ambiguous
                # so for amplicon data if only one target gene out
                # of 16S, 18S, or ITS1/2 appears then set this as the
                # prep_type and target_gene. To enable this, we'll
                # parse then tokenize the abstract and description
                seq_methods = parse_details(self.details)["seq_method"]

                for index, row in md.iterrows():
                    sample_name = index
                    prep_type = format_prep_type(row, index, seq_methods)

                    layout = row["library_layout"]
                    logger.info("Layout is: " + str(layout))
                    if prep_type not in sample_count_dict:
                        sample_count_dict[prep_type] = {
                            layout: {sample_name: 0}
                        }
                    elif layout not in sample_count_dict[prep_type]:
                        sample_count_dict[prep_type] = {
                            layout: {sample_name: 0}
                        }
                    elif (
                        sample_name
                        not in sample_count_dict[prep_type][layout]
                    ):
                        sample_count_dict[prep_type][layout][sample_name] = 0
                    else:
                        sample_count_dict[prep_type][layout][sample_name] += 1

                    # this sets the key used for splitting the files into
                    # prep_info templates. because indices can be duplicated
                    # at this stage still, need to carefully add
                    new_row = row
                    new_row["qebil_prep_file"] = (
                        layout
                        + "_"
                        + prep_type
                        + "_"
                        + str(
                            sample_count_dict[prep_type][layout][sample_name]
                        )
                    )
                    row_list.append(new_row)

                # combine the rows back together
                md = pd.DataFrame()
                for r in row_list:
                    md = md.append(r)

            for rc in _READ_COLUMNS:
                if rc not in self.metadata.columns:
                    md[rc] = "not determined"

        self.metadata = md
        self.prep_columns = (
            self.prep_columns + _QEBIL_PREP_INFO_COLUMNS + _READ_COLUMNS
        )
        self.qebil_columns = _QEBIL_COLUMNS
