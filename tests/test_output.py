import unittest

from qebil.output import (
    write_config_files,
    write_config_file,
    write_metadata_files,
    write_qebil_info_files,
)
from qebil.core import Study
from collections import OrderedDict
from os import path, remove
import glob
import pandas as pd
import xmltodict

from pandas.testing import assert_frame_equal

from qebil.tools.util import setup_output_dir

_THIS_DIR, _THIS_FILENAME = path.split(__file__)

_TEST_SUPPORT_DIR = path.join(_THIS_DIR, "support_files")

_TEST_OUTPUT_DIR = path.join(_THIS_DIR, "test_output/")

setup_output_dir(_TEST_OUTPUT_DIR)

test_study_id = "SRP283872"


class OutputTest(unittest.TestCase):
    """
    Tests for functions in the output module.
    """

    def setUp(self):
        self.maxDiff = None
        # clean up the directory at the start
        cleanup_list = glob.glob(_TEST_OUTPUT_DIR + "/*.EBI_metadata.tsv")
        for c in cleanup_list:
            remove(c)

    def tearDown(self):
        self.runner = None

    def test_stage_output(self):
        if not path.isdir(_TEST_OUTPUT_DIR):
            makedirs(_TEST_OUTPUT_DIR)
        else:
            cleanup_list = glob.glob(_TEST_OUTPUT_DIR + "*")
            for c in cleanup_list:
                remove(c)

    def test_write_config_files(self):
        test_config_file_contents = (
            "[required]\n"
            + "timeseries_type_id = 1\n"
            + "metadata_complete = True\n"
            + "mixs_compliant = True\n"
            + "principal_investigator = Qiita-EBI Import, "
            + "qiita.help@gmail.com, See study details\n"
            + "reprocess = False\n"
            + "study_alias = PRJNA660883\n"
            + "study_description = SRP283872; \n"
            + "study_abstract = 16S and metagenomic "
            + "sequences from individuals with acute "
            + "COVID19\n"
            + "efo_ids = 1\n"
            + "[optional]"
        )

        test_title_file_contents = (
            "human stool Metagenome from " + "acute COVID19+ subjects"
        )
        test_study_dict = {test_study_id: Study.from_remote(test_study_id)}
        write_config_files(test_study_dict, _TEST_OUTPUT_DIR, prefix="test")
        config_file = open(
            _TEST_OUTPUT_DIR + "test_" + test_study_id + "_study_config.txt",
            "r",
        )
        config_test_contents = config_file.read()
        config_file.close()
        title_file = open(
            _TEST_OUTPUT_DIR + "test_" + test_study_id + "_study_title.txt",
            "r",
        )
        title_test_contents = title_file.read()
        title_file.close()

        self.assertEqual(test_config_file_contents, config_test_contents)
        self.assertEqual(test_title_file_contents, title_test_contents)

    def test_write_config_file(self):
        """Tests the creation of the config and title files needed for
        automatic study creation in Qiita

        """
        test_config_file = open(
            _TEST_SUPPORT_DIR + "/test_study_config.txt", "r"
        )
        test_config_file_contents = test_config_file.read()
        test_config_file.close()
        test_title_file = open(
            _TEST_SUPPORT_DIR + "/test_study_title.txt", "r"
        )
        test_title_file_contents = test_title_file.read()
        test_title_file.close()

        test_xml_dict = OrderedDict(
            [
                (
                    "STUDY_SET",
                    OrderedDict(
                        [
                            (
                                "STUDY",
                                OrderedDict(
                                    [
                                        ("@accession", "SRP283872"),
                                        ("@alias", "PRJNA660883"),
                                        ("@center_name", "BioProject"),
                                        ("@broker_name", "NCBI"),
                                        (
                                            "IDENTIFIERS",
                                            OrderedDict(
                                                [
                                                    (
                                                        "PRIMARY_ID",
                                                        "SRP283872",
                                                    ),
                                                    (
                                                        "SECONDARY_ID",
                                                        "PRJNA660883",
                                                    ),
                                                    (
                                                        "EXTERNAL_ID",
                                                        OrderedDict(
                                                            [
                                                                (
                                                                    "@label",
                                                                    "primary",
                                                                ),
                                                                (
                                                                    "@namespace",
                                                                    "BioProject",
                                                                ),
                                                                (
                                                                    "#text",
                                                                    "PRJNA660883",
                                                                ),
                                                            ]
                                                        ),
                                                    ),
                                                    (
                                                        "SUBMITTER_ID",
                                                        OrderedDict(
                                                            [
                                                                (
                                                                    "@namespace",
                                                                    "BioProject",
                                                                ),
                                                                (
                                                                    "#text",
                                                                    "PRJNA660883",
                                                                ),
                                                            ]
                                                        ),
                                                    ),
                                                ]
                                            ),
                                        ),
                                        (
                                            "DESCRIPTOR",
                                            OrderedDict(
                                                [
                                                    (
                                                        "STUDY_TITLE",
                                                        "human stool Metagenome from acute COVID19+ subjects",
                                                    ),
                                                    (
                                                        "STUDY_TYPE",
                                                        OrderedDict(
                                                            [
                                                                (
                                                                    "@existing_study_type",
                                                                    "Metagenomics",
                                                                )
                                                            ]
                                                        ),
                                                    ),
                                                    (
                                                        "STUDY_ABSTRACT",
                                                        "16S and metagenomic sequences from individuals with acute COVID19",
                                                    ),
                                                    (
                                                        "CENTER_PROJECT_NAME",
                                                        "human stool",
                                                    ),
                                                ]
                                            ),
                                        ),
                                        (
                                            "STUDY_LINKS",
                                            OrderedDict(
                                                [
                                                    (
                                                        "STUDY_LINK",
                                                        [
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "XREF_LINK",
                                                                        OrderedDict(
                                                                            [
                                                                                (
                                                                                    "DB",
                                                                                    "ENA-SAMPLE",
                                                                                ),
                                                                                (
                                                                                    "ID",
                                                                                    "SRS7392524-SRS7392559",
                                                                                ),
                                                                            ]
                                                                        ),
                                                                    )
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "XREF_LINK",
                                                                        OrderedDict(
                                                                            [
                                                                                (
                                                                                    "DB",
                                                                                    "ENA-EXPERIMENT",
                                                                                ),
                                                                                (
                                                                                    "ID",
                                                                                    "SRX9152521-SRX9152556",
                                                                                ),
                                                                            ]
                                                                        ),
                                                                    )
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "XREF_LINK",
                                                                        OrderedDict(
                                                                            [
                                                                                (
                                                                                    "DB",
                                                                                    "ENA-RUN",
                                                                                ),
                                                                                (
                                                                                    "ID",
                                                                                    "SRR12672280-SRR12672315",
                                                                                ),
                                                                            ]
                                                                        ),
                                                                    )
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "XREF_LINK",
                                                                        OrderedDict(
                                                                            [
                                                                                (
                                                                                    "DB",
                                                                                    "ENA-SUBMISSION",
                                                                                ),
                                                                                (
                                                                                    "ID",
                                                                                    "SRA1127885",
                                                                                ),
                                                                            ]
                                                                        ),
                                                                    )
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "XREF_LINK",
                                                                        OrderedDict(
                                                                            [
                                                                                (
                                                                                    "DB",
                                                                                    "ENA-FASTQ-FILES",
                                                                                ),
                                                                                (
                                                                                    "ID",
                                                                                    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP283872&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes",
                                                                                ),
                                                                            ]
                                                                        ),
                                                                    )
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "XREF_LINK",
                                                                        OrderedDict(
                                                                            [
                                                                                (
                                                                                    "DB",
                                                                                    "ENA-SUBMITTED-FILES",
                                                                                ),
                                                                                (
                                                                                    "ID",
                                                                                    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP283872&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format",
                                                                                ),
                                                                            ]
                                                                        ),
                                                                    )
                                                                ]
                                                            ),
                                                        ],
                                                    )
                                                ]
                                            ),
                                        ),
                                        (
                                            "STUDY_ATTRIBUTES",
                                            OrderedDict(
                                                [
                                                    (
                                                        "STUDY_ATTRIBUTE",
                                                        [
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "ENA-SPOT-COUNT",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "166909868",
                                                                    ),
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "ENA-BASE-COUNT",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "50072960400",
                                                                    ),
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "ENA-FIRST-PUBLIC",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "2020-10-20",
                                                                    ),
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "ENA-LAST-UPDATE",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "2020-10-20",
                                                                    ),
                                                                ]
                                                            ),
                                                        ],
                                                    )
                                                ]
                                            ),
                                        ),
                                    ]
                                ),
                            )
                        ]
                    ),
                )
            ]
        )
        write_config_file(test_xml_dict, _TEST_OUTPUT_DIR + "test")

        # now read in output file
        config_test_file = open(
            _TEST_OUTPUT_DIR + "test_study_config.txt", "r"
        )
        config_test_contents = config_test_file.read()
        config_test_file.close()

        title_test_file = open(_TEST_OUTPUT_DIR + "test_study_title.txt", "r")
        title_test_contents = title_test_file.read()
        title_test_file.close()

        self.assertEqual(test_config_file_contents, config_test_contents)
        self.assertEqual(test_title_file_contents, title_test_contents)

    def test_write_metadata_files(self):
        test_study = Study.from_remote(test_study_id, full_details=True)
        test_study.populate_preps()
        test_study_dict = {test_study_id: test_study}
        write_metadata_files(
            test_study_dict,
            _TEST_OUTPUT_DIR,
            prefix="test_raw",
            suffix=".EBI_metadata",
            output_qiita=False,
        )
        write_metadata_files(
            test_study_dict, _TEST_OUTPUT_DIR, prefix="test_qiita"
        )
        test_raw_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/test_raw.EBI_metadata.tsv",
            sep="\t",
            header=0,
        )
        test_sample_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/test_sample_info.EBI_metadata.tsv",
            sep="\t",
            header=0,
        )
        test_prep_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/test_prep_info.EBI_metadata.tsv",
            sep="\t",
            header=0,
        )

        created_raw_df = pd.read_csv(
            _TEST_OUTPUT_DIR + "/test_raw_SRP283872.EBI_metadata.tsv",
            sep="\t",
            header=0,
        )
        created_sample_df = pd.read_csv(
            _TEST_OUTPUT_DIR + "/test_qiita_SRP283872_sample_info.tsv",
            sep="\t",
            header=0,
        )
        created_prep_df = pd.read_csv(
            _TEST_OUTPUT_DIR
            + "/test_qiita_SRP283872_prep_info_Metagenomic_0_part0.MISSING.tsv",
            sep="\t",
            header=0,
        )
        assert_frame_equal(test_raw_df, created_raw_df)
        assert_frame_equal(test_sample_df, created_sample_df)
        assert_frame_equal(test_prep_df, created_prep_df)

    def test_write_qebil_info_files(self):
        """Tests the creation of sample and prep information files"""
        test_study = Study.from_remote(test_study_id, full_details=True)
        test_study.populate_preps()
        write_qebil_info_files(
            test_study,
            _TEST_OUTPUT_DIR,
            "test_qiita2",
        )
        test_sample_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/test_sample_info.EBI_metadata.tsv",
            sep="\t",
            header=0,
        )
        test_prep_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/test_prep_info.EBI_metadata.tsv",
            sep="\t",
            header=0,
        )

        created_sample_df = pd.read_csv(
            _TEST_OUTPUT_DIR + "/test_qiita2_sample_info.tsv",
            sep="\t",
            header=0,
        )
        created_prep_df = pd.read_csv(
            _TEST_OUTPUT_DIR
            + "/test_qiita2_prep_info_PAIRED_Metagenomic_0_part0.MISSING.tsv",
            sep="\t",
            header=0,
        )

        assert_frame_equal(test_sample_df, created_sample_df)
        assert_frame_equal(test_prep_df, created_prep_df)

        cleanup_list = glob.glob(_TEST_OUTPUT_DIR + "/*.EBI_metadata.tsv")
        for c in cleanup_list:
            remove(c)


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
