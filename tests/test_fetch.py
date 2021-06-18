import unittest
from pandas.testing import assert_frame_equal
from collections import OrderedDict
from qebil.fetch import (
    fetch_ebi_info,
    fetch_fastq_files,
    fetch_fastqs,
    fetch_ebi_metadata,
)
from qebil.core import Study
import glob
from os import remove, makedirs, path
import pandas as pd
from qebil.log import setup_log

from qebil.tools.util import setup_output_dir

_THIS_DIR, _THIS_FILENAME = path.split(__file__)

_TEST_SUPPORT_DIR = path.join(_THIS_DIR, "support_files")

_TEST_OUTPUT_DIR = path.join(_THIS_DIR, "test_output/")

setup_output_dir(_TEST_OUTPUT_DIR)


def test_stage_output():
    if not path.isdir(_TEST_OUTPUT_DIR):
        makedirs(_TEST_OUTPUT_DIR)
    else:
        cleanup_list = glob.glob(_TEST_OUTPUT_DIR + "/*")
        for c in cleanup_list:
            remove(c)


class FetchTest(unittest.TestCase):
    def setUp(self):
        prefix = "test"
        suffix = "example"
        quiet = False
        test_log_file = _TEST_OUTPUT_DIR + prefix + suffix + ".log"

        if path.isfile(test_log_file):
            remove(test_log_file)

        setup_log(_TEST_OUTPUT_DIR, prefix, suffix, quiet)
        from qebil.log import logger

    def test_fetch_fastq_files(self):
        test_run_prefix = "SAMN16049500.SRR12672280"
        test_fastq_ftp_string = (
            "ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/"
            + "SRR12672280/SRR12672280_1.fastq.gz;"
            + "ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/"
            + "SRR12672280/SRR12672280_2.fastq.gz"
        )
        test_fastq_md5_string = (
            "5c1da3b86d2bbb0d09e1f05cef0107f2;"
            + "fe207ea59d80b5143e142050e37bbd11"
        )
        test_read_dict = {
            "read_1": {
                "ftp": test_fastq_ftp_string.split(";")[0],
                "md5": test_fastq_md5_string.split(";")[0],
            },
            "read_2": {
                "ftp": test_fastq_ftp_string.split(";")[1],
                "md5": test_fastq_md5_string.split(";")[1],
            },
        }
        test_expected_reads = "5740751"
        test_raw_reads = fetch_fastq_files(
            test_run_prefix, test_read_dict, _TEST_OUTPUT_DIR
        )
        self.assertEqual(str(test_raw_reads), test_expected_reads)

        # TODO: find corrupted downloads and study with > 2 reads to test error message

    def test_fetch_fastqs(self):
        test_stage_output()  # run this before any test that generates output files
        test_study_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/test_study.EBI_metadata.tsv",
            sep="\t",
            header=0,
            index_col=0,
        ).head(1)
        test_study = Study(test_study_df)
        test_expected_files = [
            "SAMN16049500.SRR12672280.R1.ebi.fastq.gz",
            "SAMN16049500.SRR12672280.R2.ebi.fastq.gz",
        ]
        test_expected_reads = "5740751"
        test_result_md = fetch_fastqs(test_study, _TEST_OUTPUT_DIR)
        print(test_result_md)
        test_result_files = [
            f.split("/")[-1]
            for f in sorted(glob.glob(_TEST_OUTPUT_DIR + "/*.fastq.gz"))
        ]
        self.assertEqual(
            str(test_result_md["qebil_raw_reads"].sum()), test_expected_reads
        )
        self.assertEqual(test_result_files, test_expected_files)

    def test_fetch_ebi_info_study(self):
        accession = "SRP283872"
        test_study_dict = OrderedDict(
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

        test_xml_dict = fetch_ebi_info(accession)
        self.assertEqual(test_xml_dict, test_study_dict)

    def test_fetch_ebi_info_sample(self):
        accession = "SRS7392559"
        test_sample_dict = OrderedDict(
            [
                (
                    "SAMPLE_SET",
                    OrderedDict(
                        [
                            (
                                "SAMPLE",
                                OrderedDict(
                                    [
                                        ("@accession", "SRS7392559"),
                                        ("@alias", "MH001880_V1"),
                                        ("@center_name", "pda|jeremiahfaith"),
                                        ("@broker_name", "NCBI"),
                                        (
                                            "IDENTIFIERS",
                                            OrderedDict(
                                                [
                                                    (
                                                        "PRIMARY_ID",
                                                        "SRS7392559",
                                                    ),
                                                    (
                                                        "EXTERNAL_ID",
                                                        OrderedDict(
                                                            [
                                                                (
                                                                    "@namespace",
                                                                    "BioSample",
                                                                ),
                                                                (
                                                                    "#text",
                                                                    "SAMN16049500",
                                                                ),
                                                            ]
                                                        ),
                                                    ),
                                                    (
                                                        "SUBMITTER_ID",
                                                        OrderedDict(
                                                            [
                                                                (
                                                                    "@label",
                                                                    "Sample name",
                                                                ),
                                                                (
                                                                    "@namespace",
                                                                    "pda|jeremiahfaith",
                                                                ),
                                                                (
                                                                    "#text",
                                                                    "MH001880_V1",
                                                                ),
                                                            ]
                                                        ),
                                                    ),
                                                ]
                                            ),
                                        ),
                                        ("TITLE", "gut metagenome"),
                                        (
                                            "SAMPLE_NAME",
                                            OrderedDict(
                                                [
                                                    ("TAXON_ID", "749906"),
                                                    (
                                                        "SCIENTIFIC_NAME",
                                                        "gut metagenome",
                                                    ),
                                                ]
                                            ),
                                        ),
                                        (
                                            "DESCRIPTION",
                                            "bristol:3,dayFromSymptomOnset:6,peakSeverity:Severe",
                                        ),
                                        (
                                            "SAMPLE_LINKS",
                                            OrderedDict(
                                                [
                                                    (
                                                        "SAMPLE_LINK",
                                                        [
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "XREF_LINK",
                                                                        OrderedDict(
                                                                            [
                                                                                (
                                                                                    "DB",
                                                                                    "ENA-STUDY",
                                                                                ),
                                                                                (
                                                                                    "ID",
                                                                                    "SRP283872",
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
                                                                                    "SRX9152556",
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
                                                                                    "SRR12672280",
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
                                                                                    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRS7392559&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes",
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
                                                                                    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRS7392559&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format",
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
                                            "SAMPLE_ATTRIBUTES",
                                            OrderedDict(
                                                [
                                                    (
                                                        "SAMPLE_ATTRIBUTE",
                                                        [
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "host",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "Homo sapiens",
                                                                    ),
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "isolation_source",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "human feces",
                                                                    ),
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "collection_date",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "2020",
                                                                    ),
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "geo_loc_name",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "USA:New York",
                                                                    ),
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "lat_lon",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "40.71 N 74.00 W",
                                                                    ),
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "subject_ID",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "MH001880_V1",
                                                                    ),
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "BioSampleModel",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "Metagenome or environmental",
                                                                    ),
                                                                ]
                                                            ),
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "ENA-SPOT-COUNT",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "5740751",
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
                                                                        "1722225300",
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
        test_xml_dict = fetch_ebi_info(accession)
        self.assertEqual(test_xml_dict, test_sample_dict)

    def test_fetch_ebi_info_expt(self):
        accession = "SRX9152556"
        test_xml_dict = fetch_ebi_info(accession)
        test_expt_dict = OrderedDict(
            [
                (
                    "EXPERIMENT_SET",
                    OrderedDict(
                        [
                            (
                                "EXPERIMENT",
                                OrderedDict(
                                    [
                                        ("@accession", "SRX9152556"),
                                        ("@alias", "MH001880_V1"),
                                        ("@center_name", "SUB8091541"),
                                        (
                                            "IDENTIFIERS",
                                            OrderedDict(
                                                [
                                                    (
                                                        "PRIMARY_ID",
                                                        "SRX9152556",
                                                    ),
                                                    (
                                                        "SUBMITTER_ID",
                                                        OrderedDict(
                                                            [
                                                                (
                                                                    "@namespace",
                                                                    "SUB8091541",
                                                                ),
                                                                (
                                                                    "#text",
                                                                    "MH001880_V1",
                                                                ),
                                                            ]
                                                        ),
                                                    ),
                                                ]
                                            ),
                                        ),
                                        (
                                            "TITLE",
                                            "Illumina HiSeq 4000 paired end sequencing; SARS-CoV-2-specific IgA and limited inflammatory cytokines are present in the stool of select patients with acute COVID-19",
                                        ),
                                        (
                                            "STUDY_REF",
                                            OrderedDict(
                                                [
                                                    (
                                                        "@accession",
                                                        "SRP283872",
                                                    ),
                                                    (
                                                        "IDENTIFIERS",
                                                        OrderedDict(
                                                            [
                                                                (
                                                                    "PRIMARY_ID",
                                                                    "SRP283872",
                                                                ),
                                                                (
                                                                    "EXTERNAL_ID",
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
                                                ]
                                            ),
                                        ),
                                        (
                                            "DESIGN",
                                            OrderedDict(
                                                [
                                                    (
                                                        "DESIGN_DESCRIPTION",
                                                        "DNA (extracted as for 16S rRNA amplicon sequencing) was sonicated and Illumina sequencing libraries generated using the NEBNext Ultra II DNA Library Prep Kit. Ligation products of 500 600 base pairs were purified using SPRIselect beads (Beckman Coulter) and enrichment PCR performed. Samples were pooled in equal proportions and size-selected using AMPure XP beads (Beckman Coulter) before sequencing with an Illumina HiSeq (paired-end 150 bp).",
                                                    ),
                                                    (
                                                        "SAMPLE_DESCRIPTOR",
                                                        OrderedDict(
                                                            [
                                                                (
                                                                    "@accession",
                                                                    "SRS7392559",
                                                                ),
                                                                (
                                                                    "IDENTIFIERS",
                                                                    OrderedDict(
                                                                        [
                                                                            (
                                                                                "PRIMARY_ID",
                                                                                "SRS7392559",
                                                                            ),
                                                                            (
                                                                                "EXTERNAL_ID",
                                                                                OrderedDict(
                                                                                    [
                                                                                        (
                                                                                            "@namespace",
                                                                                            "BioSample",
                                                                                        ),
                                                                                        (
                                                                                            "#text",
                                                                                            "SAMN16049500",
                                                                                        ),
                                                                                    ]
                                                                                ),
                                                                            ),
                                                                        ]
                                                                    ),
                                                                ),
                                                            ]
                                                        ),
                                                    ),
                                                    (
                                                        "LIBRARY_DESCRIPTOR",
                                                        OrderedDict(
                                                            [
                                                                (
                                                                    "LIBRARY_NAME",
                                                                    "MH001880_V1",
                                                                ),
                                                                (
                                                                    "LIBRARY_STRATEGY",
                                                                    "WGS",
                                                                ),
                                                                (
                                                                    "LIBRARY_SOURCE",
                                                                    "METAGENOMIC",
                                                                ),
                                                                (
                                                                    "LIBRARY_SELECTION",
                                                                    "size fractionation",
                                                                ),
                                                                (
                                                                    "LIBRARY_LAYOUT",
                                                                    OrderedDict(
                                                                        [
                                                                            (
                                                                                "PAIRED",
                                                                                None,
                                                                            )
                                                                        ]
                                                                    ),
                                                                ),
                                                            ]
                                                        ),
                                                    ),
                                                ]
                                            ),
                                        ),
                                        (
                                            "PLATFORM",
                                            OrderedDict(
                                                [
                                                    (
                                                        "ILLUMINA",
                                                        OrderedDict(
                                                            [
                                                                (
                                                                    "INSTRUMENT_MODEL",
                                                                    "Illumina HiSeq 4000",
                                                                )
                                                            ]
                                                        ),
                                                    )
                                                ]
                                            ),
                                        ),
                                        (
                                            "EXPERIMENT_LINKS",
                                            OrderedDict(
                                                [
                                                    (
                                                        "EXPERIMENT_LINK",
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
                                                                                    "SRS7392559",
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
                                                                                    "SRR12672280",
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
                                                                                    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRX9152556&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes",
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
                                                                                    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRX9152556&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format",
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
                                            "EXPERIMENT_ATTRIBUTES",
                                            OrderedDict(
                                                [
                                                    (
                                                        "EXPERIMENT_ATTRIBUTE",
                                                        [
                                                            OrderedDict(
                                                                [
                                                                    (
                                                                        "TAG",
                                                                        "ENA-SPOT-COUNT",
                                                                    ),
                                                                    (
                                                                        "VALUE",
                                                                        "5740751",
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
                                                                        "1722225300",
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

        self.assertEqual(test_xml_dict, test_expt_dict)

    def test_fetch_ebi_metadata(self):
        accession = "SRP283872"
        test_limited_fields = [
            "run_accession",
            "sample_accession",
            "library_name",
        ]
        test_study_valid_tsv = _TEST_SUPPORT_DIR + "/test_study_basic.tsv"
        test_expected_df = pd.read_csv(
            test_study_valid_tsv, header=0, sep="\t", dtype=str
        )
        test_result_df = fetch_ebi_metadata(accession)
        test_limited_result_df = fetch_ebi_metadata(
            accession, test_limited_fields
        )
        print(test_expected_df.columns)
        print(test_result_df.columns)
        assert_frame_equal(
            test_expected_df.sort_index(axis=1),
            test_result_df.sort_index(axis=1),
        )
        self.assertEqual(
            list(test_limited_result_df.columns), test_limited_fields
        )


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
