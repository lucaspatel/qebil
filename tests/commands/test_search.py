from qebil.commands.search import keyword_search, search_ebi
from qebil.commands.search import search
import unittest
from click.testing import CliRunner
from os import path, remove
import glob

_this_dir, this_filename = path.split(__file__)
_test_support_dir = path.join(_this_dir, "..", "support_files")
_test_output_dir = path.join(_this_dir, "..", "test_output/")

import pandas as pd


class TestSearch(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        # clean up the directory at the start
        cleanup_list = glob.glob(_test_output_dir + "/*.EBI_metadata.tsv")
        for c in cleanup_list:
            remove(c)

    def test_keyword_search(self):
        query = "PRJNA660883"
        test_selection_dict = {}
        test_selection_dict["library_strategy"] = [
            "amplicon",
            "other",
            "wgs",
            "rna-seq",
        ]
        test_selection_dict["instrument_platform"] = ["illumina"]
        test_selection_dict["library_selection"] = [
            "random",
            "pcr",
            "random pcr",
            "cdna_randompriming",
            "inverse rrna",
            "inverse rrna selection",
            "other",
            "unspecified",
            "size fractionation",
        ]
        test_selection_dict["library_source"] = [
            "genomic",
            "genomic single cell",
            "transcriptomic",
            "transcriptomic single cell",
            "metagenomic",
            "metatranscriptomic",
            "viral rna",
            "other",
        ]
        test_summary = [
            "study_title",
            "scientific_name",
            "library_source",
            "library_strategy",
            "library_selection",
            "instrument_platform",
        ]
        test_result_0 = keyword_search(
            query, filter_dict={}, summarize=test_summary
        )
        test_result = keyword_search(
            query, filter_dict=test_selection_dict, summarize=test_summary
        )
        self.assertEqual(len(test_result_0), 1)
        self.assertEqual(len(test_result), 1)
        self.assertCountEqual(
            test_result.columns, ["study_id", "samples"] + test_summary
        )

    def test_keyword_search_method_defaults(self):
        query = "PRJNA660883"
        test_selection_dict = {}
        test_result = keyword_search(query)
        self.assertEqual(len(test_result), 1)
        self.assertCountEqual(test_result.columns, ["study_id"])

    def test_search_ebi(self):
        query = "PRJNA660883"
        output_dir = _test_output_dir
        prefix = "test"
        strategy = ["amplicon", "other", "wgs", "rna-seq"]
        platform = ["illumina"]
        selection = [
            "random",
            "pcr",
            '"random pcr"',
            "cdna_randompriming",
            '"inverse rrna"',
            '"inverse rrna selection"',
            "other",
            "unspecified",
            '"size fractionation"',
        ]
        source = [
            "genomic",
            '"genomic single cell"',
            "transcriptomic",
            '"transcriptomic single cell"',
            "metagenomic",
            "metatranscriptomic",
            '"viral rna"',
            "other",
        ]
        summarize = ["study_title"]
        scientific_name = ['"gut metagenome"']
        search_args = [
            "ebi",
            "--query",
            query,
            "--output-dir",
            output_dir,
            "--prefix",
            prefix,
        ]
        for st in strategy:
            search_args += ["--strategy", st]
        for p in platform:
            search_args += ["--platform", p]
        for se in selection:
            search_args += ["--selection", se]
        for so in source:
            search_args += ["--source", so]
        for su in summarize:
            search_args += ["--summarize", su]
        for n in scientific_name:
            search_args += ["--scientific-name", n]

        search_args = " ".join(search_args)
        test_search = CliRunner().invoke(search, search_args)
        self.assertEqual(test_search.exit_code, 0)
        test_result = pd.read_table(
            output_dir + "/" + prefix + "_EBI_search_results.tsv"
        )
        self.assertEqual(len(test_result), 1)
        self.assertCountEqual(
            test_result.columns, ["study_id", "samples"] + summarize
        )

    def test_search_ebi_no_summary(self):
        query = "PRJNA660883"
        output_dir = _test_output_dir
        prefix = "test"
        strategy = ["amplicon", "other", "wgs", "rna-seq"]
        platform = ["illumina"]
        selection = [
            "random",
            "pcr",
            '"random pcr"',
            "cdna_randompriming",
            '"inverse rrna"',
            '"inverse rrna selection"',
            "other",
            "unspecified",
            '"size fractionation"',
        ]
        source = [
            "genomic",
            '"genomic single cell"',
            "transcriptomic",
            '"transcriptomic single cell"',
            "metagenomic",
            "metatranscriptomic",
            '"viral rna"',
            "other",
        ]
        summarize = []
        scientific_name = ['"gut metagenome"']
        search_args = [
            "ebi",
            "--query",
            query,
            "--output-dir",
            output_dir,
            "--prefix",
            prefix,
        ]
        for st in strategy:
            search_args += ["--strategy", st]
        for p in platform:
            search_args += ["--platform", p]
        for se in selection:
            search_args += ["--selection", se]
        for so in source:
            search_args += ["--source", so]
        for su in summarize:
            search_args += ["--summarize", su]
        for n in scientific_name:
            search_args += ["--scientific-name", n]

        search_args = " ".join(search_args)
        test_search = CliRunner().invoke(search, search_args)
        self.assertEqual(test_search.exit_code, 0)
        test_result = pd.read_table(
            output_dir + "/" + prefix + "_EBI_search_results.tsv"
        )

        test_result = pd.read_table(
            output_dir + prefix + "_EBI_search_results.tsv"
        )
        self.assertEqual(len(test_result), 1)
        self.assertCountEqual(
            test_result.columns,
            [
                "study_id",
                "samples",
                "study_title",
                "scientific_name",
                "library_source",
                "library_strategy",
                "library_selection",
                "instrument_platform",
            ],
        )


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
