from qebil.commands.fetch import (
    fetch_remote_studies,
    check_existing_metadata,
    fetch,
)
import unittest
from os import path, remove
import glob
from click.testing import CliRunner
from qebil.tools.util import setup_output_dir

_THIS_DIR, _THIS_FILENAME = path.split(__file__)

_TEST_SUPPORT_DIR = path.join(_THIS_DIR, "..", "support_files")

_TEST_OUTPUT_DIR = path.join(_THIS_DIR, "..", "test_output/")

setup_output_dir(_TEST_OUTPUT_DIR)


class fetchTest(unittest.TestCase):
    """
    Tests for functions in the fetch module.
    """

    def setUp(self):
        self.maxDiff = None
        # clean up the directory at the start
        cleanup_list = glob.glob(_TEST_OUTPUT_DIR + "/*.EBI_metadata.tsv")
        for c in cleanup_list:
            remove(c)

    def tearDown(self):
        pass  # TODO

    def test_fetch_remote_studies(self):
        """Standard test that should work"""
        test_study_list = [
            "PRJNA660883",
            "SRP283872",
        ]  # note, these should resolve to the same study
        test_full_details = True
        test_max_samples = ""
        test_max_samples_num_str = "2"
        test_max_samples_num = 2
        test_random_subsample = False

        result_dict = fetch_remote_studies(
            test_study_list,
            test_full_details,
            test_max_samples,
            test_random_subsample,
            _TEST_OUTPUT_DIR
        )
        self.assertEqual(
            result_dict[test_study_list[0]].proj_id, test_study_list[0]
        )
        self.assertEqual(len(result_dict[test_study_list[0]].metadata), 36)
        self.assertEqual(
            result_dict[test_study_list[1]].study_id, test_study_list[1]
        )
        self.assertEqual(len(result_dict[test_study_list[1]].metadata), 36)
        self.assertNotEqual(
            result_dict[test_study_list[0]].proj_id, test_study_list[1]
        )

        results_dict_1 = fetch_remote_studies(
            test_study_list, test_full_details,output_dir=_TEST_OUTPUT_DIR
        )
        self.assertCountEqual(result_dict.keys(), results_dict_1.keys())

        results_dict_2 = fetch_remote_studies(
            test_study_list,
            test_full_details,
            test_max_samples_num,
            test_random_subsample,
            _TEST_OUTPUT_DIR
        )
        results_dict_3 = fetch_remote_studies(
            test_study_list,
            test_full_details,
            test_max_samples_num_str,
            test_random_subsample,
            _TEST_OUTPUT_DIR
        )

        self.assertCountEqual(results_dict_2.keys(), results_dict_3.keys())

        results_dict_4 = fetch_remote_studies(
            test_study_list, test_full_details, 4, True, _TEST_OUTPUT_DIR
        )

        self.assertCountEqual(results_dict_2.keys(), results_dict_4.keys())

    def test_fetch_project(self):
        test_ebi_id = "PRJNA660883"
        test_project_file = _TEST_SUPPORT_DIR + "/test_project_file.tsv"
        test_metadata_file = (
            _TEST_SUPPORT_DIR + "/test_study.EBI_metadata.tsv"
        )
        test_output_dir = _TEST_OUTPUT_DIR
        test_prefix = "test_fetch_project"
        test_prep_max_num_str = "250"
        test_cpus_num_str = "4"
        test_max_samples_num_str = "2"
        fetch_args = [
            "--output-dir",
            test_output_dir,
            "--prefix",
            test_prefix,
            "--prep-max",
            test_prep_max_num_str,
            "--cpus",
            test_cpus_num_str,
            "--max-samples",
            test_max_samples_num_str,
        ]
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
        for st in strategy:
            fetch_args += ["--strategy", st]
        for p in platform:
            fetch_args += ["--platform", p]
        for se in selection:
            fetch_args += ["--selection", se]
        for so in source:
            fetch_args += ["--source", so]
        for n in scientific_name:
            fetch_args += ["--scientific-name", n]

        fetch_proj_args = ["project", "--ebi-id", test_ebi_id] + fetch_args
        fetch_file_args = [
            "project",
            "--project-file",
            test_project_file,
        ] + fetch_args
        fetch_md_args = [
            "project",
            "--metadata-file",
            test_metadata_file,
        ] + fetch_args

        fetch_proj_arg_string = " ".join(fetch_proj_args)
        print(fetch_proj_arg_string)
        test_fetch_proj = CliRunner().invoke(fetch, fetch_proj_arg_string)
        self.assertEqual(test_fetch_proj.exit_code, 0)

        fetch_file_arg_string = " ".join(fetch_file_args)
        print(fetch_file_arg_string)
        test_fetch_file = CliRunner().invoke(fetch, fetch_file_arg_string)
        self.assertEqual(test_fetch_file.exit_code, 0)

        fetch_md_arg_string = " ".join(fetch_md_args)
        test_fetch_md = CliRunner().invoke(fetch, fetch_md_arg_string)
        self.assertEqual(test_fetch_md.exit_code, 0)
        # TODO more comprehensive checks

        # TODO add tests for these
        test_metadata_base = (
            _TEST_SUPPORT_DIR + "/SRP116878_sample_info.EBI_metadata.tsv"
        )
        test_merge_metadata = _TEST_SUPPORT_DIR + "/test_merger_metadata.tsv"
        test_merge_col = "library_name"
        fetch_augment_args = (
            ["project", "--metadata-file", test_metadata_base]
            + fetch_proj_args
            + [
                "--add_metadata_file",
                test_merge_metadata,
                "--merge_column",
                test_merge_col,
            ]
        )

        # additional valid configurations to test
        test_quiet_yes = "--quiet"
        test_qebil_yes = "--qiita"
        test_qebil_no = "--raw"
        test_download_fastq_yes = "--download-fastq"
        test_human_removal_yes = "--human-removal"
        test_keep_files_yes = "--keep-files"
        test_no_filter_yes = "--no-filter"
        test_random_subsample_yes = "random-subsample"
        test_overwrite = "--overwrite"
        test_output_dir_int = 2
        test_output_dir_no_slash = _TEST_OUTPUT_DIR
        test_prefix_int = 2
        test_prep_max = 250
        test_cpus_int = 4
        test_max_samples_num = 2
        test_add_emp = "--emp-protocol"

        # error modes to handle
        test_prep_max_str = "error"
        test_cpus_str = "error"
        test_max_samples = ""


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
