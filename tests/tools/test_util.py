from qebil.tools.util import (
    load_project_file,
    retrieve_ftp_file,
    check_download_integrity,
    scrape_ebi_ids,
    parse_document,
)
import unittest
from os import path, remove
from shutil import copy
import glob

from qebil.tools.util import setup_output_dir

_THIS_DIR, _THIS_FILENAME = path.split(__file__)

_TEST_SUPPORT_DIR = path.join(_THIS_DIR, "..", "support_files")

_TEST_OUTPUT_DIR = path.join(_THIS_DIR, "..", "test_output/")

setup_output_dir(_TEST_OUTPUT_DIR)


class utilTest(unittest.TestCase):
    """
    Tests for functions in the util module.
    """

    def setUp(self):
        self.maxDiff = None
        # clean up the directory at the start
        cleanup_list = glob.glob(_TEST_OUTPUT_DIR + "/*.EBI_metadata.tsv")
        for c in cleanup_list:
            remove(c)

    def test_load_project_file(self):
        test_project_file = _TEST_SUPPORT_DIR + "/test_project_file.tsv"
        test_project_file_no_header = (
            _TEST_SUPPORT_DIR + "/test_project_file_no_header.tsv"
        )
        test_load_project_list = load_project_file(test_project_file)
        test_load_project_list_no_header = load_project_file(
            test_project_file_no_header
        )

        self.assertEqual(len(test_load_project_list), 2)
        self.assertEqual(len(test_load_project_list_no_header), 2)
        self.assertEqual(
            test_load_project_list, test_load_project_list_no_header
        )

    def test_retrieve_ftp_file(self):
        test_ftp_path = "ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/071/SRR13874871/SRR13874871.fastq.gz"
        test_local_fastq_path = _TEST_OUTPUT_DIR + "/SRR13874871.fastq.gz"
        test_checksum = "06445ed5341e3779ac1d5230c787c538"
        test_corrupt_fastq_path = (
            _TEST_SUPPORT_DIR + "/corrupt_fastq1.fastq.gz"
        )
        test_overwrite = True

        test_checksum_1 = retrieve_ftp_file(
            test_ftp_path, test_local_fastq_path, test_checksum
        )

        self.assertTrue(path.isfile(test_local_fastq_path))
        self.assertEqual(test_checksum, test_checksum_1)

        copy(test_corrupt_fastq_path, test_local_fastq_path)

        test_checksum_2 = retrieve_ftp_file(
            test_ftp_path, test_local_fastq_path, test_checksum
        )
        self.assertTrue(path.isfile(test_local_fastq_path))
        self.assertEqual(test_checksum, test_checksum_2)

    def test_check_download_integrity(self):
        test_local_fastq_path = _TEST_SUPPORT_DIR + "/SRR13874871.fastq.gz"
        expected_md5 = "06445ed5341e3779ac1d5230c787c538"
        self.assertTrue(path.isfile(test_local_fastq_path))
        md5checksum_res = check_download_integrity(
            test_local_fastq_path, expected_md5
        )
        self.assertEqual(expected_md5, md5checksum_res)

    def test_parse_document(self):
        tokens_html = parse_document(
            "https://www.nature.com/articles/s41598-020-60564-8"
        )
        tokens_pdf = parse_document(
            "https://www.nature.com/articles/s41598-021-83922-6.pdf"
        )

        self.assertEqual(len(tokens_html), 84361) # this changes a lot...
        self.assertEqual(len(tokens_pdf), 15771)

    def scrape_ebi_ids(self):
        test_paths = [
            "https://www.nature.com/articles/s41598-021-83922-6.pdf",
            "https://www.nature.com/articles/s41598-020-60564-8",
            _TEST_SUPPORT_DIR + "/test.pdf",
        ]
        test_ids = [["PRJNA61366"], [], ["PRJNA61366"]]
        test_ids_retrieved = []

        for f in test_paths:
            parsed_doc = parse_document(f)
            test_ids_retrieved.append(scrape_ebi_ids(parsed_doc))
            self.assertEqual(test_ids, test_ids_retrieved)


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
