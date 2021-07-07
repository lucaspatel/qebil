from os import path, remove
import unittest

from qebil.tools.util import setup_output_dir, compare_checksum

from qebil.tools.fastq import (
    get_read_count,
    unpack_fastq_ftp,
    check_valid_fastq,
    remove_index_read_file,
    get_fastq_head,
    fastq_to_fasta,
    blast_for_type
)


_THIS_DIR, _THIS_FILENAME = path.split(__file__)


_TEST_SUPPORT_DIR = path.join(_THIS_DIR, "..", "support_files")


_TEST_OUTPUT_DIR = path.join(_THIS_DIR, "..", "test_output/")


setup_output_dir(_TEST_OUTPUT_DIR)


_CLEANUP_LIST = []


class fastqTest(unittest.TestCase):
    """
    Tests for functions in the fastq module.
    """

    def tearDown(self):
        self.maxDiff = None
        # clean up the directory at the end
        for c in _CLEANUP_LIST:
            if path.isfile(c):
                remove(c)

    def test_check_valid_fastq(self):
        test_fastq_path_r1 = _TEST_SUPPORT_DIR + "/SRR13874871.fastq.gz"
        test_corrupt_fastq_path = (
            _TEST_SUPPORT_DIR + "/corrupt_fastq1.fastq.gz"
        )

        self.assertTrue(check_valid_fastq(test_fastq_path_r1))
        self.assertFalse(check_valid_fastq(test_corrupt_fastq_path))

    def test_get_read_count(self):
        test_fastq_path_r1 = _TEST_SUPPORT_DIR + "/SRR13874871.fastq.gz"
        test_fastq_path_r2 = _TEST_SUPPORT_DIR + "/SRR13874871.fastq.gz"
        test_read_count_expected = "1125866"

        test_read_count = get_read_count(
            test_fastq_path_r1, test_fastq_path_r2
        )

        self.assertEqual(test_read_count, test_read_count_expected)

    def test_unpack_fastq_ftp(self):
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
        test_unpack_result = unpack_fastq_ftp(
            test_fastq_ftp_string, test_fastq_md5_string
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
        self.assertEqual(test_unpack_result, test_read_dict)


def test_remove_index_read_file(self):
    # set up files
    test_fastq_path = _TEST_SUPPORT_DIR + "/SRR13874871.fastq.gz"
    test_list = []
    for r in range(1, 4):
        r_clone = test_fastq_path.replace(
            ".fastq.gz", "_R" + str(r) + ".fastq.gz"
        )
        get_fastq_head(test_fastq_path, r_clone, r * 10)
        test_list.append(r_clone)

    # run test
    remove_index_read_file(test_list)

    # confirm behavior
    self.assertTrue(path.isfile(test_list[0]))
    self.assertEqual(get_read_count(test_list[0], 20))
    self.assertTrue(path.isfile(test_list[1]))
    self.assertEqual(get_read_count(test_list[1], 30))
    self.assertFalse(path.isfile(test_list[2]))

    # cleanup
    _CLEANUP_LIST += test_list


def test_get_fastq_head(self):
    # set up files
    test_fastq_path = _TEST_SUPPORT_DIR + "/SRR13874871.fastq.gz"
    test_fastq_path_10 = _TEST_SUPPORT_DIR + "/test_fastq_head_10.fastq.gz"
    test_fastq_path_100 = _TEST_SUPPORT_DIR + "/test_fastq_head_100.fastq.gz"
    test_fastq_path_100_md5 = "2d6c47eeca6ade6e1497e70109c7bf99"

    # run tests
    head_fq_path = get_fastq_head(test_fastq_path)
    self.assertEqual(
        head_fq_path, test_fastq_path.replace(".fastq.gz", ".head.fastq.gz")
    )
    self.assertTrue(path.isfile(head_fq_path))
    self.assertEqual(get_read_count(head_fq_path), 100)
    self.assertEqual(
        test_fastq_path_100_md5,
        compare_checksum(head_fq_path, test_fastq_path_100_md5),
    )
    self.assertFalse(
        compare_checksum(test_fastq_path_10, test_fastq_path_100_md5)
    )

    # cleanup
    _CLEANUP_LIST.append(head_fq_path)


    def test_fastq_to_fasta(self):
        test_fastq_path_10 = _TEST_SUPPORT_DIR + "/test_fastq_head_10.fastq.gz"
        test_fasta_path = _TEST_SUPPORT_DIR + "/test_fastq_head_10.fasta"
        test_fasta_path_md5 = "18a56513d251aeb59a8af3a427afcfbe"

        # run tests
        fa_path = fastq_to_fasta(test_fastq_path_10)
        self.assertEqual(fa_path, test_fasta_path)
        self.assertTrue(path.isfile(fa_path))
        self.assertEqual(get_read_count(fa_path), "fqtools error")
        self.assertEqual(
            test_fasta_path_md5, compare_checksum(fa_path, test_fasta_path_md5)
        )

        # cleanup
        _CLEANUP_LIST.append(fa_path)

    def test_blast_for_type(fastq_file, db_dict={}):
        test_local_fastq_path = (
            _TEST_SUPPORT_DIR + "/SRR13874871.fastq.gz"
        )  # AMBIGUOUS
        test_16S_fq_file = (
            _TEST_SUPPORT_DIR + "/test_16S_file.fastq.gz"
        )  # 16S result
        test_ITS_fq_file = (
            _TEST_SUPPORT_DIR + "/test_ITS_file.fastq.gz"
        )  # ITS resul

        self.assertEqual("16S", blast_for_type(test_16S_fq_file))
        self.assertEqual("ITS", blast_for_type(test_ITS_fq_file))
        self.assertEqual("AMBIGUOUS", blast_for_type(test_local_fastq_path))

if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
