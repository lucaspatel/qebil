from qebil.tools.fastq import get_read_count, unpack_fastq_ftp, check_valid_fastq, check_fastq_tail
import unittest
from os import path

this_dir, this_filename = path.split(__file__)
_test_support_dir = path.join(this_dir, "..", "support_files")
_test_output_dir = path.join(this_dir, "..","test_output/")

class fastqTest(unittest.TestCase):
    """
    Tests for functions in the fastq module.
    """
    def setUp(self):
        self.maxDiff = None
        # clean up the directory at the start
        cleanup_list = glob.glob(_test_output_dir + '/*.EBI_metadata.tsv')
        for c in cleanup_list:
            remove(c)
    
    def test_check_valid_fastq(self):
        test_fastq_path_r1 = _test_support_dir + '/SAMN16049500.SRR12672280.R1.ebi.fastq.gz'
        test_corrupt_fastq_path = _test_support_dir + '/corrupt_fastq1.fastq.gz'
        
        self.assertTrue(check_valid_fastq(test_fastq_path_r1))
        self.assertFalse(check_valid_fastq(test_corrupt_fastq_path))
    
    def test_get_read_count(self):
        test_fastq_path_r1 = _test_support_dir + '/SAMN16049500.SRR12672280.R1.ebi.fastq.gz'
        test_fastq_path_r2 = _test_support_dir + '/SAMN16049500.SRR12672280.R2.ebi.fastq.gz'
        test_read_count_expected = 11481502
        
        test_read_count = get_read_count(test_fastq_path_r1, test_fastq_path_r2)
        
        self.assertEqual(test_read_count, test_read_count_expected)

    
    def test_unpack_fastq_ftp(self):
        test_fastq_ftp_string = ('ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/'
                                +'SRR12672280/SRR12672280_1.fastq.gz;'
                                +'ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/'
                                + 'SRR12672280/SRR12672280_2.fastq.gz'
                                )
        test_fastq_md5_string = ('5c1da3b86d2bbb0d09e1f05cef0107f2;'
                                + 'fe207ea59d80b5143e142050e37bbd11')
        test_unpack_result = unpack_fastq_ftp(test_fastq_ftp_string, test_fastq_md5_string)
        test_read_dict = {'read_1':{'ftp': test_fastq_ftp_string.split(';')[0],
                                    'md5': test_fastq_md5_string.split(';')[0]},
                          'read_2':{'ftp': test_fastq_ftp_string.split(';')[1],
                                    'md5': test_fastq_md5_string.split(';')[1]}
                          }
        self.assertEqual(test_unpack_result,test_read_dict)
        
    def test_check_fastq_tail(self):
        test_fastq_path_r1 = _test_support_dir + '/SAMN16049500.SRR12672280.R1.ebi.fastq.gz'
        test_corrupt_fastq_path = _test_support_dir + '/corrupt_fastq1.fastq.gz'
        
        self.assertTrue(check_fastq_tail(test_fastq_path_r1))
        self.assertFalse(check_fastq_tail(test_corrupt_fastq_path))
                          
if __name__ == '__main__':
    # begin the unittest.main()
    unittest.main()