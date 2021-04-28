from qebil.process import deplete_on_the_fly, run_fastp, run_host_depletion
import unittest
from os import path
from shutil import copy
import pandas as pd

this_dir, this_filename = path.split(__file__)
_test_support_dir = path.join(this_dir, "support_files")
_test_output_dir = path.join(this_dir, "test_output/")
test_study_id='SRP283872'

class ProcessTest(unittest.TestCase):
    
    # https://click.palletsprojects.com/en/7.x/testing/#basic-testing
    def test_deplete_on_the_fly(self):
        """Tests on-the-fly host depletion and recovery

        """            
        result_md = pd.DataFrame()
        # need to create a study and supply the test_output dir, with and without keeping files
        # and make sure that the returned result dictionary has been populated with the expected
        # read counts and fractions of quality and microbial reads
        # for now, will read in expected final dataframe to compare, and run from remote with
        # small study. Need to also load from local partial files with only the first file finished
        # and second file still to do as well as with only a partially complete first file
        # need to also make sure that this behaves for mixed data types so that only metaG an metaT
        # files are processed, not amplicon. To do this, we'll load from a local file that is a mix
        # of the short metaG study and the short amplicon study used elsewhere in the tests
        raise NotImplementedError


    def test_run_fastp(self):
        """ Tests the running of fastp"""

        # need to supply a run_prefix, a list of raw_reads, instrument model and an output location
        # test keeping the raw files and removing them
        #raise NotImplementedError
        copy('test file','temp test file')
        paired_reads_prefix = ''
        read_1_only_prefix = ''
        too_many_reads_prefix = ''
        paired_read_count = 1
        read_1_only_count = 1
        expected_filtered_paired_read_count = 1
        expected_filtered_read_1_only_count = 1
        corrupt_file_prefix = ''

        self.assertEqual(None, run_fastp(too_many_reads_prefix,0,keep_files=True))
        self.assertEqual(None, run_fastp(paired_reads_prefix,0,keep_files=True))
        self.assertEqual('not determined',run_fastp(corrupt_file_prefix,0,keep_files=True))

        self.assertEqual(expected_filtered_paired_read_count,run_fastp(paired_reads_prefix,paired_read_count,keep_raw=True))
        self.assertEqual(expected_filtered_read_1_only_count,run_fastp(read_1_only_prefix,read_1_only_count,keep_raw=True))    
        self.assertTrue(path.isdir(_test_output_dir + "fastp_reports/"))
        self.assertTrue(path.isfile('temp test file'))
        self.assertTrue(path.isfile('temp test file'))

        self.assertEqual(expected_filtered_paired_read_count,run_fastp(paired_reads_prefix,paired_read_count))
        self.assertEqual(expected_filtered_read_1_only_count,run_fastp(read_1_only_prefix,read_1_only_count))
        self.assertFalse(path.isfile('temp test file'))
        self.assertFalse(path.isfile('temp test file'))


    def test_run_host_depletion(self):
        """Tests host depletion of fastq files
        """
        #raise NotImplementedError
        copy('test file','temp test file')
        paired_reads_prefix = ''
        read_1_only_prefix = ''
        too_many_reads_prefix = ''
        paired_read_count = 1
        read_1_only_count = 1
        expected_filtered_paired_read_count = 1
        expected_filtered_read_1_only_count = 1
        corrupt_file_prefix = ''

        self.assertEqual(None, run_host_depletion(too_many_reads_prefix,0))
        self.assertEqual(None, run_host_depletion(paired_reads_prefix,0))
        self.assertEqual('not determined',run_host_depletion(corrupt_file_prefix,0))
        self.assertTrue(path.isdir(_test_output_dir + "fastp_reports/"))

        self.assertEqual(expected_filtered_paired_read_count,run_host_depletion(paired_reads_prefix,paired_read_count, keep=True))
        self.assertEqual(expected_filtered_read_1_only_count,run_host_depletion(read_1_only_prefix,read_1_only_count, keep=True))              
        self.assertTrue(path.isfile('temp test file'))
        self.assertTrue(path.isfile('temp test file'))
        self.assertFalse(path.isfile(paired_reads_prefix + "_minimap2.start"))

        self.assertEqual(expected_filtered_paired_read_count,run_fastp(paired_reads_prefix,paired_read_count))
        self.assertEqual(expected_filtered_read_1_only_count,run_fastp(read_1_only_prefix,read_1_only_count))
        self.assertFalse(path.isfile('temp test file'))
        self.assertFalse(path.isfile('temp test file'))
                     
        
if __name__ == '__main__':
    # begin the unittest.main()
    unittest.main()
