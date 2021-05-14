import unittest
import socket
from os import path, remove
from qebil.log import (
    setup_log,
    setup_logging,
    host_log_adapter,
    get_timestamp,
)
import logging

from qebil.tools.util import setup_output_dir

_THIS_DIR, _THIS_FILENAME = path.split(__file__)

_TEST_SUPPORT_DIR = path.join(_THIS_DIR, "support_files")

_TEST_OUTPUT_DIR = path.join(_THIS_DIR, "test_output/")

setup_output_dir(_TEST_OUTPUT_DIR)


class TestSetupLogging(unittest.TestCase):
    def test_setup_log(self):
        """Tests whether a log file with the excpected contents is created"""

        prefix = "test"
        suffix = "example"
        quiet = False
        test_log_file = _TEST_OUTPUT_DIR + prefix + suffix + ".log"

        if path.isfile(test_log_file):
            remove(test_log_file)

        setup_log(_TEST_OUTPUT_DIR, prefix, suffix, quiet)
        from qebil.log import logger

        test_log_line = "This is a test."
        logger.info(test_log_line)

        self.assertEqual(path.isfile(test_log_file), True)

        read_test_log = open(test_log_file, "r")
        test_log_lines = read_test_log.readlines()
        read_test_log.close()
        print("Log contents: " + str(test_log_lines))

        self.assertEqual(len(test_log_lines), 1)
        self.assertTrue(test_log_line in test_log_lines[-1])

    def test_setup_logging_defaults(self):
        """Test empty setup of logging"""
        # setup up logger
        setup_logging()
        from qebil.log import logger

        self.assertEqual(
            logger.getEffectiveLevel(),
            logging.LoggerAdapter(
                logging.getLogger(), {"hostname": socket.gethostname()}
            ).getEffectiveLevel(),
        )

    def test_setup_logging_file_handler(self):
        """Test setup of the logging file handler with a dummy log

        # TODO: add test for level not on list and ERROR

        """
        level_list = ["DEBUG", "INFO", "WARNING"]

        for l in level_list:
            test_path = _TEST_OUTPUT_DIR + l.lower() + ".log"

            if path.isfile(test_path):
                remove(test_path)

            setup_logging(test_path, l)
            from qebil.log import logger

            logger.error(l.lower())
            test_log_expected_contents = l.lower() + "\n"

            # need to get timestamp format to update in string to compare
            test_log_file = open(test_path, "r")
            test_log_contents = test_log_file.read()
            test_log_file.close()
            self.assertEqual(test_log_expected_contents, test_log_contents)


def test_get_timestamp(self):
    test_time = get_timestamp()
    self.assertEqual(len(test_time), 15)


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
