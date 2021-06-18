import unittest
from qebil.tools.util import setup_output_dir
from os import path


_THIS_DIR, _THIS_FILENAME = path.split(__file__)


_TEST_SUPPORT_DIR = path.join(_THIS_DIR, "..", "support_files")


_TEST_OUTPUT_DIR = path.join(_THIS_DIR, "..", "test_output/")


setup_output_dir(_TEST_OUTPUT_DIR)


class normalizeTest(unittest.TestCase):
    """
    Tests for functions in the normalize module.
    """

    @classmethod
    def setUpClass(cls):
        pass  # TODO

    @classmethod
    def tearDownClass(cls):
        pass  # TODO

    def setUp(self):
        pass  # TODO

    def tearDown(self):
        pass  # TODO

    def test_normalize(self):
        raise NotImplementedError()  # TODO: test normalize

    def test_normalize_metadata(self):
        raise NotImplementedError()  # TODO: test normalize_metadata


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
