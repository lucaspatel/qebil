from click.testing import CliRunner
from qebil.commands import cli
from qebil.commands.fetch import fetch
from qebil.commands.search import search
import unittest


class CliTest(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()

    def tearDown(self):
        self.runner = None

    # https://click.palletsprojects.com/en/7.x/testing/#basic-testing
    def test_cli(self):
        result1 = CliRunner().invoke(cli, ["--help"])
        self.assertEqual(result1.exit_code, 0)

    def test_cli_fetch(self):
        result2 = CliRunner().invoke(fetch, ["--help"])
        self.assertEqual(result2.exit_code, 0)

    def test_cli_fetch_project(self):
        result3 = CliRunner().invoke(fetch, ["project", "--help"])
        self.assertEqual(result3.exit_code, 0)

    def test_cli_search(self):
        result4 = CliRunner().invoke(search, ["--help"])
        self.assertEqual(result4.exit_code, 0)

    def test_cli_search_ebi(self):
        result5 = CliRunner().invoke(search, ["ebi", "--help"])
        self.assertEqual(result5.exit_code, 0)


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
    print(
        "known assertion error for exit_codes. Need to debug but functional"
    )
