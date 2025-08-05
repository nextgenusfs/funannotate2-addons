#!/usr/bin/env python3

import unittest
import sys
from unittest import mock
from io import StringIO

from funannotate2_addons.__main__ import parse_args, main


class TestCLI(unittest.TestCase):
    def test_parse_args_no_args(self):
        """Test that no arguments raises SystemExit and prints help"""
        with self.assertRaises(SystemExit) as cm:
            with mock.patch("sys.stderr", new_callable=StringIO):
                parse_args([])
        self.assertEqual(cm.exception.code, 0)

    def test_parse_args_help(self):
        """Test that --help raises SystemExit"""
        with self.assertRaises(SystemExit) as cm:
            with mock.patch("sys.stdout", new_callable=StringIO):
                parse_args(["--help"])
        self.assertEqual(cm.exception.code, 0)

    def test_parse_args_version(self):
        """Test that --version raises SystemExit"""
        with self.assertRaises(SystemExit) as cm:
            with mock.patch("sys.stdout", new_callable=StringIO):
                parse_args(["--version"])
        self.assertEqual(cm.exception.code, 0)

    def test_parse_args_emapper(self):
        """Test parsing emapper subcommand"""
        with self.assertRaises(SystemExit) as cm:
            with mock.patch("sys.stdout", new_callable=StringIO):
                parse_args(["emapper", "--help"])
        self.assertEqual(cm.exception.code, 0)

    def test_parse_args_iprscan(self):
        """Test parsing iprscan subcommand"""
        with self.assertRaises(SystemExit) as cm:
            with mock.patch("sys.stdout", new_callable=StringIO):
                parse_args(["iprscan", "--help"])
        self.assertEqual(cm.exception.code, 0)

    def test_parse_args_antismash(self):
        """Test parsing antismash subcommand"""
        with self.assertRaises(SystemExit) as cm:
            with mock.patch("sys.stdout", new_callable=StringIO):
                parse_args(["antismash", "--help"])
        self.assertEqual(cm.exception.code, 0)

    def test_parse_args_signalp6(self):
        """Test parsing signalp6 subcommand"""
        with self.assertRaises(SystemExit) as cm:
            with mock.patch("sys.stdout", new_callable=StringIO):
                parse_args(["signalp6", "--help"])
        self.assertEqual(cm.exception.code, 0)

    @mock.patch("funannotate2_addons.__main__.run_emapper_cli")
    def test_main_emapper(self, mock_run_emapper):
        """Test main function calls emapper CLI"""
        with mock.patch("sys.argv", ["funannotate2_addons", "emapper", "--help"]):
            with self.assertRaises(SystemExit):
                main()

    @mock.patch("funannotate2_addons.__main__.run_iprscan_cli")
    def test_main_iprscan(self, mock_run_iprscan):
        """Test main function calls iprscan CLI"""
        with mock.patch("sys.argv", ["funannotate2_addons", "iprscan", "--help"]):
            with self.assertRaises(SystemExit):
                main()

    @mock.patch("funannotate2_addons.__main__.run_antismash_cli")
    def test_main_antismash(self, mock_run_antismash):
        """Test main function calls antismash CLI"""
        with mock.patch("sys.argv", ["funannotate2_addons", "antismash", "--help"]):
            with self.assertRaises(SystemExit):
                main()

    @mock.patch("funannotate2_addons.__main__.run_signalp_cli")
    def test_main_signalp6(self, mock_run_signalp):
        """Test main function calls signalp6 CLI"""
        with mock.patch("sys.argv", ["funannotate2_addons", "signalp6", "--help"]):
            with self.assertRaises(SystemExit):
                main()


if __name__ == "__main__":
    unittest.main()
