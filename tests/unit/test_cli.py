#!/usr/bin/env python3

import argparse
import unittest
import sys
from unittest import mock
from io import StringIO

from funannotate2_addons.__main__ import __version__, parse_args, main


class TestCLI(unittest.TestCase):
    def _assert_main_dispatch(self, subparser_name, patch_target):
        args = argparse.Namespace(subparser_name=subparser_name)

        with mock.patch("funannotate2_addons.__main__.parse_args", return_value=args):
            with mock.patch(patch_target) as mock_runner:
                main()

        mock_runner.assert_called_once_with(args)

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
        """Test that --version prints the current package version"""
        with self.assertRaises(SystemExit) as cm:
            with mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
                parse_args(["--version"])
        self.assertEqual(cm.exception.code, 0)
        self.assertIn("funannotate2_addons", mock_stdout.getvalue())
        self.assertIn(__version__, mock_stdout.getvalue())

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
        self._assert_main_dispatch("emapper", "funannotate2_addons.__main__.run_emapper_cli")

    @mock.patch("funannotate2_addons.__main__.run_iprscan_cli")
    def test_main_iprscan(self, mock_run_iprscan):
        """Test main function calls iprscan CLI"""
        self._assert_main_dispatch("iprscan", "funannotate2_addons.__main__.run_iprscan_cli")

    @mock.patch("funannotate2_addons.__main__.run_antismash_cli")
    def test_main_antismash(self, mock_run_antismash):
        """Test main function calls antismash CLI"""
        self._assert_main_dispatch(
            "antismash", "funannotate2_addons.__main__.run_antismash_cli"
        )

    @mock.patch("funannotate2_addons.__main__.run_signalp_cli")
    def test_main_signalp6(self, mock_run_signalp):
        """Test main function calls signalp6 CLI"""
        self._assert_main_dispatch("signalp6", "funannotate2_addons.__main__.run_signalp_cli")


if __name__ == "__main__":
    unittest.main()
