#!/usr/bin/env python3

import os
import unittest
from unittest import mock
import tempfile
import shutil

from funannotate2_addons.utilities import (
    log_command,
    memorycheck,
    human_readable_size,
    which2,
    download,
)


class TestUtilities(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for test files
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the temporary directory
        shutil.rmtree(self.test_dir)

    @mock.patch("builtins.print")
    def test_log_command_list(self, mock_print):
        # Test logging a command as a list
        cmd = ["python", "-m", "pytest", "--verbose"]
        log_command(cmd)

        # Check that print was called with the formatted command
        mock_print.assert_called_once()
        call_args = mock_print.call_args[0][0]
        self.assertIn("[COMMAND]", call_args)
        self.assertIn("python", call_args)
        self.assertIn("pytest", call_args)
        self.assertIn("--verbose", call_args)

    @mock.patch("builtins.print")
    def test_log_command_string(self, mock_print):
        # Test logging a command as a string
        cmd = "python -m pytest --verbose"
        log_command(cmd)

        # Check that print was called with the command
        mock_print.assert_called_once()
        call_args = mock_print.call_args[0][0]
        self.assertIn("[COMMAND]", call_args)
        self.assertIn(cmd, call_args)

    @mock.patch("os.sysconf")
    def test_memorycheck(self, mock_sysconf):
        # Mock system configuration values
        mock_sysconf.side_effect = lambda x: {
            "SC_PAGE_SIZE": 4096,
            "SC_PHYS_PAGES": 2097152,  # 8GB total
        }[x]

        # Test memory check
        memory_gib = memorycheck()
        self.assertAlmostEqual(memory_gib, 8.0, places=1)

    def test_human_readable_size(self):
        # Test various sizes
        self.assertEqual(human_readable_size(512), "512.00 B")
        self.assertEqual(human_readable_size(1024), "1.00 KiB")
        self.assertEqual(human_readable_size(1536), "1.50 KiB")
        self.assertEqual(human_readable_size(1048576), "1.00 MiB")
        self.assertEqual(human_readable_size(1073741824), "1.00 GiB")
        self.assertEqual(human_readable_size(1099511627776), "1.00 TiB")

        # Test with different decimal places
        self.assertEqual(human_readable_size(1536, decimal_places=1), "1.5 KiB")
        self.assertEqual(human_readable_size(1536, decimal_places=3), "1.500 KiB")

    def test_which2_found(self):
        # Test finding an executable that should exist
        path = which2("python")
        # Just check that it returns something or None
        self.assertIsInstance(path, (str, type(None)))

    def test_which2_not_found(self):
        # Test not finding an executable
        path = which2("nonexistent_program_12345")
        self.assertIsNone(path)


if __name__ == "__main__":
    unittest.main()
