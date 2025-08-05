#!/usr/bin/env python3

import unittest
import tempfile
import shutil
import os
import json
import subprocess
from unittest import mock

from funannotate2_addons.signalp6 import (
    get_version,
    parse_signalp,
    signalp_to_json,
    run_signalp,
)


class TestSignalp6(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for test files
        self.test_dir = tempfile.mkdtemp()

        # Create a sample SignalP6 prediction_results.txt file
        self.sample_results_file = os.path.join(self.test_dir, "prediction_results.txt")
        with open(self.sample_results_file, "w") as f:
            f.write("# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 20250804212443\n")
            f.write("# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n")
            f.write("protein1\tOTHER\t1.000053\t0.000000\t\n")
            f.write("protein2\tSP\t0.000229\t0.999747\tCS pos: 17-18. Pr: 0.9751\n")
            f.write("protein3\tSP\t0.007633\t0.992322\tCS pos: 23-24. Pr: 0.6512\n")
            f.write("protein4\tOTHER\t1.000045\t0.000000\t\n")

        # Create a sample protein FASTA file
        self.sample_fasta_file = os.path.join(self.test_dir, "proteins.fasta")
        with open(self.sample_fasta_file, "w") as f:
            f.write(">protein1\n")
            f.write("MHACMHACMHACMHACMHAC\n")
            f.write(">protein2\n")
            f.write(
                "MKLLVVLSLVLAFAVAFPQFVAQIHENLKKQGNFGEVFCLFEGSQNMQEQNKQVQYLQEQKLISEEDL\n"
            )
            f.write(">protein3\n")
            f.write(
                "MKLKWSLLLVLGLVLPVLGAAEVQVTGPGPQPQPQPQPQPQPQPQPQPQPQPQPQPQPQPQPQPQPQ\n"
            )
            f.write(">protein4\n")
            f.write("VRTVRTVRTVRTVRTVRT\n")

    def tearDown(self):
        # Remove the temporary directory
        shutil.rmtree(self.test_dir)

    @mock.patch("subprocess.run")
    def test_get_version(self, mock_run):
        # Mock the subprocess.run to return a version
        mock_process = mock.Mock()
        mock_process.stdout = "SignalP version 6.0g"
        mock_process.returncode = 0
        mock_run.return_value = mock_process

        # Test the function
        version = get_version("signalp6")
        self.assertEqual(version, "6.0g")

    @mock.patch("subprocess.run")
    def test_get_version_error(self, mock_run):
        # Mock the subprocess.run to raise an exception
        mock_run.side_effect = subprocess.SubprocessError("Command failed")

        # Test the function
        version = get_version("signalp6")
        self.assertIsNone(version)

    def test_parse_signalp(self):
        # Test parsing the sample prediction results file
        predictions = parse_signalp(self.sample_results_file)

        # Check that predictions were parsed correctly
        self.assertEqual(len(predictions), 4)
        self.assertIn("protein1", predictions)
        self.assertIn("protein2", predictions)
        self.assertIn("protein3", predictions)
        self.assertIn("protein4", predictions)

        # Check protein1 (no signal peptide)
        protein1 = predictions["protein1"]
        self.assertEqual(protein1["prediction"], "OTHER")
        self.assertFalse(protein1["has_signal_peptide"])
        self.assertEqual(protein1["probability"], 1.000053)

        # Check protein2 (has signal peptide with cleavage site)
        protein2 = predictions["protein2"]
        self.assertEqual(protein2["prediction"], "SP")
        self.assertTrue(protein2["has_signal_peptide"])
        self.assertEqual(protein2["probability"], 0.999747)
        self.assertEqual(protein2["cleavage_site"], "17-18")
        self.assertEqual(protein2["cleavage_prob"], 0.9751)

        # Check protein3 (has signal peptide with cleavage site)
        protein3 = predictions["protein3"]
        self.assertEqual(protein3["prediction"], "SP")
        self.assertTrue(protein3["has_signal_peptide"])
        self.assertEqual(protein3["cleavage_site"], "23-24")
        self.assertEqual(protein3["cleavage_prob"], 0.6512)

    def test_parse_signalp_with_output_file(self):
        # Test parsing with output file generation
        output_file = os.path.join(self.test_dir, "signalp.annotations.txt")
        predictions = parse_signalp(self.sample_results_file, output_file=output_file)

        # Check that the output file was created
        self.assertTrue(os.path.exists(output_file))

        # Check the content of the output file
        with open(output_file, "r") as f:
            lines = f.readlines()

        # Should have header + 2 signal peptide predictions
        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[0].startswith("#gene_id"))

        # Check annotation format
        self.assertIn("protein2\tnote\tSignalP:1-17", lines[1])
        self.assertIn("protein3\tnote\tSignalP:1-23", lines[2])

    def test_signalp_to_json(self):
        # Test converting to JSON
        json_file = os.path.join(self.test_dir, "output.json")
        result = signalp_to_json(self.sample_results_file, json_file)

        # Check that the function returned True
        self.assertTrue(result)

        # Check that the JSON file was created
        self.assertTrue(os.path.exists(json_file))

        # Check the content of the JSON file
        with open(json_file, "r") as f:
            data = json.load(f)

        # Check that the JSON data is correct
        self.assertIn("protein1", data)
        self.assertIn("protein2", data)
        self.assertIn("protein3", data)
        self.assertIn("protein4", data)

        # Check protein2 data (has signal peptide)
        protein2_data = data["protein2"]
        self.assertEqual(protein2_data["prediction"], "SP")
        self.assertTrue(protein2_data["has_signal_peptide"])
        self.assertEqual(protein2_data["cleavage_site"], "17-18")

    @mock.patch("subprocess.run")
    @mock.patch("funannotate2_addons.signalp6.log_command")
    def test_run_signalp(self, mock_log_command, mock_run):
        # Mock successful subprocess run
        mock_process = mock.Mock()
        mock_process.returncode = 0
        mock_run.return_value = mock_process

        # Create output directory
        output_dir = os.path.join(self.test_dir, "signalp_output")
        os.makedirs(output_dir)

        # Create expected output file
        expected_output = os.path.join(output_dir, "prediction_results.txt")
        with open(expected_output, "w") as f:
            f.write("# SignalP-6.0\tOrganism: Eukarya\n")
            f.write("# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n")
            f.write("protein1\tOTHER\t1.0\t0.0\t\n")

        # Test the function
        result = run_signalp(
            input_file=self.sample_fasta_file,
            output_dir=output_dir,
            signalp_path="signalp6",
        )

        # Check that the function returned the expected output file
        self.assertEqual(result, expected_output)

        # Check that subprocess.run was called with correct arguments
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]  # Get the command list
        self.assertIn("signalp6", call_args)
        self.assertIn("--fastafile", call_args)
        self.assertIn(self.sample_fasta_file, call_args)
        self.assertIn("--output_dir", call_args)
        self.assertIn(output_dir, call_args)


if __name__ == "__main__":
    unittest.main()
