#!/usr/bin/env python3

import os
import unittest
from unittest import mock
import tempfile
import shutil
import json

from funannotate2_addons.signalp import (
    get_version,
    parse_signalp,
    parse_signalp_json,
    signalp_to_json,
)


class TestSignalP(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for test files
        self.test_dir = tempfile.mkdtemp()
        
        # Create a sample SignalP output file
        self.sample_signalp_file = os.path.join(self.test_dir, "sample_summary.signalp6")
        with open(self.sample_signalp_file, "w") as f:
            f.write("# SignalP-6.0 predictions\n")
            f.write("# ID\tPrediction\tPROB\tCS_pos\tProtein length\n")
            f.write("protein1\tSP(Sec/SPI)\t0.9876\t24\t300\n")
            f.write("protein2\tLIPO(Sec/SPII)\t0.8765\t18\t250\n")
            f.write("protein3\tTAT(Tat/SPI)\t0.7654\t30\t400\n")
            f.write("protein4\tOTHER\t0.1234\t-\t200\n")
        
        # Create a sample SignalP JSON output file
        self.sample_json_file = os.path.join(self.test_dir, "sample.json")
        with open(self.sample_json_file, "w") as f:
            json.dump([
                {"ID": "protein1", "PREDICTION": "SP(Sec/SPI)", "PROB": 0.9876, "CS_POS": "24", "PROTEIN_LENGTH": 300},
                {"ID": "protein2", "PREDICTION": "LIPO(Sec/SPII)", "PROB": 0.8765, "CS_POS": "18", "PROTEIN_LENGTH": 250},
                {"ID": "protein3", "PREDICTION": "TAT(Tat/SPI)", "PROB": 0.7654, "CS_POS": "30", "PROTEIN_LENGTH": 400},
                {"ID": "protein4", "PREDICTION": "OTHER", "PROB": 0.1234, "CS_POS": "-", "PROTEIN_LENGTH": 200}
            ], f)
    
    def tearDown(self):
        # Remove the temporary directory and its contents
        shutil.rmtree(self.test_dir)
    
    @mock.patch("subprocess.run")
    def test_get_version(self, mock_run):
        # Mock the subprocess.run to return a version
        mock_process = mock.Mock()
        mock_process.stdout = "SignalP 6.0"
        mock_run.return_value = mock_process
        
        # Test the function
        version = get_version("signalp")
        self.assertEqual(version, "6.0")
        
        # Test with different output format
        mock_process.stdout = "This is SignalP version 6.0"
        version = get_version("signalp")
        self.assertEqual(version, "6.0")
    
    def test_parse_signalp(self):
        # Test parsing the sample file
        predictions = parse_signalp(self.sample_signalp_file)
        
        # Check that predictions were parsed correctly
        self.assertIn("protein1", predictions)
        self.assertEqual(predictions["protein1"]["prediction"], "SP(Sec/SPI)")
        self.assertEqual(predictions["protein1"]["probability"], 0.9876)
        self.assertEqual(predictions["protein1"]["cleavage_site"], "24")
        self.assertTrue(predictions["protein1"]["has_signal_peptide"])
        
        self.assertIn("protein2", predictions)
        self.assertEqual(predictions["protein2"]["prediction"], "LIPO(Sec/SPII)")
        self.assertTrue(predictions["protein2"]["has_signal_peptide"])
        
        self.assertIn("protein3", predictions)
        self.assertEqual(predictions["protein3"]["prediction"], "TAT(Tat/SPI)")
        self.assertTrue(predictions["protein3"]["has_signal_peptide"])
        
        self.assertIn("protein4", predictions)
        self.assertEqual(predictions["protein4"]["prediction"], "OTHER")
        self.assertFalse(predictions["protein4"]["has_signal_peptide"])
        
        # Test writing to output file
        output_file = os.path.join(self.test_dir, "output.txt")
        parse_signalp(self.sample_signalp_file, output_file=output_file)
        
        # Check that the output file was created
        self.assertTrue(os.path.isfile(output_file))
        
        # Read the output file and check its contents
        with open(output_file, "r") as f:
            lines = f.readlines()
        
        # Check header
        self.assertEqual(lines[0].strip(), "#gene_id\tannotation_type\tannotation_value")
        
        # Check that only proteins with signal peptides are included
        protein_lines = [line.split("\t")[0] for line in lines[1:]]
        self.assertIn("protein1", protein_lines)
        self.assertIn("protein2", protein_lines)
        self.assertIn("protein3", protein_lines)
        self.assertNotIn("protein4", protein_lines)
        
        # Check that protein1 has the correct annotation
        protein1_line = next(line for line in lines if line.startswith("protein1"))
        self.assertIn("SignalP:SP(Sec/SPI)", protein1_line)
        self.assertIn("prob=0.988", protein1_line)
        self.assertIn("cleavage_site=24", protein1_line)
    
    def test_parse_signalp_json(self):
        # Test parsing the sample JSON file
        predictions = parse_signalp_json(self.sample_json_file)
        
        # Check that predictions were parsed correctly
        self.assertIn("protein1", predictions)
        self.assertEqual(predictions["protein1"]["prediction"], "SP(Sec/SPI)")
        self.assertEqual(predictions["protein1"]["probability"], 0.9876)
        self.assertEqual(predictions["protein1"]["cleavage_site"], "24")
        self.assertTrue(predictions["protein1"]["has_signal_peptide"])
        
        self.assertIn("protein2", predictions)
        self.assertEqual(predictions["protein2"]["prediction"], "LIPO(Sec/SPII)")
        self.assertTrue(predictions["protein2"]["has_signal_peptide"])
        
        self.assertIn("protein3", predictions)
        self.assertEqual(predictions["protein3"]["prediction"], "TAT(Tat/SPI)")
        self.assertTrue(predictions["protein3"]["has_signal_peptide"])
        
        self.assertIn("protein4", predictions)
        self.assertEqual(predictions["protein4"]["prediction"], "OTHER")
        self.assertFalse(predictions["protein4"]["has_signal_peptide"])
        
        # Test writing to output file
        output_file = os.path.join(self.test_dir, "output.txt")
        parse_signalp_json(self.sample_json_file, output_file=output_file)
        
        # Check that the output file was created
        self.assertTrue(os.path.isfile(output_file))
    
    def test_signalp_to_json(self):
        # Test converting standard format to JSON
        json_file = os.path.join(self.test_dir, "output.json")
        result = signalp_to_json(self.sample_signalp_file, json_file)
        
        # Check that the function returned True
        self.assertTrue(result)
        
        # Check that the JSON file was created
        self.assertTrue(os.path.isfile(json_file))
        
        # Read the JSON file and check its contents
        with open(json_file, "r") as f:
            data = json.load(f)
        
        self.assertIn("protein1", data)
        self.assertEqual(data["protein1"]["prediction"], "SP(Sec/SPI)")
        self.assertEqual(data["protein1"]["probability"], 0.9876)
        
        # Test with JSON input
        json_file2 = os.path.join(self.test_dir, "output2.json")
        result = signalp_to_json(self.sample_json_file, json_file2)
        
        # Check that the function returned True
        self.assertTrue(result)
        
        # Check that the JSON file was created
        self.assertTrue(os.path.isfile(json_file2))


if __name__ == "__main__":
    unittest.main()
