#!/usr/bin/env python3

import os
import unittest
from unittest import mock
import tempfile
import shutil
import json

from funannotate2_addons.emapper import (
    get_version,
    parse_emapper_annotations,
    emapper_to_json,
)


class TestEmapper(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for test files
        self.test_dir = tempfile.mkdtemp()
        
        # Create a sample eggnog-mapper output file
        self.sample_emapper_file = os.path.join(self.test_dir, "sample.emapper.annotations")
        with open(self.sample_emapper_file, "w") as f:
            f.write("##emapper-2.1.9\n")
            f.write("#query_name\tseed_ortholog\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\ttaxonomic scope\tOG\tbestOG\tCOG\tAnnotation\n")
            f.write("gene1\t1034452.ENOG4105E7Z\t1034452.ENOG4105E7Z|1|1e-40,2759.ENOG41052UI|1|8e-40\t1034452\tJ\tRNA helicase\tgeneA\tGO:0003723,GO:0005524\t3.6.4.13\tK03257\tko03040\t\t\t\t\t\t\t\t1034452\t1034452.ENOG4105E7Z\t1034452.ENOG4105E7Z\tCOG0513\tRNA helicase\n")
            f.write("gene2\t1034452.ENOG4105QR3\t1034452.ENOG4105QR3|1|1e-30,2759.ENOG41052XC|1|6e-30\t1034452\tK\tTranscription factor\tgeneB\tGO:0003700,GO:0006355\t\tK09095\tko03022\t\t\t\t\t\t\t\t1034452\t1034452.ENOG4105QR3\t1034452.ENOG4105QR3\tCOG0789\tTranscription factor\n")
            f.write("gene3\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n")
        
        # Create a sample annotations file
        self.sample_annotations_file = os.path.join(self.test_dir, "sample.emapper.annotations.txt")
        
    def tearDown(self):
        # Remove the temporary directory
        shutil.rmtree(self.test_dir)
    
    @mock.patch("subprocess.run")
    def test_get_version(self, mock_run):
        # Mock the subprocess.run to return a version
        mock_process = mock.Mock()
        mock_process.stdout = "emapper-2.1.9"
        mock_process.returncode = 0
        mock_run.return_value = mock_process
        
        # Test the function
        version = get_version("emapper.py")
        self.assertEqual(version, "2.1.9")
    
    def test_parse_emapper_annotations(self):
        # Test parsing the sample file
        annotations = parse_emapper_annotations(self.sample_emapper_file)
        
        # Check that annotations were parsed correctly
        self.assertIn("gene1", annotations)
        self.assertIn("gene2", annotations)
        
        # Check specific annotations
        self.assertEqual(annotations["gene1"]["name"], ["geneA"])
        self.assertEqual(annotations["gene1"]["product"], ["RNA helicase"])
        self.assertEqual(annotations["gene1"]["ec_number"], ["3.6.4.13"])
        self.assertEqual(annotations["gene1"]["go_terms"], ["GO:0003723", "GO:0005524"])
        
        self.assertEqual(annotations["gene2"]["name"], ["geneB"])
        self.assertEqual(annotations["gene2"]["product"], ["Transcription factor"])
        self.assertEqual(annotations["gene2"]["go_terms"], ["GO:0003700", "GO:0006355"])
    
    def test_emapper_to_json(self):
        # Test converting to JSON
        json_file = os.path.join(self.test_dir, "output.json")
        result = emapper_to_json(self.sample_emapper_file, json_file)
        
        # Check that the function returned True
        self.assertTrue(result)
        
        # Check that the JSON file was created
        self.assertTrue(os.path.exists(json_file))
        
        # Check the content of the JSON file
        with open(json_file, "r") as f:
            data = json.load(f)
        
        # Check that the JSON data is correct
        self.assertIn("gene1", data)
        self.assertIn("gene2", data)
        self.assertEqual(data["gene1"]["name"], "geneA")
        self.assertEqual(data["gene1"]["product"], "RNA helicase")
        self.assertEqual(data["gene1"]["ec_number"], "3.6.4.13")
        self.assertEqual(data["gene2"]["name"], "geneB")
        self.assertEqual(data["gene2"]["product"], "Transcription factor")


if __name__ == "__main__":
    unittest.main()
