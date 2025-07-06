#!/usr/bin/env python3

import os
import unittest
from unittest import mock
import tempfile
import shutil
import json
import subprocess

from funannotate2_addons.antismash import (
    get_version,
    parse_antismash_gbk,
    antismash_to_json,
)


class TestAntismash(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for test files
        self.test_dir = tempfile.mkdtemp()

        # Create a sample antiSMASH GenBank file
        self.sample_gbk_file = os.path.join(self.test_dir, "sample.gbk")
        with open(self.sample_gbk_file, "w") as f:
            f.write(
                "LOCUS       sample_region001        5000 bp    DNA     linear   UNK 01-JAN-1980\n"
            )
            f.write("DEFINITION  sample_region001\n")
            f.write("ACCESSION   sample_region001\n")
            f.write("VERSION     sample_region001\n")
            f.write("KEYWORDS    .\n")
            f.write("SOURCE      .\n")
            f.write("  ORGANISM  .\n")
            f.write("            .\n")
            f.write("FEATURES             Location/Qualifiers\n")
            f.write("     source          1..5000\n")
            f.write('                     /organism="Test organism"\n')
            f.write("     region          1..5000\n")
            f.write('                     /region_number="1"\n')
            f.write('                     /product="T1PKS"\n')
            f.write('                     /candidate_cluster_number="1"\n')
            f.write(
                '                     /detection_rule="(PKS_KS and (PKS_AT or ene_KS))"\n'
            )
            f.write("     CDS             100..1500\n")
            f.write('                     /gene="gene1"\n')
            f.write('                     /locus_tag="gene1"\n')
            f.write('                     /product="polyketide synthase"\n')
            f.write(
                '                     /sec_met_domain="PKS_KS (E-value: 1.2e-162, bitscore:\n'
            )
            f.write(
                '                     534.7, seeds: 2284, tool: rule-based-clusters)"\n'
            )
            f.write('                     /translation="MTEST..."\n')
            f.write("     CDS             2000..3000\n")
            f.write('                     /gene="gene2"\n')
            f.write('                     /locus_tag="gene2"\n')
            f.write('                     /product="hypothetical protein"\n')
            f.write('                     /translation="MHYPO..."\n')
            f.write("ORIGIN\n")
            f.write(
                "        1 atgcgatcga tcgatcgatc gatcgatcga tcgatcgatc gatcgatcga\n"
            )
            f.write("//\n")

    def tearDown(self):
        # Remove the temporary directory
        shutil.rmtree(self.test_dir)

    @mock.patch("subprocess.run")
    def test_get_version(self, mock_run):
        # Mock the subprocess.run to return a version
        mock_process = mock.Mock()
        mock_process.stdout = "antiSMASH version 7.1.0"
        mock_process.returncode = 0
        mock_run.return_value = mock_process

        # Test the function
        version = get_version("antismash")
        self.assertEqual(version, "7.1.0")

    @mock.patch("subprocess.run")
    def test_get_version_error(self, mock_run):
        # Mock subprocess error
        mock_run.side_effect = subprocess.SubprocessError("Command failed")

        # Test the function
        version = get_version("antismash")
        self.assertIsNone(version)

    def test_parse_antismash_gbk(self):
        # Test parsing the sample GenBank file
        clusters = parse_antismash_gbk(self.sample_gbk_file)

        # Check that clusters were parsed correctly
        self.assertEqual(len(clusters), 1)

        cluster = clusters[0]
        self.assertEqual(cluster["region_number"], "1")
        self.assertEqual(cluster["product"], "T1PKS")
        self.assertEqual(cluster["candidate_cluster_number"], "1")
        self.assertIn("detection_rule", cluster)

        # Check genes in cluster
        self.assertEqual(len(cluster["genes"]), 2)

        gene1 = cluster["genes"][0]
        self.assertEqual(gene1["locus_tag"], "gene1")
        self.assertEqual(gene1["product"], "polyketide synthase")
        self.assertIn("sec_met_domain", gene1)

        gene2 = cluster["genes"][1]
        self.assertEqual(gene2["locus_tag"], "gene2")
        self.assertEqual(gene2["product"], "hypothetical protein")

    def test_antismash_to_json(self):
        # Test converting to JSON
        json_file = os.path.join(self.test_dir, "output.json")
        result = antismash_to_json(self.sample_gbk_file, json_file)

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

        # Check gene1 annotations
        gene1_data = data["gene1"]
        self.assertEqual(gene1_data["product"], "polyketide synthase")
        self.assertIn("antiSMASH", gene1_data["note"])

        # Check gene2 annotations (should not have antiSMASH annotations)
        gene2_data = data["gene2"]
        self.assertEqual(gene2_data["product"], "hypothetical protein")


if __name__ == "__main__":
    unittest.main()
