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
                "LOCUS       sample_region001          50 bp    DNA     linear   UNK 01-JAN-1980\n"
            )
            f.write("DEFINITION  sample_region001\n")
            f.write("ACCESSION   sample_region001\n")
            f.write("VERSION     sample_region001\n")
            f.write("KEYWORDS    .\n")
            f.write("SOURCE      .\n")
            f.write("  ORGANISM  .\n")
            f.write("            .\n")
            f.write("FEATURES             Location/Qualifiers\n")
            f.write("     source          1..50\n")
            f.write('                     /organism="Test organism"\n')
            f.write("     region          1..50\n")
            f.write('                     /region_number="1"\n')
            f.write('                     /product="T1PKS"\n')
            f.write('                     /candidate_cluster_number="1"\n')
            f.write(
                '                     /detection_rule="(PKS_KS and (PKS_AT or ene_KS))"\n'
            )
            f.write("     CDS             10..30\n")
            f.write('                     /gene="gene1"\n')
            f.write('                     /locus_tag="gene1"\n')
            f.write('                     /product="polyketide synthase"\n')
            f.write(
                '                     /sec_met_domain="PKS_KS (E-value: 1.2e-162, bitscore:\n'
            )
            f.write(
                '                     534.7, seeds: 2284, tool: rule-based-clusters)"\n'
            )
            f.write('                     /translation="MTEST"\n')
            f.write("     CDS             35..45\n")
            f.write('                     /gene="gene2"\n')
            f.write('                     /locus_tag="gene2"\n')
            f.write('                     /product="hypothetical protein"\n')
            f.write('                     /translation="MHYPO"\n')
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
        backbone_domains, backbone_subtype, backbone_enzymes, cluster_genes = (
            parse_antismash_gbk(self.sample_gbk_file)
        )

        # Check that the function returned valid data structures
        self.assertIsInstance(backbone_domains, dict)
        self.assertIsInstance(backbone_subtype, dict)
        self.assertIsInstance(backbone_enzymes, dict)
        self.assertIsInstance(cluster_genes, dict)

        # Check that gene1 has domain annotations
        self.assertIn("gene1", backbone_domains)
        gene1_domains = backbone_domains["gene1"]
        self.assertIsInstance(gene1_domains, list)
        self.assertTrue(len(gene1_domains) > 0)

        # Check that the domain annotation contains PKS_KS
        domain_annotation = gene1_domains[0]
        self.assertIn("PKS_KS", domain_annotation)
        self.assertIn("E-value: 1.2e-162", domain_annotation)

        # Check cluster genes structure
        self.assertIn("cluster_1", cluster_genes)
        cluster_1_genes = cluster_genes["cluster_1"]
        self.assertIn("gene1", cluster_1_genes)
        self.assertIn("gene2", cluster_1_genes)

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

        # Check that the JSON data has the expected structure
        self.assertIn("clusters", data)
        self.assertIn("backbone_domains", data)
        self.assertIn("backbone_enzymes", data)

        # Check that backbone_domains contains gene1
        backbone_domains = data["backbone_domains"]
        self.assertIn("gene1", backbone_domains)

        # Check that gene1 has domain annotations
        gene1_domains = backbone_domains["gene1"]
        self.assertIsInstance(gene1_domains, list)
        self.assertTrue(len(gene1_domains) > 0)

        # Check cluster data structure
        clusters = data["clusters"]
        self.assertIsInstance(clusters, dict)


if __name__ == "__main__":
    unittest.main()
