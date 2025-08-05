#!/usr/bin/env python3

import os
import unittest
from unittest import mock
import tempfile
import shutil
import json
import xml.etree.ElementTree as ET
import subprocess

from funannotate2_addons.iprscan import (
    get_version,
    parse_iprscan_tsv,
    parse_iprscan_xml,
    iprscan_to_json,
)


class TestIprscan(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for test files
        self.test_dir = tempfile.mkdtemp()

        # Create a sample InterProScan TSV file
        self.sample_tsv_file = os.path.join(self.test_dir, "sample.tsv")
        with open(self.sample_tsv_file, "w") as f:
            f.write(
                "protein1\tMD5\t100\tPfam\tPF00001\tDNA binding domain\t10\t90\t1.2E-15\tT\t01-01-2023\tIPR000001\tDNA binding\tGO:0003677\n"
            )
            f.write(
                "protein1\tMD5\t100\tSMART\tSM00001\tHelix-turn-helix\t15\t85\t2.3E-10\tT\t01-01-2023\tIPR000001\tDNA binding\tGO:0003677\n"
            )
            f.write(
                "protein2\tMD5\t200\tPfam\tPF00002\tKinase domain\t20\t180\t5.4E-20\tT\t01-01-2023\tIPR000002\tProtein kinase\tGO:0004672\n"
            )
            f.write(
                "protein3\tMD5\t150\tPfam\tPF00003\tUnknown function\t5\t145\t1.0E-5\tT\t01-01-2023\t-\t-\t-\n"
            )

        # Create a sample InterProScan XML file
        self.sample_xml_file = os.path.join(self.test_dir, "sample.xml")
        xml_content = """<?xml version="1.0" encoding="UTF-8"?>
<protein-matches xmlns="http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5">
    <protein>
        <sequence md5="abc123">protein1</sequence>
        <matches>
            <pfam-match>
                <signature ac="PF00001" name="DNA_binding">
                    <entry ac="IPR000001" desc="DNA binding" name="DNA_binding" type="Domain"/>
                    <locations>
                        <pfam-location start="10" end="90" score="15.2" evalue="1.2E-15"/>
                    </locations>
                    <go-xrefs>
                        <go-xref id="GO:0003677" name="DNA binding"/>
                    </go-xrefs>
                </signature>
            </pfam-match>
        </matches>
    </protein>
    <protein>
        <sequence md5="def456">protein2</sequence>
        <matches>
            <pfam-match>
                <signature ac="PF00002" name="Kinase">
                    <entry ac="IPR000002" desc="Protein kinase" name="Kinase" type="Domain"/>
                    <locations>
                        <pfam-location start="20" end="180" score="20.4" evalue="5.4E-20"/>
                    </locations>
                    <go-xrefs>
                        <go-xref id="GO:0004672" name="protein kinase activity"/>
                    </go-xrefs>
                </signature>
            </pfam-match>
        </matches>
    </protein>
</protein-matches>"""
        with open(self.sample_xml_file, "w") as f:
            f.write(xml_content)

    def tearDown(self):
        # Remove the temporary directory
        shutil.rmtree(self.test_dir)

    @mock.patch("subprocess.run")
    def test_get_version(self, mock_run):
        # Mock the subprocess.run to return a version
        mock_process = mock.Mock()
        mock_process.stdout = "InterProScan version 5.61-93.0"
        mock_process.returncode = 0
        mock_run.return_value = mock_process

        # Test the function
        version = get_version("interproscan.sh")
        self.assertEqual(version, "5.61-93.0")

    @mock.patch("subprocess.run")
    def test_get_version_error(self, mock_run):
        # Mock subprocess error
        mock_run.side_effect = subprocess.SubprocessError("Command failed")

        # Test the function
        version = get_version("interproscan.sh")
        self.assertIsNone(version)

    def test_parse_iprscan_tsv(self):
        # Test parsing the sample TSV file
        annotations = parse_iprscan_tsv(self.sample_tsv_file)

        # Check that annotations were parsed correctly
        self.assertIn("protein1", annotations)
        self.assertIn("protein2", annotations)
        self.assertIn("protein3", annotations)

        # Check protein1 annotations
        protein1 = annotations["protein1"]
        self.assertEqual(len(protein1["signatures"]), 2)  # PF00001 and SM00001
        self.assertEqual(len(protein1["interpro_domains"]), 1)
        self.assertEqual(len(protein1["go_terms"]), 1)

        # Check signatures (Pfam and SMART)
        signature_ids = [s["id"] for s in protein1["signatures"]]
        self.assertIn("PF00001", signature_ids)
        self.assertIn("SM00001", signature_ids)

        # Check InterPro domain
        self.assertEqual(protein1["interpro_domains"][0]["id"], "IPR000001")

        # Check GO terms
        self.assertEqual(protein1["go_terms"][0]["id"], "GO:0003677")

        # Check protein2 annotations
        protein2 = annotations["protein2"]
        self.assertEqual(len(protein2["signatures"]), 1)
        signature_ids = [s["id"] for s in protein2["signatures"]]
        self.assertIn("PF00002", signature_ids)
        self.assertEqual(protein2["go_terms"][0]["id"], "GO:0004672")

        # Check protein3 annotations (no InterPro or GO terms)
        protein3 = annotations["protein3"]
        self.assertEqual(len(protein3["signatures"]), 1)
        self.assertEqual(len(protein3["interpro_domains"]), 0)
        self.assertEqual(len(protein3["go_terms"]), 0)
        signature_ids = [s["id"] for s in protein3["signatures"]]
        self.assertIn("PF00003", signature_ids)

    def test_parse_iprscan_xml(self):
        # Test parsing the sample XML file
        annotations = parse_iprscan_xml(self.sample_xml_file)

        # Check that annotations were parsed correctly
        self.assertIsInstance(annotations, dict)
        self.assertTrue(len(annotations) > 0)

        # Check that we have protein entries
        protein_ids = list(annotations.keys())
        self.assertTrue(len(protein_ids) > 0)

        # Check the structure of the first protein
        first_protein = annotations[protein_ids[0]]
        self.assertIn("signatures", first_protein)
        self.assertIn("interpro_domains", first_protein)
        self.assertIn("go_terms", first_protein)
        self.assertIn("pathways", first_protein)

        # Check that signatures are present
        self.assertIsInstance(first_protein["signatures"], list)
        self.assertIsInstance(first_protein["interpro_domains"], list)
        self.assertIsInstance(first_protein["go_terms"], list)
        self.assertIsInstance(first_protein["pathways"], list)

    def test_iprscan_to_json(self):
        # Test converting TSV to JSON
        json_file = os.path.join(self.test_dir, "output.json")
        result = iprscan_to_json(self.sample_tsv_file, json_file)

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

        # Check protein1 data
        protein1_data = data["protein1"]
        self.assertIn("signatures", protein1_data)
        self.assertIn("interpro_domains", protein1_data)
        self.assertIn("go_terms", protein1_data)


if __name__ == "__main__":
    unittest.main()
