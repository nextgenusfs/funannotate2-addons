#!/usr/bin/env python3

import argparse
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
    run_emapper_cli,
)


class TestEmapper(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for test files
        self.test_dir = tempfile.mkdtemp()

        # Create a sample eggnog-mapper output file
        self.sample_emapper_file = os.path.join(
            self.test_dir, "sample.emapper.annotations"
        )
        with open(self.sample_emapper_file, "w") as f:
            f.write("##emapper-2.1.9\n")
            f.write(
                "#query\tseed_ortholog\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\ttaxonomic scope\tOG\tbestOG\tCOG\tAnnotation\n"
            )
            f.write(
                "gene1\t1034452.ENOG4105E7Z\t1034452.ENOG4105E7Z|1|1e-40,2759.ENOG41052UI|1|8e-40\t1034452\tJ\tRNA helicase\tgeneA\tGO:0003723,GO:0005524\t3.6.4.13\tK03257\tko03040\t\t\t\t\t\t\t\t1034452\t1034452.ENOG4105E7Z\t1034452.ENOG4105E7Z\tCOG0513\tRNA helicase\n"
            )
            f.write(
                "gene2\t1034452.ENOG4105QR3\t1034452.ENOG4105QR3|1|1e-30,2759.ENOG41052XC|1|6e-30\t1034452\tK\tTranscription factor\tgeneB\tGO:0003700,GO:0006355\t\tK09095\tko03022\t\t\t\t\t\t\t\t1034452\t1034452.ENOG4105QR3\t1034452.ENOG4105QR3\tCOG0789\tTranscription factor\n"
            )
            f.write(
                "gene3\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n"
            )

        # Create a sample annotations file
        self.sample_annotations_file = os.path.join(
            self.test_dir, "sample.emapper.annotations.txt"
        )

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
        self.assertEqual(annotations["gene1"]["gene_name"], "geneA")
        self.assertEqual(annotations["gene1"]["description"], "RNA helicase")
        self.assertEqual(annotations["gene1"]["ec_number"], "3.6.4.13")

        self.assertEqual(annotations["gene2"]["gene_name"], "geneB")
        self.assertEqual(annotations["gene2"]["description"], "Transcription factor")

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
        self.assertEqual(data["gene1"]["gene_name"], "geneA")
        self.assertEqual(data["gene1"]["description"], "RNA helicase")
        self.assertEqual(data["gene1"]["ec_number"], "3.6.4.13")
        self.assertEqual(data["gene2"]["gene_name"], "geneB")
        self.assertEqual(data["gene2"]["description"], "Transcription factor")

    @mock.patch("funannotate2_addons.emapper.logger.error")
    @mock.patch("funannotate2_addons.emapper.get_input_file", return_value=None)
    def test_run_emapper_cli_missing_input_file(self, mock_get_input_file, mock_logger):
        args = mock.Mock(input="predict_dir", file=None, parse=None)

        run_emapper_cli(args)

        mock_get_input_file.assert_called_once_with(args, "proteins")
        mock_logger.assert_called_once_with("No protein FASTA file found")

    @mock.patch("funannotate2_addons.emapper.startLogging")
    @mock.patch("funannotate2_addons.emapper.parse_emapper_annotations")
    @mock.patch("funannotate2_addons.emapper.run_emapper")
    @mock.patch("funannotate2_addons.emapper.get_input_file")
    def test_run_emapper_cli_uses_stable_output_prefix(
        self,
        mock_get_input_file,
        mock_run_emapper,
        mock_parse_annotations,
        mock_start_logging,
    ):
        input_file = os.path.join(self.test_dir, "proteins.fa")
        with open(input_file, "w") as f:
            f.write(">gene1\nMSTNPKPQRIT\n")

        output_dir = os.path.join(self.test_dir, "output")
        expected_prefix = os.path.join(output_dir, "emapper")
        mock_get_input_file.return_value = input_file
        mock_run_emapper.return_value = f"{expected_prefix}.emapper.annotations"
        mock_parse_annotations.return_value = {"gene1": {"gene_name": "geneA"}}
        mock_start_logging.return_value = mock.Mock()

        args = argparse.Namespace(
            input="predict_dir",
            file=None,
            parse=None,
            output=output_dir,
            emapper_path="emapper.py",
            cpus=4,
            database=None,
            data_dir=None,
            temp_dir=None,
            resume=True,
            override=False,
            mode="diamond",
            tax_scope="auto",
            target_orthologs="all",
            sensmode="sensitive",
            no_annot=False,
            no_search=False,
            dmnd_db=None,
            seed_orthologs=None,
            report_orthologs=False,
            go_evidence=None,
            pfam_realign=False,
            seed_ortholog_score=None,
            seed_ortholog_evalue=None,
            query_cover=None,
            subject_cover=None,
            matrix=None,
            gapopen=None,
            gapextend=None,
            evalue=None,
            pident=None,
            query_sgcov=None,
            subject_sgcov=None,
            json=False,
        )

        run_emapper_cli(args)

        self.assertTrue(os.path.isdir(output_dir))
        self.assertEqual(mock_run_emapper.call_args.kwargs["output_prefix"], expected_prefix)
        self.assertTrue(mock_run_emapper.call_args.kwargs["resume"])
        mock_parse_annotations.assert_called_once_with(
            f"{expected_prefix}.emapper.annotations",
            output_file=f"{expected_prefix}.annotations.txt",
        )


if __name__ == "__main__":
    unittest.main()
