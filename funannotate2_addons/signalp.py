#!/usr/bin/env python3

import os
import subprocess
import argparse
import logging
import json
import csv
from pathlib import Path
from Bio import SeqIO
from .log import startLogging

# Set up logger
logger = logging.getLogger(__name__)


def get_version(signalp_path):
    """
    Get the version of SignalP.

    Args:
        signalp_path (str): Path to the signalp executable

    Returns:
        str: Version of SignalP
    """
    try:
        cmd = [signalp_path, "--version"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Parse version from output
        for line in proc.stdout.split("\n"):
            if "SignalP" in line and "version" in line:
                return line.split("version")[1].strip()
        return None
    except subprocess.SubprocessError as e:
        logger.error(f"Error getting SignalP version: {e}")
        return None


def run_signalp(
    input_file,
    output_prefix,
    signalp_path="signalp6",
    organism="euk",
    format="short",
    plot=False,
    output_dir=None,
    force=False,
):
    """
    Run SignalP 6.0 on a protein FASTA file.

    Args:
        input_file (str): Path to input protein FASTA file
        output_prefix (str): Prefix for output files
        signalp_path (str): Path to signalp6 executable
        organism (str): Organism group (euk, arch, gram+, gram-)
        format (str): Prediction format (short, long)
        plot (bool): Generate plots
        output_dir (str): Directory for output files
        force (bool): Force overwrite of existing files

    Returns:
        str: Path to the output file
    """
    # Check if input file exists
    if not os.path.isfile(input_file):
        logger.error(f"Input file not found: {input_file}")
        return None

    # Create output directory if specified and it doesn't exist
    if output_dir and not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Build command
    cmd = [
        signalp_path,
        "--fastafile",
        input_file,
        "--prefix",
        output_prefix,
    ]

    # Add optional arguments
    if organism:
        cmd.extend(["--organism", organism])
    if format == "long":
        cmd.append("--verbose")
    if plot:
        cmd.append("--plot")
    if output_dir:
        cmd.extend(["--output", output_dir])
    if force:
        cmd.append("--force")

    # Run command
    logger.info(f"Running SignalP: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)

        # Determine output file based on format
        if format == "short":
            output_file = f"{output_prefix}_summary.signalp6"  # SignalP 6 format
        else:
            output_file = f"{output_prefix}_summary.signalp6"  # SignalP 6 format

        if os.path.isfile(output_file):
            logger.info(f"SignalP completed successfully: {output_file}")
            return output_file
        else:
            logger.error(f"SignalP output file not found: {output_file}")
            return None
    except subprocess.SubprocessError as e:
        logger.error(f"Error running SignalP: {e}")
        return None


def parse_signalp(input_file, output_file=None, gene_dict=None):
    """
    Parse SignalP 6.0 output and convert to a standardized format.

    Args:
        input_file (str): Path to SignalP 6.0 output file
        output_file (str, optional): Path to output file
        gene_dict (dict, optional): Dictionary of gene information

    Returns:
        dict: Dictionary of signal peptide predictions keyed by protein ID
    """
    # Check if input file exists
    if not os.path.isfile(input_file):
        logger.error(f"Input file not found: {input_file}")
        return None

    # Parse SignalP output
    predictions = {}

    try:
        with open(input_file, "r") as f:
            # Skip header lines
            header_seen = False
            for line in f:
                if line.startswith("#"):
                    if "ID" in line and "PREDICTION" in line and "PROB" in line:
                        header_seen = True
                    continue

                # Skip lines until we see the header
                if not header_seen and not line.strip():
                    continue

                # Parse prediction line
                fields = line.strip().split()

                # Skip if not enough fields
                if len(fields) < 3:
                    continue

                # Get values
                protein_id = fields[0]
                prediction = fields[
                    1
                ]  # SP(Sec/SPI), LIPO(Sec/SPII), TAT(Tat/SPI), TATLIPO(Tat/SPII), OTHER
                probability = float(fields[2])

                # Store prediction
                predictions[protein_id] = {
                    "prediction": prediction,
                    "probability": probability,
                    "has_signal_peptide": prediction != "OTHER",
                }

                # Add cleavage site if available
                if len(fields) > 3 and fields[3] != "-":
                    predictions[protein_id]["cleavage_site"] = fields[3]

        # Write to output file if specified
        if output_file:
            with open(output_file, "w") as out:
                # Write header
                out.write("#gene_id\tannotation_type\tannotation_value\n")

                for protein_id, pred in predictions.items():
                    if pred["has_signal_peptide"]:
                        out.write(
                            f"{protein_id}\tnote\tSignalP:{pred['prediction']} prob={pred['probability']:.3f}"
                        )
                        if "cleavage_site" in pred:
                            out.write(f" cleavage_site={pred['cleavage_site']}")
                        out.write("\n")

        return predictions

    except Exception as e:
        logger.error(f"Error parsing SignalP output: {e}")
        return None


def parse_signalp_json(input_file, output_file=None, gene_dict=None):
    """
    Parse SignalP JSON output and convert to a standardized format.

    Args:
        input_file (str): Path to SignalP JSON output file
        output_file (str, optional): Path to output file
        gene_dict (dict, optional): Dictionary of gene information

    Returns:
        dict: Dictionary of signal peptide predictions keyed by protein ID
    """
    # Check if input file exists
    if not os.path.isfile(input_file):
        logger.error(f"Input file not found: {input_file}")
        return None

    # Parse SignalP JSON output
    predictions = {}

    try:
        with open(input_file, "r") as f:
            data = json.load(f)

            for entry in data:
                protein_id = entry["ID"]
                prediction = entry["PREDICTION"]
                probability = float(entry["PROB"])

                # Store prediction
                predictions[protein_id] = {
                    "prediction": prediction,
                    "probability": probability,
                    "has_signal_peptide": prediction != "OTHER",
                }

                # Add cleavage site if available
                if "CS_POS" in entry and entry["CS_POS"] != "-":
                    predictions[protein_id]["cleavage_site"] = entry["CS_POS"]

        # Write to output file if specified
        if output_file:
            with open(output_file, "w") as out:
                # Write header
                out.write("#gene_id\tannotation_type\tannotation_value\n")

                for protein_id, pred in predictions.items():
                    if pred["has_signal_peptide"]:
                        out.write(
                            f"{protein_id}\tnote\tSignalP:{pred['prediction']} prob={pred['probability']:.3f}"
                        )
                        if "cleavage_site" in pred:
                            out.write(f" cleavage_site={pred['cleavage_site']}")
                        out.write("\n")

        return predictions

    except Exception as e:
        logger.error(f"Error parsing SignalP JSON output: {e}")
        return None


def signalp_to_json(input_file, output_file):
    """
    Convert SignalP predictions to JSON format.

    Args:
        input_file (str): Path to SignalP output file
        output_file (str): Path to output JSON file

    Returns:
        bool: True if successful, False otherwise
    """
    # Check if input is already JSON
    if input_file.endswith(".json"):
        predictions = parse_signalp_json(input_file)
    else:
        predictions = parse_signalp(input_file)

    if predictions:
        with open(output_file, "w") as f:
            json.dump(predictions, f, indent=2)
        return True
    return False


def signalp_subparser(subparsers):
    """
    Add argparse subparser for SignalP functionality.

    Args:
        subparsers: argparse subparsers object
    """
    parser = subparsers.add_parser(
        "signalp",
        description="Run SignalP and parse results",
        help="Run SignalP and parse results",
    )

    # Input/output options
    input_group = parser.add_argument_group("Input/Output")
    input_group.add_argument(
        "-i", "--input", required=True, help="Input protein FASTA file"
    )
    input_group.add_argument("-o", "--output", required=True, help="Output prefix")
    input_group.add_argument(
        "--signalp_path", default="signalp", help="Path to signalp executable"
    )

    # SignalP options
    signalp_group = parser.add_argument_group("SignalP options")
    signalp_group.add_argument(
        "--format", default="short", choices=["short", "long"], help="Prediction format"
    )
    signalp_group.add_argument(
        "--organism",
        default="euk",
        choices=["euk", "gram+", "gram-"],
        help="Organism group",
    )
    signalp_group.add_argument(
        "--min-prob",
        type=float,
        default=0.5,
        help="Minimum probability for signal peptide prediction",
    )
    signalp_group.add_argument(
        "--output-format",
        default="standard",
        choices=["standard", "none"],
        help="Output format",
    )
    signalp_group.add_argument(
        "--no-fasta", action="store_false", dest="fasta", help="Disable FASTA output"
    )
    signalp_group.add_argument(
        "--mature", action="store_true", help="Output mature sequences"
    )
    signalp_group.add_argument(
        "--gff3", action="store_true", help="Output in GFF3 format"
    )
    signalp_group.add_argument(
        "--json-output", action="store_true", help="Output in JSON format"
    )

    # Output format options
    format_group = parser.add_argument_group("Output format")
    format_group.add_argument("--json", help="Output predictions in JSON format")

    # Set the function to call
    parser.set_defaults(func=run_signalp_cli)


def run_signalp_cli(args):
    """
    Command-line interface for running SignalP.

    Args:
        args: argparse arguments
    """
    # Set up logging
    log_file = f"{args.output}.signalp.log"
    startLogging(log_file=log_file)

    # Run SignalP
    output_file = run_signalp(
        input_file=args.input,
        output_prefix=args.output,
        signalp_path=args.signalp_path,
        format=args.format,
        organism=args.organism,
        min_prob=args.min_prob,
        output_format=args.output_format,
        fasta=args.fasta,
        mature=args.mature,
        gff3=args.gff3,
        json_output=args.json_output,
    )

    if output_file:
        # Parse predictions
        funannotate_file = f"{args.output}.signalp.annotations.txt"

        # Determine parser based on output format
        if args.json_output and output_file.endswith(".json"):
            predictions = parse_signalp_json(output_file, output_file=funannotate_file)
        else:
            predictions = parse_signalp(output_file, output_file=funannotate_file)

        if predictions:
            # Count signal peptides
            signal_count = sum(
                1 for pred in predictions.values() if pred["has_signal_peptide"]
            )

            logger.info(f"Parsed {len(predictions)} predictions from {output_file}")
            logger.info(f"Found {signal_count} proteins with signal peptides")
            logger.info(f"Wrote annotations to {funannotate_file}")

            # Convert to JSON if requested
            if args.json:
                if signalp_to_json(output_file, args.json):
                    logger.info(f"Wrote JSON predictions to {args.json}")
                else:
                    logger.error(f"Failed to write JSON predictions to {args.json}")
        else:
            logger.error(f"Failed to parse predictions from {output_file}")
    else:
        logger.error("Failed to run SignalP")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run SignalP and parse results")
    subparsers = parser.add_subparsers(dest="command")
    signalp_subparser(subparsers)
    args = parser.parse_args()

    if args.command == "signalp":
        run_signalp_cli(args)
    else:
        parser.print_help()
