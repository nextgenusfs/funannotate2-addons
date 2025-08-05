#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from .log import startLogging
from .utils import (
    parse_funannotate_predict_dir,
    get_default_output_dir,
    get_input_file,
    log_command,
)

# Initialize logger at module level
logger = startLogging()


def get_version(iprscan_path):
    """
    Get the version of InterProScan.

    Args:
        iprscan_path (str): Path to the interproscan.sh script

    Returns:
        str: Version of InterProScan
    """
    try:
        cmd = [iprscan_path, "--version"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Parse version from output
        for line in proc.stdout.split("\n"):
            if "InterProScan version" in line:
                return line.split("version")[1].strip()
        return None
    except subprocess.SubprocessError as e:
        logger.error(f"Error getting InterProScan version: {e}")
        return None


def run_iprscan(
    input_file,
    output_prefix,
    iprscan_path="interproscan.sh",
    cpu=1,
    applications=None,
    goterms=True,
    pathways=True,
    iprlookup=True,
    appl_goterms=False,
    tempdir=None,
    disable_precalc=False,
    disable_residue_annot=False,
    seqtype="p",
    mode="standalone",
    format="XML",
):
    """
    Run InterProScan on a protein FASTA file.

    Args:
        input_file (str): Path to input protein FASTA file
        output_prefix (str): Prefix for output files
        iprscan_path (str): Path to interproscan.sh script
        cpu (int): Number of CPU cores to use
        applications (list): List of applications to run
        goterms (bool): Include GO terms
        pathways (bool): Include pathway annotations
        iprlookup (bool): Use InterPro lookup service
        appl_goterms (bool): Use application-specific GO terms
        tempdir (str): Directory for temporary files
        disable_precalc (bool): Disable precalculated match lookup
        disable_residue_annot (bool): Disable residue-level annotations
        seqtype (str): Sequence type (p: protein, n: nucleotide)
        mode (str): Mode (standalone, binary, convert)
        format (str): Output format (XML, JSON, TSV, GFF3, HTML)

    Returns:
        str: Path to the output file
    """
    # Check if input file exists
    if not os.path.isfile(input_file):
        logger.error(f"Input file not found: {input_file}")
        return None

    # Build command
    cmd = [
        iprscan_path,
        "--input",
        input_file,
        "-o",
        os.path.basename(f"{output_prefix}.{format.lower()}"),
        "-f",
        format,
        "-cpu",
        str(cpu),
    ]

    # Add optional arguments
    if applications:
        cmd.extend(["-appl", ",".join(applications)])
    if goterms:
        cmd.append("-goterms")
    if pathways:
        cmd.append("-pa")
    if iprlookup:
        cmd.append("-iprlookup")
    if appl_goterms:
        cmd.append("-appl-goterms")
    if tempdir:
        cmd.extend(["-T", tempdir])
    if disable_precalc:
        cmd.append("-dp")
    if disable_residue_annot:
        cmd.append("-dra")
    if seqtype:
        cmd.extend(["-t", seqtype])
    if mode:
        cmd.extend(["-mode", mode])

    # Run command and log it to the terminal
    log_command(cmd)
    try:
        # Redirect stdout and stderr to the log file
        log_file = f"{output_prefix}.iprscan.log"
        with open(log_file, "w") as log_fh:
            logger.info(f"Redirecting InterProScan output to {log_file}")
            subprocess.run(cmd, check=True, stdout=log_fh, stderr=log_fh)

        # Return path to output file
        output_file = f"{output_prefix}.{format.lower()}"
        if os.path.isfile(output_file):
            logger.info(f"InterProScan completed successfully: {output_file}")
            return output_file
        else:
            logger.error(f"InterProScan output file not found: {output_file}")
            return None
    except subprocess.SubprocessError as e:
        logger.error(f"Error running InterProScan: {e}")
        return None


def parse_iprscan_xml(input_file, output_file=None, gene_dict=None):
    """
    Parse InterProScan XML output and convert to a standardized format.

    Args:
        input_file (str): Path to InterProScan XML output file
        output_file (str, optional): Path to output file
        gene_dict (dict, optional): Dictionary of gene information

    Returns:
        dict: Dictionary of annotations keyed by protein ID
    """
    # Check if input file exists
    if not os.path.isfile(input_file):
        logger.error(f"Input file not found: {input_file}")
        return None

    # Parse XML
    try:
        tree = ET.parse(input_file)
        root = tree.getroot()
    except ET.ParseError as e:
        logger.error(f"Error parsing XML file: {e}")
        return None

    # Define namespaces
    ns = {"ns": "http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5"}

    # Parse annotations
    annotations = {}

    # Iterate through protein matches
    for protein in root.findall(".//ns:protein", ns):
        protein_id = protein.attrib.get("id", "")

        if protein_id not in annotations:
            annotations[protein_id] = {
                "interpro_domains": [],
                "go_terms": [],
                "pathways": [],
                "signatures": [],
            }

        # Get matches
        for match in protein.findall(".//ns:match", ns):
            signature = match.find(".//ns:signature", ns)
            if signature is not None:
                sig_acc = signature.attrib.get("ac", "")
                sig_desc = signature.attrib.get("desc", "")
                sig_name = signature.attrib.get("name", "")

                # Add to signatures
                if sig_acc and sig_acc not in [
                    s["id"] for s in annotations[protein_id]["signatures"]
                ]:
                    annotations[protein_id]["signatures"].append(
                        {"id": sig_acc, "name": sig_name, "description": sig_desc}
                    )

                # Get InterPro domains
                entry = match.find(".//ns:entry", ns)
                if entry is not None:
                    entry_acc = entry.attrib.get("ac", "")
                    entry_desc = entry.attrib.get("desc", "")
                    entry_name = entry.attrib.get("name", "")

                    # Add to InterPro domains
                    if entry_acc and entry_acc not in [
                        d["id"] for d in annotations[protein_id]["interpro_domains"]
                    ]:
                        annotations[protein_id]["interpro_domains"].append(
                            {
                                "id": entry_acc,
                                "name": entry_name,
                                "description": entry_desc,
                            }
                        )

                # Get GO terms
                for go_term in match.findall(".//ns:go-term", ns):
                    go_id = go_term.attrib.get("id", "")
                    go_name = go_term.attrib.get("name", "")
                    go_category = go_term.attrib.get("category", "")

                    # Add to GO terms
                    if go_id and go_id not in [
                        g["id"] for g in annotations[protein_id]["go_terms"]
                    ]:
                        annotations[protein_id]["go_terms"].append(
                            {"id": go_id, "name": go_name, "category": go_category}
                        )

                # Get pathways
                for pathway in match.findall(".//ns:pathway", ns):
                    pathway_id = pathway.attrib.get("id", "")
                    pathway_name = pathway.attrib.get("name", "")
                    pathway_db = pathway.attrib.get("db", "")

                    # Add to pathways
                    if pathway_id and pathway_id not in [
                        p["id"] for p in annotations[protein_id]["pathways"]
                    ]:
                        annotations[protein_id]["pathways"].append(
                            {
                                "id": pathway_id,
                                "name": pathway_name,
                                "database": pathway_db,
                            }
                        )

    # Write to output file if specified
    if output_file:
        with open(output_file, "w") as out:
            # Write header
            out.write("#gene_id\tannotation_type\tannotation_value\n")

            for protein_id, annot in annotations.items():
                # Format for funannotate

                # Add InterPro domains
                for domain in annot["interpro_domains"]:
                    out.write(
                        f"{protein_id}\tnote\tInterPro:{domain['id']} {domain['name']}\n"
                    )

                # Add signatures
                for sig in annot["signatures"]:
                    out.write(f"{protein_id}\tnote\t{sig['id']} {sig['name']}\n")

                # Add GO terms
                for go in annot["go_terms"]:
                    out.write(f"{protein_id}\tgo_term\t{go['id']}\n")

                # Add pathways
                for pathway in annot["pathways"]:
                    out.write(
                        f"{protein_id}\tnote\tPathway:{pathway['database']}:{pathway['id']} {pathway['name']}\n"
                    )

    return annotations


def parse_iprscan_tsv(input_file, output_file=None, gene_dict=None):
    """
    Parse InterProScan TSV output and convert to a standardized format.

    Args:
        input_file (str): Path to InterProScan TSV output file
        output_file (str, optional): Path to output file
        gene_dict (dict, optional): Dictionary of gene information

    Returns:
        dict: Dictionary of annotations keyed by protein ID
    """
    # Check if input file exists
    if not os.path.isfile(input_file):
        logger.error(f"Input file not found: {input_file}")
        return None

    # Parse TSV
    annotations = {}

    with open(input_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")

            # Skip if not enough fields
            if len(fields) < 11:
                continue

            # Get values
            protein_id = fields[0]
            analysis = fields[3]
            signature_acc = fields[4]
            signature_desc = fields[5]
            start = fields[6]
            end = fields[7]
            evalue = fields[8] if fields[8] != "-" else None
            status = fields[9]
            date = fields[10]

            # InterPro annotations
            interpro_acc = (
                fields[11] if len(fields) > 11 and fields[11] != "-" else None
            )
            interpro_desc = (
                fields[12] if len(fields) > 12 and fields[12] != "-" else None
            )

            # GO terms
            go_terms = []
            if len(fields) > 13 and fields[13] != "-":
                go_terms = fields[13].split("|")

            # Pathways
            pathways = []
            if len(fields) > 14 and fields[14] != "-":
                pathways = fields[14].split("|")

            # Initialize protein entry if not exists
            if protein_id not in annotations:
                annotations[protein_id] = {
                    "interpro_domains": [],
                    "go_terms": [],
                    "pathways": [],
                    "signatures": [],
                }

            # Add signature
            if signature_acc and signature_acc not in [
                s["id"] for s in annotations[protein_id]["signatures"]
            ]:
                annotations[protein_id]["signatures"].append(
                    {
                        "id": signature_acc,
                        "name": analysis,
                        "description": signature_desc,
                    }
                )

            # Add InterPro domain
            if interpro_acc and interpro_acc not in [
                d["id"] for d in annotations[protein_id]["interpro_domains"]
            ]:
                annotations[protein_id]["interpro_domains"].append(
                    {"id": interpro_acc, "name": "", "description": interpro_desc}
                )

            # Add GO terms
            for go_term in go_terms:
                # GO terms should be kept as whole (e.g., GO:0003677)
                # Don't split GO terms on colon
                if go_term and go_term not in [
                    g["id"] for g in annotations[protein_id]["go_terms"]
                ]:
                    annotations[protein_id]["go_terms"].append(
                        {"id": go_term, "name": "", "category": ""}
                    )

            # Add pathways
            for pathway in pathways:
                pathway_parts = pathway.split(":", 1)
                if len(pathway_parts) > 1:
                    pathway_db, pathway_id_name = pathway_parts
                    pathway_id, pathway_name = (
                        pathway_id_name.split(" - ", 1)
                        if " - " in pathway_id_name
                        else (pathway_id_name, "")
                    )

                    if pathway_id and pathway_id not in [
                        p["id"] for p in annotations[protein_id]["pathways"]
                    ]:
                        annotations[protein_id]["pathways"].append(
                            {
                                "id": pathway_id,
                                "name": pathway_name,
                                "database": pathway_db,
                            }
                        )

    # Write to output file if specified
    if output_file:
        with open(output_file, "w") as out:
            # Write header
            out.write("#gene_id\tannotation_type\tannotation_value\n")

            for protein_id, annot in annotations.items():
                # Format for funannotate

                # Add InterPro domains
                for domain in annot["interpro_domains"]:
                    out.write(
                        f"{protein_id}\tnote\tInterPro:{domain['id']} {domain['description']}\n"
                    )

                # Add signatures
                for sig in annot["signatures"]:
                    out.write(f"{protein_id}\tnote\t{sig['id']} {sig['description']}\n")

                # Add GO terms
                for go in annot["go_terms"]:
                    out.write(f"{protein_id}\tgo_term\t{go['id']}\n")

                # Add pathways
                for pathway in annot["pathways"]:
                    out.write(
                        f"{protein_id}\tnote\tPathway:{pathway['database']}:{pathway['id']} {pathway['name']}\n"
                    )

    return annotations


def parse_iprscan(input_file, output_file=None, gene_dict=None):
    """
    Parse InterProScan output and convert to a standardized format.

    Args:
        input_file (str): Path to InterProScan output file
        output_file (str, optional): Path to output file
        gene_dict (dict, optional): Dictionary of gene information

    Returns:
        dict: Dictionary of annotations keyed by protein ID
    """
    # Determine file format based on extension
    if input_file.endswith(".xml"):
        return parse_iprscan_xml(input_file, output_file, gene_dict)
    elif input_file.endswith(".tsv"):
        return parse_iprscan_tsv(input_file, output_file, gene_dict)
    else:
        logger.error(f"Unsupported file format: {input_file}")
        return None


def iprscan_to_json(input_file, output_file):
    """
    Convert InterProScan annotations to JSON format.

    Args:
        input_file (str): Path to InterProScan output file
        output_file (str): Path to output JSON file

    Returns:
        bool: True if successful, False otherwise
    """
    annotations = parse_iprscan(input_file)
    if annotations:
        with open(output_file, "w") as f:
            json.dump(annotations, f, indent=2)
        return True
    return False


def iprscan_subparser(subparsers):
    """
    Add argparse subparser for InterProScan functionality.

    Args:
        subparsers: argparse subparsers object
    """
    parser = subparsers.add_parser(
        "iprscan",
        description="Run InterProScan and parse results",
        help="Run InterProScan and parse results",
    )

    # Input/output options
    input_group = parser.add_argument_group("Input/Output")
    input_group.add_argument(
        "-i",
        "--input",
        help="Path to funannotate2 predict output directory",
    )
    input_group.add_argument(
        "-f",
        "--file",
        help="Path to protein FASTA file (alternative to --input)",
    )
    input_group.add_argument(
        "--parse",
        help="Path to pre-computed InterProScan output file (XML or TSV) to parse (skips running InterProScan)",
    )
    input_group.add_argument(
        "-o", "--output", help="Output directory (default: input_dir/annotate_misc)"
    )
    input_group.add_argument(
        "--iprscan_path",
        default="interproscan.sh",
        help="Path to interproscan.sh script",
    )

    # InterProScan options
    iprscan_group = parser.add_argument_group("InterProScan options")
    iprscan_group.add_argument(
        "--cpus", type=int, default=1, help="Number of CPU cores to use"
    )
    iprscan_group.add_argument(
        "--applications", help="Comma-separated list of applications to run"
    )
    iprscan_group.add_argument(
        "--no-goterms",
        action="store_false",
        dest="goterms",
        help="Disable GO term lookup",
    )
    iprscan_group.add_argument(
        "--no-pathways",
        action="store_false",
        dest="pathways",
        help="Disable pathway lookup",
    )
    iprscan_group.add_argument(
        "--no-iprlookup",
        action="store_false",
        dest="iprlookup",
        help="Disable InterPro lookup",
    )
    iprscan_group.add_argument(
        "--appl-goterms", action="store_true", help="Use application-specific GO terms"
    )
    iprscan_group.add_argument("--tempdir", help="Directory for temporary files")
    iprscan_group.add_argument(
        "--disable-precalc",
        action="store_true",
        help="Disable precalculated match lookup",
    )
    iprscan_group.add_argument(
        "--disable-residue-annot",
        action="store_true",
        help="Disable residue-level annotations",
    )
    iprscan_group.add_argument(
        "--seqtype",
        default="p",
        choices=["p", "n"],
        help="Sequence type (p: protein, n: nucleotide)",
    )
    iprscan_group.add_argument(
        "--mode",
        default="standalone",
        choices=["standalone", "binary", "convert"],
        help="Mode",
    )
    iprscan_group.add_argument(
        "--format",
        default="XML",
        choices=["XML", "JSON", "TSV", "GFF3", "HTML"],
        help="Output format",
    )

    # Output format options
    format_group = parser.add_argument_group("Output format")
    format_group.add_argument(
        "--json", action="store_true", help="Output annotations in JSON format"
    )

    # Set the function to call
    parser.set_defaults(func=run_iprscan_cli)


def run_iprscan_cli(args):
    """
    Command-line interface for running InterProScan.

    Args:
        args: argparse arguments
    """
    # Handle parse-only mode
    if args.parse:
        if not os.path.isfile(args.parse):
            logger.error(f"Parse file not found: {args.parse}")
            return

        # Get output directory
        if args.output:
            output_dir = args.output
        else:
            # Use the directory of the parse file
            output_dir = os.path.dirname(args.parse) or os.getcwd()

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        logger.info(f"Parsing pre-computed InterProScan results: {args.parse}")

        # Parse InterProScan results
        annotations_file = os.path.join(output_dir, "iprscan.annotations.txt")
        annotations = parse_iprscan(args.parse, output_file=annotations_file)

        if annotations:
            logger.info(f"Parsed {len(annotations)} annotations from {args.parse}")
            logger.info(f"Wrote annotations to {annotations_file}")

            # Convert to JSON if requested
            if args.json:
                json_file = os.path.join(output_dir, "iprscan.json")
                if iprscan_to_json(args.parse, json_file):
                    logger.info(f"Wrote JSON annotations to {json_file}")
                else:
                    logger.error(f"Failed to write JSON annotations to {json_file}")
        else:
            logger.error(f"Failed to parse annotations from {args.parse}")

        return

    # Check if at least one input option is provided
    if not args.input and not args.file:
        logger.error("Either --input, --file, or --parse must be specified")
        return

    # Get input file
    input_file = get_input_file(args, "proteins") if args.input else args.file

    if not input_file:
        logger.error("No protein FASTA file found")
        return

    # Get output directory
    if args.output:
        output_dir = args.output
    elif args.input:  # Only use default output dir if using funannotate2 predict dir
        output_dir = get_default_output_dir(args.input)
    else:
        # If using direct file input without output dir, use current directory
        output_dir = (
            os.path.dirname(args.file) if os.path.dirname(args.file) else os.getcwd()
        )

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Set output prefix
    output_prefix = os.path.join(output_dir, "iprscan")

    # We'll use the module-level logger which is already initialized

    # Parse applications if provided
    applications = args.applications.split(",") if args.applications else None

    # Run InterProScan
    output_file = run_iprscan(
        input_file=input_file,
        output_prefix=output_prefix,
        iprscan_path=args.iprscan_path,
        cpu=args.cpus,
        applications=applications,
        goterms=args.goterms,
        pathways=args.pathways,
        iprlookup=args.iprlookup,
        appl_goterms=args.appl_goterms,
        tempdir=args.tempdir,
        disable_precalc=args.disable_precalc,
        disable_residue_annot=args.disable_residue_annot,
        seqtype=args.seqtype,
        mode=args.mode,
        format=args.format,
    )

    if output_file:
        # Parse annotations
        funannotate_file = f"{output_prefix}.annotations.txt"
        annotations = parse_iprscan(output_file, output_file=funannotate_file)

        if annotations:
            logger.info(f"Parsed {len(annotations)} annotations from {output_file}")
            logger.info(f"Wrote annotations to {funannotate_file}")

            # Convert to JSON if requested
            if args.json:
                json_file = os.path.join(output_dir, "iprscan.json")
                if iprscan_to_json(output_file, json_file):
                    logger.info(f"Wrote JSON annotations to {json_file}")
                else:
                    logger.error(f"Failed to write JSON annotations to {json_file}")
        else:
            logger.error(f"Failed to parse annotations from {output_file}")
    else:
        logger.error("Failed to run InterProScan")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run InterProScan and parse results")
    subparsers = parser.add_subparsers(dest="command")
    iprscan_subparser(subparsers)
    args = parser.parse_args()

    if args.command == "iprscan":
        run_iprscan_cli(args)
    else:
        parser.print_help()
