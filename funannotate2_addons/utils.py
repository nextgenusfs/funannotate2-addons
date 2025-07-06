#!/usr/bin/env python3

import os
import re
import logging
import json
import shlex
from pathlib import Path
from typing import List, Union

# Set up logger
logger = logging.getLogger(__name__)


def log_command(cmd: Union[List[str], str]) -> None:
    """
    Log a command to the terminal using the logger.

    Args:
        cmd: The command to log, either as a list of arguments or a string.
    """
    if isinstance(cmd, list):
        # Format the command for display
        formatted_cmd = " ".join(shlex.quote(str(arg)) for arg in cmd)
    else:
        formatted_cmd = str(cmd)

    # Use the logger to output the command
    logger.info(f"Executing command: {formatted_cmd}")


def find_funannotate_predict_results(input_dir):
    """
    Find funannotate2 predict results in the given directory.

    Args:
        input_dir (str): Path to funannotate2 output directory

    Returns:
        dict: Dictionary with paths to relevant files
    """
    results = {}

    # Check if input_dir exists
    if not os.path.isdir(input_dir):
        logger.error(f"Input directory not found: {input_dir}")
        return results

    # Check if predict_results directory exists
    predict_results = os.path.join(input_dir, "predict_results")
    if not os.path.isdir(predict_results):
        logger.error(f"predict_results directory not found in {input_dir}")
        return results

    # Find protein FASTA file
    protein_files = [
        f
        for f in os.listdir(predict_results)
        if f.endswith(".proteins.fa") or f.endswith(".proteins.fasta")
    ]
    if protein_files:
        results["proteins"] = os.path.join(predict_results, protein_files[0])

    # Find GenBank file
    gbk_files = [
        f
        for f in os.listdir(predict_results)
        if f.endswith(".gbk") or f.endswith(".gb")
    ]
    if gbk_files:
        results["genbank"] = os.path.join(predict_results, gbk_files[0])

    # Find genome FASTA file
    genome_files = [
        f
        for f in os.listdir(predict_results)
        if f.endswith(".scaffolds.fa") or f.endswith(".scaffolds.fasta")
    ]
    if genome_files:
        results["genome"] = os.path.join(predict_results, genome_files[0])

    # Find transcripts FASTA file
    transcript_files = [
        f
        for f in os.listdir(predict_results)
        if f.endswith(".transcripts.fa") or f.endswith(".transcripts.fasta")
    ]
    if transcript_files:
        results["transcripts"] = os.path.join(predict_results, transcript_files[0])

    # Check if we found the necessary files
    if "proteins" not in results:
        logger.warning("Protein FASTA file not found in predict_results")

    return results


def parse_gff3_for_transcript_ids(gff3_file):
    """
    Parse GFF3 file to extract transcript IDs and their corresponding protein IDs.

    Args:
        gff3_file (str): Path to GFF3 file

    Returns:
        dict: Dictionary mapping protein IDs to transcript IDs
    """
    protein_to_transcript = {}

    # Check if GFF3 file exists
    if not os.path.isfile(gff3_file):
        logger.error(f"GFF3 file not found: {gff3_file}")
        return protein_to_transcript

    # Parse GFF3 file
    try:
        with open(gff3_file, "r") as f:
            current_mrna_id = None

            for line in f:
                # Skip comment lines
                if line.startswith("#"):
                    continue

                # Skip empty lines
                if not line.strip():
                    continue

                # Parse GFF3 line
                fields = line.strip().split("\t")

                # Skip if not enough fields
                if len(fields) < 9:
                    continue

                feature_type = fields[2]
                attributes = fields[8]

                # Extract mRNA ID
                if feature_type == "mRNA":
                    id_match = re.search(r"ID=([^;]+)", attributes)
                    if id_match:
                        current_mrna_id = id_match.group(1)

                # Extract protein ID and map to mRNA ID
                elif feature_type == "CDS":
                    protein_id_match = re.search(r"protein_id=([^;]+)", attributes)
                    if protein_id_match and current_mrna_id:
                        protein_id = protein_id_match.group(1)
                        protein_to_transcript[protein_id] = current_mrna_id

    except Exception as e:
        logger.error(f"Error parsing GFF3 file: {e}")

    return protein_to_transcript


def get_default_output_dir(input_dir):
    """
    Get the default output directory for funannotate2-addons.

    Args:
        input_dir (str): Path to funannotate2 predict output directory

    Returns:
        str: Path to default output directory
    """
    # If input_dir is a predict_results directory, use its parent
    if os.path.basename(input_dir) == "predict_results":
        base_dir = os.path.dirname(input_dir)
    else:
        base_dir = input_dir

    # Create annotate_misc directory if it doesn't exist
    output_dir = os.path.join(base_dir, "annotate_misc")
    os.makedirs(output_dir, exist_ok=True)

    return output_dir


def write_standardized_output(annotations, output_file, protein_to_transcript=None):
    """
    Write annotations to a standardized 3-column TSV file.

    Args:
        annotations (dict): Dictionary of annotations
        output_file (str): Path to output file
        protein_to_transcript (dict, optional): Dictionary mapping protein IDs to transcript IDs

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        with open(output_file, "w") as out:
            # Write header
            out.write("#gene_id\tannotation_type\tannotation_value\n")

            # Write annotations
            for entity_id, annot in annotations.items():
                # Map protein ID to transcript ID if mapping is provided
                if protein_to_transcript and entity_id in protein_to_transcript:
                    entity_id = protein_to_transcript[entity_id]

                # Handle different annotation formats
                if isinstance(annot, dict):
                    # Complex annotation structure (e.g., from emapper, iprscan)
                    if "gene_name" in annot and annot["gene_name"]:
                        out.write(f"{entity_id}\tname\t{annot['gene_name']}\n")

                    if "description" in annot and annot["description"]:
                        out.write(f"{entity_id}\tnote\t{annot['description']}\n")

                    if "cog_category" in annot and annot["cog_category"]:
                        out.write(f"{entity_id}\tnote\tCOG:{annot['cog_category']}\n")

                    if "eggnog_ogs" in annot and annot["eggnog_ogs"]:
                        out.write(f"{entity_id}\tnote\teggNOG:{annot['eggnog_ogs']}\n")

                    if "ec_number" in annot and annot["ec_number"]:
                        for ec in annot["ec_number"].split(","):
                            if ec.strip():
                                out.write(f"{entity_id}\tEC_number\t{ec.strip()}\n")

                    if "interpro_domains" in annot:
                        for domain in annot["interpro_domains"]:
                            out.write(
                                f"{entity_id}\tnote\tInterPro:{domain['id']} {domain.get('name', '')}\n"
                            )

                    if "go_terms" in annot:
                        for go in annot["go_terms"]:
                            out.write(f"{entity_id}\tgo_term\t{go['id']}\n")

                    if "pathways" in annot:
                        for pathway in annot["pathways"]:
                            out.write(
                                f"{entity_id}\tnote\tPathway:{pathway['database']}:{pathway['id']} {pathway.get('name', '')}\n"
                            )

                    if "prediction" in annot and annot.get("has_signal_peptide", False):
                        out.write(
                            f"{entity_id}\tnote\tSignalP:{annot['prediction']} prob={annot['probability']:.3f}"
                        )
                        if "cleavage_site" in annot:
                            out.write(f" cleavage_site={annot['cleavage_site']}")
                        out.write("\n")

                elif isinstance(annot, list):
                    # List of annotations (e.g., from antismash)
                    for item in annot:
                        if isinstance(item, dict):
                            # Handle dictionary items
                            if "type" in item and "value" in item:
                                out.write(
                                    f"{entity_id}\t{item['type']}\t{item['value']}\n"
                                )
                        else:
                            # Handle simple items
                            out.write(f"{entity_id}\tnote\t{item}\n")

                elif isinstance(annot, str):
                    # Simple string annotation
                    out.write(f"{entity_id}\tnote\t{annot}\n")

        return True

    except Exception as e:
        logger.error(f"Error writing standardized output: {e}")
        return False


def parse_funannotate_predict_dir(input_dir):
    """
    Parse funannotate2 predict output directory and return paths to relevant files.

    Args:
        input_dir (str): Path to funannotate2 output directory or predict_results directory

    Returns:
        dict: Dictionary with paths to relevant files
    """
    # Check if input_dir is a predict_results directory
    if os.path.basename(input_dir) == "predict_results" and os.path.isdir(input_dir):
        predict_results = input_dir
        results = {}
    else:
        # Try to find predict_results directory
        results = find_funannotate_predict_results(input_dir)
        predict_results = os.path.join(input_dir, "predict_results")

    # If we didn't find any results, return empty dictionary
    if not results and os.path.isdir(predict_results):
        # Find protein FASTA file
        protein_files = [
            f
            for f in os.listdir(predict_results)
            if f.endswith(".proteins.fa") or f.endswith(".proteins.fasta")
        ]
        if protein_files:
            results["proteins"] = os.path.join(predict_results, protein_files[0])

        # Find GenBank file
        gbk_files = [
            f
            for f in os.listdir(predict_results)
            if f.endswith(".gbk") or f.endswith(".gb")
        ]
        if gbk_files:
            results["genbank"] = os.path.join(predict_results, gbk_files[0])

        # Find GFF3 file
        gff_files = [f for f in os.listdir(predict_results) if f.endswith(".gff3")]
        if gff_files:
            results["gff3"] = os.path.join(predict_results, gff_files[0])

        # Find genome FASTA file
        genome_files = [
            f
            for f in os.listdir(predict_results)
            if f.endswith(".scaffolds.fa") or f.endswith(".scaffolds.fasta")
        ]
        if genome_files:
            results["genome"] = os.path.join(predict_results, genome_files[0])

        # Find transcripts FASTA file
        transcript_files = [
            f
            for f in os.listdir(predict_results)
            if f.endswith(".transcripts.fa") or f.endswith(".transcripts.fasta")
        ]
        if transcript_files:
            results["transcripts"] = os.path.join(predict_results, transcript_files[0])

    return results


def get_input_file(args, file_type):
    """
    Get input file from either direct file input or funannotate2 predict directory.

    Args:
        args: argparse arguments
        file_type (str): Type of file to get ('proteins', 'genbank', 'gff3', 'genome', 'transcripts')

    Returns:
        str: Path to input file
    """
    # Check if direct file input is provided
    if hasattr(args, "file") and args.file:
        return args.file

    # Otherwise, parse funannotate2 predict directory
    if hasattr(args, "input") and args.input:
        predict_files = parse_funannotate_predict_dir(args.input)
        if file_type in predict_files:
            return predict_files[file_type]

    return None
