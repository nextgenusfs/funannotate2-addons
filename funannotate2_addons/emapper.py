#!/usr/bin/env python3

import os
import sys
import uuid
import subprocess
import argparse
import json
from packaging import version
from .log import startLogging
from .utils import (
    get_default_output_dir,
    get_input_file,
    log_command,
)

# Initialize logger at module level
logger = startLogging()


def get_version(emapper_path):
    """
    Get the version of eggnog-mapper.

    Args:
        emapper_path (str): Path to the emapper.py executable

    Returns:
        str: Version of eggnog-mapper
    """
    try:
        cmd = [emapper_path, "--version"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Parse version from output
        for line in proc.stdout.split("\n"):
            if line.startswith("emapper-"):
                return line.split("emapper-")[1].strip()
        return None
    except subprocess.SubprocessError as e:
        logger.error(f"Error getting eggnog-mapper version: {e}")
        return None


def run_emapper(
    input_file,
    output_prefix,
    emapper_path="emapper.py",
    cpu=1,
    database=None,
    data_dir=None,
    temp_dir=None,
    resume=False,
    override=False,
    mode="diamond",
    tax_scope=None,
    target_orthologs=None,
    sensmode="default",
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
):
    """
    Run eggnog-mapper on a protein FASTA file.

    Args:
        input_file (str): Path to input protein FASTA file
        output_prefix (str): Prefix for output files
        emapper_path (str): Path to emapper.py executable
        cpu (int): Number of CPU cores to use
        database (str): Specify the target database for sequence searches
        data_dir (str): Directory with eggnog-mapper databases
        temp_dir (str): Where temporary files are created
        resume (bool): Continue previous run
        override (bool): Override existing output files
        mode (str): Search mode (diamond, hmmer)
        tax_scope (str): Taxonomic scope for orthology assignment
        target_orthologs (str): Type of orthologs to use for functional transfer
        sensmode (str): Diamond sensitivity mode
        no_annot (bool): Skip functional annotation, reporting only orthologs
        no_search (bool): Skip sequence search, use existing hits file
        dmnd_db (str): Path to diamond database
        seed_orthologs (str): Path to seed orthologs file
        report_orthologs (bool): Include NOG alignments in output
        go_evidence (str): Filter GO terms by evidence
        pfam_realign (bool): Realign PFAM domains
        seed_ortholog_score (float): Min score for seed orthologs
        seed_ortholog_evalue (float): Max E-value for seed orthologs
        query_cover (float): Min query coverage
        subject_cover (float): Min subject coverage
        matrix (str): Scoring matrix
        gapopen (int): Gap open penalty
        gapextend (int): Gap extend penalty
        evalue (float): Max E-value
        pident (float): Min percent identity
        query_sgcov (float): Min query segment coverage
        subject_sgcov (float): Min subject segment coverage

    Returns:
        str: Path to the annotations output file
    """
    # Check if input file exists
    if not os.path.isfile(input_file):
        logger.error(f"Input file not found: {input_file}")
        return None

    # Build command
    cmd = [emapper_path, "-i", input_file, "--output", output_prefix, "--cpu", str(cpu)]

    # Add optional arguments
    if database:
        cmd.extend(["--database", database])
    if data_dir:
        cmd.extend(["--data_dir", data_dir])
    if temp_dir:
        cmd.extend(["--temp_dir", temp_dir])
    if resume:
        cmd.append("--resume")
    if override:
        cmd.append("--override")
    if mode:
        cmd.extend(["-m", mode])
    if tax_scope:
        cmd.extend(["--tax_scope", tax_scope])
    if target_orthologs:
        cmd.extend(["--target_orthologs", target_orthologs])
    if sensmode:
        cmd.extend(["--sensmode", sensmode])
    if no_annot:
        cmd.append("--no_annot")
    if no_search:
        cmd.append("--no_search")
    if dmnd_db:
        cmd.extend(["--dmnd_db", dmnd_db])
    if seed_orthologs:
        cmd.extend(["--seed_orthologs", seed_orthologs])
    if report_orthologs:
        cmd.append("--report_orthologs")
    if go_evidence:
        cmd.extend(["--go_evidence", go_evidence])
    if pfam_realign:
        cmd.append("--pfam_realign")
    if seed_ortholog_score is not None:
        cmd.extend(["--seed_ortholog_score", str(seed_ortholog_score)])
    if seed_ortholog_evalue is not None:
        cmd.extend(["--seed_ortholog_evalue", str(seed_ortholog_evalue)])
    if query_cover is not None:
        cmd.extend(["--query_cover", str(query_cover)])
    if subject_cover is not None:
        cmd.extend(["--subject_cover", str(subject_cover)])
    if matrix:
        cmd.extend(["--matrix", matrix])
    if gapopen is not None:
        cmd.extend(["--gapopen", str(gapopen)])
    if gapextend is not None:
        cmd.extend(["--gapextend", str(gapextend)])
    if evalue is not None:
        cmd.extend(["--evalue", str(evalue)])
    if pident is not None:
        cmd.extend(["--pident", str(pident)])
    if query_sgcov is not None:
        cmd.extend(["--query_sgcov", str(query_sgcov)])
    if subject_sgcov is not None:
        cmd.extend(["--subject_sgcov", str(subject_sgcov)])

    # Run command and log it to the terminal
    log_command(cmd)
    try:
        # Redirect stdout and stderr to the log file
        log_file = f"{output_prefix}.emapper.log"
        with open(log_file, "w") as log_fh:
            logger.info(f"Redirecting eggnog-mapper output to {log_file}")
            subprocess.run(cmd, check=True, stdout=log_fh, stderr=log_fh)

        # Return path to annotations file
        annotations_file = f"{output_prefix}.emapper.annotations"
        if os.path.isfile(annotations_file):
            logger.info(f"eggnog-mapper completed successfully: {annotations_file}")
            return annotations_file
        else:
            logger.error(f"eggnog-mapper output file not found: {annotations_file}")
            return None
    except subprocess.SubprocessError as e:
        logger.error(f"Error running eggnog-mapper: {e}")
        return None


def get_eggnog_headers(input_file):
    """
    Get the column indices for different data types in eggnog-mapper output.

    Args:
        input_file (str): Path to eggnog-mapper annotations file

    Returns:
        tuple: Indices for query_name, seed_ortholog, ortholog_group, gene_name,
               COG_category, description, EC_number
    """
    # Default indices (None means not found)
    query_idx, seed_ortholog_idx, og_idx, gene_idx = None, None, None, None
    cog_idx, desc_idx, ec_idx = None, None, None

    # Get eggnog-mapper version from file
    emapper_version = None
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith(("# emapper-", "## emapper-")):
                emapper_version = line.split("emapper-")[1].split()[0].strip()
                break

    # Parse header line
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith("#query") or line.startswith("# query"):
                # This is the header line
                # Handle both formats: "#query" and "# query"
                header = line.strip()
                if header.startswith("# "):
                    header = header[2:]
                elif header.startswith("#"):
                    header = header[1:]
                header = header.split("\t")

                # Log the header for debugging
                logger.debug(f"Found header: {header}")

                # Find indices for each column
                for i, col in enumerate(header):
                    if col == "query":
                        query_idx = i
                    elif col == "seed_ortholog":
                        seed_ortholog_idx = i
                    elif col == "eggNOG_OGs":
                        og_idx = i
                    elif col == "Preferred_name":
                        gene_idx = i
                    elif col == "COG_category":
                        cog_idx = i
                    elif col == "Description":
                        desc_idx = i
                    elif col == "EC":
                        ec_idx = i
                break

    # If no header found, try to guess based on version
    if query_idx is None and emapper_version:
        logger.warning(
            f"No header found in {input_file}, guessing based on version {emapper_version}"
        )

        # Try to read the header line again with more flexible matching
        with open(input_file, "r") as f:
            for line in f:
                if line.startswith("#") and "\t" in line and not line.startswith("##"):
                    # This might be the header line
                    header = line.strip()
                    if header.startswith("# "):
                        header = header[2:]
                    elif header.startswith("#"):
                        header = header[1:]
                    header = header.split("\t")

                    logger.debug(f"Trying alternative header: {header}")

                    # Try to find columns with more flexible matching
                    for i, col in enumerate(header):
                        if col.lower() == "query" or "query" in col.lower():
                            query_idx = i
                            logger.debug(f"Found query column at index {i}: {col}")
                        elif (
                            col.lower() == "seed_ortholog"
                            or "seed" in col.lower()
                            or "ortholog" in col.lower()
                        ):
                            seed_ortholog_idx = i
                            logger.debug(
                                f"Found seed_ortholog column at index {i}: {col}"
                            )
                        elif "eggnog" in col.lower() or "og" in col.lower():
                            og_idx = i
                            logger.debug(f"Found eggNOG_OGs column at index {i}: {col}")
                        elif "preferred" in col.lower() or "name" in col.lower():
                            gene_idx = i
                            logger.debug(
                                f"Found Preferred_name column at index {i}: {col}"
                            )
                        elif "cog" in col.lower() and "category" in col.lower():
                            cog_idx = i
                            logger.debug(
                                f"Found COG_category column at index {i}: {col}"
                            )
                        elif col.lower() == "description" or "desc" in col.lower():
                            desc_idx = i
                            logger.debug(
                                f"Found Description column at index {i}: {col}"
                            )
                        elif col.lower() == "ec" or "ec_" in col.lower():
                            ec_idx = i
                            logger.debug(f"Found EC column at index {i}: {col}")

                    if query_idx is not None:
                        break

        # If still no header found, use default indices based on version
        if query_idx is None:
            logger.warning(
                f"Still no header found, using default indices for version {emapper_version}"
            )

            # Web-based eggnog mapper has no header
            if version.parse(emapper_version) < version.parse("2.0.0"):
                # Old web-based format
                query_idx, seed_ortholog_idx = 0, 1
                og_idx, gene_idx = 4, 5
                cog_idx, desc_idx = 6, 7
                ec_idx = None
            elif version.parse(emapper_version) < version.parse("2.1.2"):
                # Version < 2.1.2
                query_idx, seed_ortholog_idx = 0, 1
                og_idx, gene_idx = 4, 5
                cog_idx, desc_idx = 6, 7
                ec_idx = 13
            else:
                # Version >= 2.1.2
                query_idx, seed_ortholog_idx = 0, 1
                og_idx, gene_idx = 4, 5
                cog_idx, desc_idx = 6, 7
                ec_idx = 13

    return query_idx, seed_ortholog_idx, og_idx, gene_idx, cog_idx, desc_idx, ec_idx


def parse_emapper_annotations(input_file, output_file=None, gene_dict=None):
    """
    Parse eggnog-mapper annotations file and convert to a standardized format.

    Args:
        input_file (str): Path to eggnog-mapper annotations file
        output_file (str, optional): Path to output file
        gene_dict (dict, optional): Dictionary of gene information

    Returns:
        dict: Dictionary of annotations keyed by gene ID
    """
    # Get column indices
    query_idx, seed_ortholog_idx, og_idx, gene_idx, cog_idx, desc_idx, ec_idx = (
        get_eggnog_headers(input_file)
    )

    # Log the indices for debugging
    logger.debug(
        f"Column indices: query={query_idx}, seed_ortholog={seed_ortholog_idx}, og={og_idx}, gene={gene_idx}, cog={cog_idx}, desc={desc_idx}, ec={ec_idx}"
    )

    if query_idx is None:
        logger.error(f"Could not determine column indices for {input_file}")
        # Try to read the first few lines of the file to help with debugging
        try:
            with open(input_file, "r") as f:
                first_lines = [next(f) for _ in range(10) if f]
                logger.debug(
                    f"First few lines of {input_file}:\n{''.join(first_lines)}"
                )
        except Exception as e:
            logger.debug(f"Error reading file for debugging: {e}")
        return None

    # Try to determine eggNOG version from file
    emapper_version = None
    prefix = "ENOG50"  # Default prefix if we can't determine it
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith(("# emapper-", "## emapper-")):
                emapper_version = line.split("emapper-")[1].split()[0].strip()
                break

    logger.debug(f"Detected eggNOG-mapper version: {emapper_version}")

    # Determine version for parsing logic
    if emapper_version:
        from packaging import version as pkg_version

        version_obj = pkg_version.parse(emapper_version)
        version_lt_2 = version_obj < pkg_version.parse("2.0.0")
        version_lt_212 = version_obj < pkg_version.parse("2.1.2")
    else:
        # Default to latest version behavior if we can't determine
        version_lt_2 = False
        version_lt_212 = False

    logger.debug(f"Version < 2.0.0: {version_lt_2}, Version < 2.1.2: {version_lt_212}")

    # Parse annotations
    annotations = {}
    definitions = {}  # Store NOG definitions
    line_count = 0
    parsed_count = 0
    error_count = 0

    with open(input_file, "r") as f:
        for line_num, line in enumerate(f, 1):
            line_count += 1

            # Skip comment lines
            if line.startswith("#"):
                continue

            # Skip empty lines
            if not line.strip():
                continue

            try:
                fields = line.strip().split("\t")
                # Replace "-" with empty string for easier handling
                fields = ["" if x == "-" else x for x in fields]

                # Log the line for debugging if it's one of the first few data lines
                if parsed_count < 3:
                    logger.debug(f"Line {line_num}: {fields}")

                # Skip if not enough fields
                max_idx = max(
                    filter(
                        None,
                        [
                            query_idx,
                            seed_ortholog_idx,
                            og_idx,
                            gene_idx,
                            cog_idx,
                            desc_idx,
                            ec_idx,
                        ],
                    )
                )

                if len(fields) <= max_idx:
                    logger.debug(
                        f"Line {line_num} has {len(fields)} fields, but need at least {max_idx + 1}"
                    )
                    error_count += 1
                    continue

                # Get values with additional error checking
                try:
                    query = fields[query_idx] if query_idx is not None else None
                    if not query:
                        logger.debug(f"Line {line_num}: Invalid query ID: {query}")
                        error_count += 1
                        continue

                    seed_ortholog = (
                        fields[seed_ortholog_idx]
                        if seed_ortholog_idx is not None
                        and seed_ortholog_idx < len(fields)
                        else None
                    )
                    og_field = (
                        fields[og_idx]
                        if og_idx is not None and og_idx < len(fields)
                        else None
                    )
                    gene_name = (
                        fields[gene_idx]
                        if gene_idx is not None and gene_idx < len(fields)
                        else None
                    )
                    cog_category = (
                        fields[cog_idx]
                        if cog_idx is not None and cog_idx < len(fields)
                        else None
                    )
                    description = (
                        fields[desc_idx]
                        if desc_idx is not None and desc_idx < len(fields)
                        else None
                    )
                    ec_number = (
                        fields[ec_idx]
                        if ec_idx is not None
                        and ec_idx < len(fields)
                        and ec_idx < len(fields)
                        else None
                    )

                    # Take only the first sentence of the description if it contains multiple sentences
                    if description and ". " in description:
                        description = description.split(". ")[0]

                    # Parse orthologous groups more carefully based on version
                    nog_id = None
                    db = None

                    if version_lt_2:  # version < 2.0.0
                        if og_field:
                            db = seed_ortholog.split("[")[0] if seed_ortholog else None
                            ogs = og_field.split(",")
                            for x in ogs:
                                if db and db in x:
                                    nog_id = f"{prefix}{x.split('@')[0]}"

                        # No EC number in older versions
                        ec_number = None

                    elif version_lt_212:  # version < 2.1.2
                        if og_field:
                            try:
                                nog, db = og_field.split("@")
                                nog_id = f"{prefix}{nog}"
                            except ValueError:
                                # If we can't split, try to get from seed_ortholog
                                if seed_ortholog and "@" in seed_ortholog:
                                    ogs = seed_ortholog.split(",")
                                    if len(ogs) > 1:
                                        nog, db = ogs[-2].split("@")
                                        nog_id = f"{prefix}{nog}"

                        # Format EC numbers
                        if ec_number and "," in ec_number:
                            # Use common prefix approach
                            ec_parts = ec_number.split(",")
                            import os

                            common_prefix = os.path.commonprefix(ec_parts).rstrip(".")
                            if common_prefix:
                                ec_number = common_prefix

                    else:  # version >= 2.1.2
                        if og_field:
                            db = og_field
                            ogs = seed_ortholog.split(",") if seed_ortholog else []
                            for ogx in ogs:
                                if "@" in ogx:
                                    nog_acc, taxname = ogx.split("@")
                                    if taxname == db:
                                        nog_id = f"{prefix}{nog_acc}"

                        # Format EC numbers
                        if ec_number and "," in ec_number:
                            # Use common prefix approach
                            ec_parts = ec_number.split(",")
                            import os

                            common_prefix = os.path.commonprefix(ec_parts).rstrip(".")
                            if common_prefix:
                                ec_number = common_prefix

                    # If we still don't have a NOG ID, try to extract it from eggNOG_OGs field
                    if not nog_id and og_field:
                        # Format is typically: OG@taxid|taxname,OG@taxid|taxname,...
                        ogs = og_field.split(",")
                        if ogs and "@" in ogs[0]:
                            nog = ogs[0].split("@", 1)[0]
                            nog_id = f"{prefix}{nog}"

                    # Format COG categories properly
                    formatted_cog = None
                    if cog_category:
                        cog_category = cog_category.replace(" ", "")
                        if len(cog_category) > 1:
                            formatted_cog = ",".join(list(cog_category))
                        else:
                            formatted_cog = cog_category

                    # Store annotations
                    if query not in annotations:
                        annotations[query] = {
                            "seed_ortholog": seed_ortholog,
                            "eggnog_ogs": og_field,  # Keep the original field
                            "nog_id": nog_id,  # Store the primary NOG ID
                            "gene_name": gene_name,
                            "cog_category": formatted_cog,
                            "description": description,
                            "ec_number": ec_number,
                        }

                        # Store NOG definition for later use
                        if nog_id and description and nog_id not in definitions:
                            definitions[nog_id] = description

                    parsed_count += 1

                    # Update gene_dict if provided
                    if gene_dict is not None and gene_name and description:
                        if (
                            not gene_name.startswith("_")
                            and "_" not in gene_name
                            and "." not in gene_name
                        ):
                            # Check if it's a valid gene name (contains at least one letter)
                            import re

                            if re.search("[a-zA-Z]", gene_name) and len(gene_name) > 2:
                                gene_id = query
                                if gene_id not in gene_dict:
                                    gene_dict[gene_id] = [
                                        {
                                            "name": gene_name,
                                            "product": description,
                                            "source": "EggNog-Mapper",
                                        }
                                    ]
                                else:
                                    gene_dict[gene_id].append(
                                        {
                                            "name": gene_name,
                                            "product": description,
                                            "source": "EggNog-Mapper",
                                        }
                                    )

                except Exception as e:
                    logger.debug(f"Error parsing line {line_num}: {e}")
                    error_count += 1
            except Exception as e:
                logger.debug(f"Error processing line {line_num}: {e}")
                error_count += 1

    logger.info(
        f"Processed {line_count} lines, parsed {parsed_count} annotations, encountered {error_count} errors"
    )

    # Write to output file if specified
    if output_file:
        with open(output_file, "w") as out:
            # Write header
            out.write("#gene_id\tannotation_type\tannotation_value\n")

            for query, annot in annotations.items():
                # Format for funannotate
                if annot["gene_name"]:
                    out.write(f"{query}\tname\t{annot['gene_name']}\n")

                if annot["description"]:
                    out.write(f"{query}\tnote\t{annot['description']}\n")

                if annot["cog_category"]:
                    out.write(f"{query}\tnote\tCOG:{annot['cog_category']}\n")

                # Write the primary NOG ID if available
                if annot["nog_id"]:
                    out.write(f"{query}\tnote\tEggNog:{annot['nog_id']}\n")
                # If no primary NOG ID but we have the original field, use that
                elif annot["eggnog_ogs"]:
                    out.write(f"{query}\tnote\teggNOG:{annot['eggnog_ogs']}\n")

                if annot["ec_number"]:
                    for ec in annot["ec_number"].split(","):
                        if ec.strip():
                            out.write(f"{query}\tEC_number\t{ec.strip()}\n")

    return annotations


def emapper_to_json(input_file, output_file):
    """
    Convert eggnog-mapper annotations to JSON format.

    Args:
        input_file (str): Path to eggnog-mapper annotations file
        output_file (str): Path to output JSON file

    Returns:
        bool: True if successful, False otherwise
    """
    annotations = parse_emapper_annotations(input_file)
    if annotations:
        # Create a simplified version for JSON output
        json_output = {}
        for query, annot in annotations.items():
            json_output[query] = {
                "seed_ortholog": annot.get("seed_ortholog"),
                "eggnog_ogs": annot.get("eggnog_ogs"),
                "nog_id": annot.get("nog_id"),
                "gene_name": annot.get("gene_name"),
                "cog_category": annot.get("cog_category"),
                "description": annot.get("description"),
                "ec_number": annot.get("ec_number"),
            }

        with open(output_file, "w") as f:
            json.dump(json_output, f, indent=2)
        return True
    return False


def emapper_subparser(subparsers):
    """
    Add argparse subparser for eggnog-mapper functionality.

    Args:
        subparsers: argparse subparsers object
    """
    parser = subparsers.add_parser(
        "emapper",
        description="Run eggnog-mapper and parse results",
        help="Run eggnog-mapper and parse results",
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
        help="Path to pre-computed eggNOG-mapper annotations file to parse (skips running eggnog-mapper)",
    )
    input_group.add_argument(
        "-o", "--output", help="Output directory (default: input_dir/annotate_misc)"
    )
    input_group.add_argument(
        "--emapper_path", default="emapper.py", help="Path to emapper.py executable"
    )

    # EggNOG-mapper options
    emapper_group = parser.add_argument_group("EggNOG-mapper options")
    emapper_group.add_argument(
        "--cpus", type=int, default=1, help="Number of CPU cores to use"
    )
    emapper_group.add_argument(
        "--database", help="Specify the target database for sequence searches"
    )
    emapper_group.add_argument(
        "--data_dir", help="Directory with eggnog-mapper databases"
    )
    emapper_group.add_argument("--temp_dir", help="Where temporary files are created")
    emapper_group.add_argument(
        "--resume", action="store_true", help="Continue previous run"
    )
    emapper_group.add_argument(
        "--override", action="store_true", help="Override existing output files"
    )
    emapper_group.add_argument(
        "--mode", default="diamond", choices=["diamond", "hmmer"], help="Search mode"
    )
    emapper_group.add_argument(
        "--tax_scope", default="auto", help="Taxonomic scope for orthology assignment"
    )
    emapper_group.add_argument(
        "--target_orthologs",
        default="all",
        choices=["one2one", "many2one", "one2many", "many2many", "all"],
        help="Type of orthologs to use for functional transfer",
    )
    emapper_group.add_argument(
        "--sensmode",
        default="sensitive",
        choices=[
            "default",
            "fast",
            "mid-sensitive",
            "sensitive",
            "more-sensitive",
            "very-sensitive",
            "ultra-sensitive",
        ],
        help="Diamond sensitivity mode",
    )
    emapper_group.add_argument(
        "--no_annot",
        action="store_true",
        help="Skip functional annotation, reporting only orthologs",
    )
    emapper_group.add_argument(
        "--no_search",
        action="store_true",
        help="Skip sequence search, use existing hits file",
    )
    emapper_group.add_argument("--dmnd_db", help="Path to diamond database")
    emapper_group.add_argument("--seed_orthologs", help="Path to seed orthologs file")
    emapper_group.add_argument(
        "--report_orthologs",
        action="store_true",
        help="Include NOG alignments in output",
    )
    emapper_group.add_argument("--go_evidence", help="Filter GO terms by evidence")
    emapper_group.add_argument(
        "--pfam_realign", action="store_true", help="Realign PFAM domains"
    )

    # Advanced options
    advanced_group = parser.add_argument_group("Advanced options")
    advanced_group.add_argument(
        "--seed_ortholog_score", type=float, help="Min score for seed orthologs"
    )
    advanced_group.add_argument(
        "--seed_ortholog_evalue", type=float, help="Max E-value for seed orthologs"
    )
    advanced_group.add_argument("--query_cover", type=float, help="Min query coverage")
    advanced_group.add_argument(
        "--subject_cover", type=float, help="Min subject coverage"
    )
    advanced_group.add_argument("--matrix", help="Scoring matrix")
    advanced_group.add_argument("--gapopen", type=int, help="Gap open penalty")
    advanced_group.add_argument("--gapextend", type=int, help="Gap extend penalty")
    advanced_group.add_argument("--evalue", type=float, help="Max E-value")
    advanced_group.add_argument("--pident", type=float, help="Min percent identity")
    advanced_group.add_argument(
        "--query_sgcov", type=float, help="Min query segment coverage"
    )
    advanced_group.add_argument(
        "--subject_sgcov", type=float, help="Min subject segment coverage"
    )

    # Output format options
    format_group = parser.add_argument_group("Output format")
    format_group.add_argument(
        "--json", action="store_true", help="Output annotations in JSON format"
    )

    # Set the function to call
    parser.set_defaults(func=run_emapper_cli)


def run_emapper_cli(args):
    """
    Command-line interface for running eggnog-mapper.

    Args:
        args: argparse arguments
    """
    # Check if at least one input option is provided
    if not args.input and not args.file and not args.parse:
        logger.error("Either --input, --file, or --parse must be specified")
        return

    # Handle pre-calculated annotations file
    if args.parse:
        if not os.path.isfile(args.parse):
            logger.error(f"Annotations file not found: {args.parse}")
            return

        # Get output directory
        if args.output:
            output_dir = args.output
        else:
            # Use the directory of the annotations file
            output_dir = os.path.dirname(args.parse) or os.getcwd()

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Set output prefix based on the annotations file name
        base_name = os.path.basename(args.parse)
        if base_name.endswith(".emapper.annotations"):
            prefix = base_name[: -len(".emapper.annotations")]
        else:
            prefix = os.path.splitext(base_name)[0]

        output_prefix = os.path.join(output_dir, prefix)

        # Set up logging
        log_file = f"{output_prefix}.log"
        # Use the custom formatter via startLogging
        log = startLogging(logfile=log_file)

        log.info(f"Parsing pre-calculated annotations file: {args.parse}")

        # Parse annotations
        funannotate_file = f"{output_prefix}.annotations.txt"
        annotations = parse_emapper_annotations(
            args.parse, output_file=funannotate_file
        )

        if annotations:
            log.info(f"Parsed {len(annotations)} annotations from {args.parse}")
            log.info(f"Wrote annotations to {funannotate_file}")

            # Convert to JSON if requested
            if args.json:
                json_file = os.path.join(output_dir, f"{prefix}.json")
                if emapper_to_json(args.parse, json_file):
                    log.info(f"Wrote JSON annotations to {json_file}")
                else:
                    log.error(f"Failed to write JSON annotations to {json_file}")
        else:
            log.error(f"Failed to parse annotations from {args.parse}")

        return

    # Regular eggnog-mapper run with protein FASTA input
    # Get input file
    input_file = get_input_file(args, "proteins") if args.input else args.file

    if not input_file:
        log.error("No protein FASTA file found")
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
    output_prefix = os.path.join(output_dir, f"emapper_{uuid.uuid4().hex[:6]}")

    # Set up logging
    log_file = f"{output_prefix}.log"
    # Use the custom formatter via startLogging
    log = startLogging(logfile=log_file)

    log.info(f"Setting up eggnog-mapper run with input file: {input_file}")

    # Run eggnog-mapper
    annotations_file = run_emapper(
        input_file=input_file,
        output_prefix=output_prefix,
        emapper_path=args.emapper_path,
        cpu=args.cpus,
        database=args.database,
        data_dir=args.data_dir,
        temp_dir=args.temp_dir,
        resume=args.resume,
        override=args.override,
        mode=args.mode,
        tax_scope=args.tax_scope,
        target_orthologs=args.target_orthologs,
        sensmode=args.sensmode,
        no_annot=args.no_annot,
        no_search=args.no_search,
        dmnd_db=args.dmnd_db,
        seed_orthologs=args.seed_orthologs,
        report_orthologs=args.report_orthologs,
        go_evidence=args.go_evidence,
        pfam_realign=args.pfam_realign,
        seed_ortholog_score=args.seed_ortholog_score,
        seed_ortholog_evalue=args.seed_ortholog_evalue,
        query_cover=args.query_cover,
        subject_cover=args.subject_cover,
        matrix=args.matrix,
        gapopen=args.gapopen,
        gapextend=args.gapextend,
        evalue=args.evalue,
        pident=args.pident,
        query_sgcov=args.query_sgcov,
        subject_sgcov=args.subject_sgcov,
    )

    if annotations_file:
        # Parse annotations
        funannotate_file = f"{output_prefix}.annotations.txt"
        annotations = parse_emapper_annotations(
            annotations_file, output_file=funannotate_file
        )

        if annotations:
            log.info(f"Parsed {len(annotations)} annotations from {annotations_file}")
            log.info(f"Wrote annotations to {funannotate_file}")

            # Convert to JSON if requested
            if args.json:
                json_file = os.path.join(output_dir, "emapper.json")
                if emapper_to_json(annotations_file, json_file):
                    log.info(f"Wrote JSON annotations to {json_file}")
                else:
                    log.error(f"Failed to write JSON annotations to {json_file}")
        else:
            log.error(f"Failed to parse annotations from {annotations_file}")
    else:
        log.error("Failed to run eggnog-mapper")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run eggnog-mapper and parse results")
    subparsers = parser.add_subparsers(dest="command")
    emapper_subparser(subparsers)
    args = parser.parse_args()

    if args.command == "emapper":
        run_emapper_cli(args)
    else:
        parser.print_help()
