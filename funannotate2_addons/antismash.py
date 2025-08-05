#!/usr/bin/env python3

import os
import subprocess
import argparse
import json
import uuid
import gb_io
from .log import startLogging
from .utils import (
    parse_funannotate_predict_dir,
    get_default_output_dir,
    get_input_file,
    log_command,
)

# Initialize logger at module level
logger = startLogging()


def get_version(antismash_path):
    """
    Get the version of antiSMASH.

    Args:
        antismash_path (str): Path to the antismash executable

    Returns:
        str: Version of antiSMASH
    """
    try:
        cmd = [antismash_path, "--version"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Parse version from output
        for line in proc.stdout.split("\n"):
            if "antiSMASH" in line and "version" in line:
                return line.split("version")[1].strip()
        return None
    except subprocess.SubprocessError as e:
        logger.error(f"Error getting antiSMASH version: {e}")
        return None


def run_antismash(
    input_file,
    output_dir,
    antismash_path="antismash",
    cpu=1,
    taxon="fungi",
    fullhmmer=False,
    cassis=False,
    clusterhmmer=False,
    tigrfam=False,
    asf=False,
    cc_mibig=False,
    cb_general=False,
    cb_subclusters=False,
    cb_knownclusters=False,
    pfam2go=False,
    rre=False,
    smcog_trees=False,
    tfbs=False,
    html_description="antiSMASH results",
):
    """
    Run antiSMASH on a GenBank file.

    Args:
        input_file (str): Path to input GenBank file (.gbk or .gbff)
        output_dir (str): Directory for output files
        antismash_path (str): Path to antismash executable
        cpu (int): Number of CPU cores to use
        taxon (str): Taxonomy of the input genome (bacteria, fungi)
        fullhmmer (bool): Run a whole-genome HMMer analysis using Pfam profiles
        cassis (bool): Motif based prediction of SM gene cluster regions
        clusterhmmer (bool): Run a cluster-limited HMMer analysis using Pfam profiles
        tigrfam (bool): Annotate clusters using TIGRFam profiles
        asf (bool): Run active site finder analysis
        cc_mibig (bool): Run a comparison against the MIBiG dataset
        cb_general (bool): Compare identified clusters against a database of antiSMASH-predicted clusters
        cb_subclusters (bool): Compare identified clusters against known subclusters
        cb_knownclusters (bool): Compare identified clusters against known gene clusters from the MIBiG database
        pfam2go (bool): Run Pfam to Gene Ontology mapping module
        rre (bool): Run RREFinder precision mode on all RiPP gene clusters
        smcog_trees (bool): Generate phylogenetic trees of sec. met. cluster orthologous groups
        tfbs (bool): Run TFBS finder on all gene clusters
        html_description (str): Description for HTML output (HTML is always enabled)

    Returns:
        str: Path to the output GenBank file
    """
    # Check if input file exists
    if not os.path.isfile(input_file):
        logger.error(f"Input file not found: {input_file}")
        return None

    # Create a unique output directory for antiSMASH (it requires an empty directory)
    # antiSMASH will fail if the output directory already exists and contains files
    unique_id = str(uuid.uuid4())[:8]  # Use first 8 characters of UUID
    antismash_output_dir = os.path.join(output_dir, f"antismash_{unique_id}")

    # Create the parent output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Create the unique antiSMASH output directory (should be empty)
    os.makedirs(antismash_output_dir, exist_ok=False)
    logger.info(f"Created unique antiSMASH output directory: {antismash_output_dir}")

    # Build command
    cmd = [
        antismash_path,
        input_file,
        "--output-dir",
        antismash_output_dir,
        "--cpus",
        str(cpu),
    ]

    # Add optional arguments
    if taxon:
        cmd.extend(["--taxon", taxon])
    # Always set genefinding to none since we're using GenBank files with existing gene models
    cmd.extend(["--genefinding-tool", "none"])

    # Additional analysis options (all disabled by default to match antiSMASH v8.0.0)
    if fullhmmer:
        cmd.append("--fullhmmer")
    if cassis:
        cmd.append("--cassis")
    if clusterhmmer:
        cmd.append("--clusterhmmer")
    if tigrfam:
        cmd.append("--tigrfam")
    if asf:
        cmd.append("--asf")
    if cc_mibig:
        cmd.append("--cc-mibig")
    if cb_general:
        cmd.append("--cb-general")
    if cb_subclusters:
        cmd.append("--cb-subclusters")
    if cb_knownclusters:
        cmd.append("--cb-knownclusters")
    if pfam2go:
        cmd.append("--pfam2go")
    if rre:
        cmd.append("--rre")
    if smcog_trees:
        cmd.append("--smcog-trees")
    if tfbs:
        cmd.append("--tfbs")

    # HTML output (enabled by default in antiSMASH v8.0.0, no --html flag needed)
    if html_description:
        cmd.extend(["--html-description", html_description])

    # Run command and log it to the terminal
    log_command(cmd)
    try:
        # Redirect stdout and stderr to the log file (in parent output dir to keep antiSMASH dir empty)
        log_file = os.path.join(output_dir, "antismash.log")
        with open(log_file, "w") as log_fh:
            logger.info(f"Redirecting antiSMASH output to {log_file}")
            subprocess.run(cmd, check=True, stdout=log_fh, stderr=log_fh)

        # Find the main output GenBank file (not region files)
        gbk_files = []
        for root, _, files in os.walk(antismash_output_dir):
            for file in files:
                if (
                    file.endswith(".gbk") or file.endswith(".gb")
                ) and ".region" not in file:
                    gbk_files.append(os.path.join(root, file))

        if gbk_files:
            # Prefer the main genome file over region files
            main_file = None
            for gbk_file in gbk_files:
                # Look for the main genome file (usually the input filename or organism name)
                if ".region" not in gbk_file and "cluster" not in gbk_file.lower():
                    main_file = gbk_file
                    break

            result_file = main_file if main_file else gbk_files[0]
            logger.info(f"antiSMASH completed successfully: {result_file}")
            return result_file
        else:
            logger.error(
                f"antiSMASH output GenBank file not found in {antismash_output_dir}"
            )
            return None
    except subprocess.SubprocessError as e:
        logger.error(f"Error running antiSMASH: {e}")
        return None


def _qualifiers_to_dict(qualifiers):
    """
    Convert gb_io qualifiers list to dictionary format.

    Args:
        qualifiers: List of gb_io.Qualifier objects

    Returns:
        dict: Dictionary with qualifier keys and lists of values
    """
    qual_dict = {}
    for qual in qualifiers:
        if qual.key not in qual_dict:
            qual_dict[qual.key] = []
        qual_dict[qual.key].append(qual.value)
    return qual_dict


def parse_antismash_gbk(
    input_file, output_dir=None, annotations_file=None, clusters_file=None
):
    """
    Parse antiSMASH GenBank output and extract cluster information.

    Args:
        input_file (str): Path to antiSMASH GenBank output file
        output_dir (str, optional): Directory to store extracted cluster information
        annotations_file (str, optional): Path to output annotations file
        clusters_file (str, optional): Path to output clusters file

    Returns:
        tuple: (backbone_domains, backbone_subtype, backbone_enzymes, cluster_genes)
    """
    # Check if input file exists
    if not os.path.isfile(input_file):
        logger.error(f"Input file not found: {input_file}")
        return None, None, None, None

    # Create output directory if specified and it doesn't exist
    if output_dir and not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Parse GenBank file
    backbone_domains = {}  # protID -> domains
    backbone_subtype = {}  # clusterID -> subtype
    backbone_enzymes = {}  # clusterID -> backbone enzymes
    cluster_genes = {}  # clusterID -> list of genes

    try:
        # Parse GenBank file using gb_io instead of BioPython
        with open(input_file, "rb") as f:
            records = list(gb_io.load(f))
        counter = 1
        for record in records:
            # Get regions as the primary cluster units (matches HTML interface)
            # antiSMASH v8.0.0 organizes clusters into regions
            clusters = [f for f in record.features if f.kind == "region"]

            for cluster in clusters:
                # Get cluster information
                # Convert qualifiers to dictionary format
                cluster_quals = _qualifiers_to_dict(cluster.qualifiers)

                # Extract region ID (now focusing on regions as primary clusters)
                if "region_number" not in cluster_quals:
                    continue

                # Create globally unique cluster ID across all records
                cluster_id = f"cluster_{counter}"

                # Get cluster type
                if "product" in cluster_quals:
                    cluster_type = cluster_quals["product"][0]
                    if cluster_id in backbone_subtype:
                        backbone_subtype[cluster_id].append(cluster_type)
                    else:
                        backbone_subtype[cluster_id] = [cluster_type]

                # Get cluster location
                # cluster_location = f"{record.id}:{cluster.location.start}:{cluster.location.end}"

                # Get genes in cluster
                cluster_genes[cluster_id] = []

                # Find CDS features within cluster boundaries
                for feature in record.features:
                    if feature.kind == "CDS":
                        if (
                            feature.location.start >= cluster.location.start
                            and feature.location.end <= cluster.location.end
                        ):
                            # Get gene ID
                            # Convert feature qualifiers to dictionary format
                            feature_quals = _qualifiers_to_dict(feature.qualifiers)

                            gene_id = None
                            if "locus_tag" in feature_quals:
                                gene_id = feature_quals["locus_tag"][0]
                            elif "protein_id" in feature_quals:
                                gene_id = feature_quals["protein_id"][0]
                            elif "gene" in feature_quals:
                                gene_id = feature_quals["gene"][0]

                            if gene_id:
                                cluster_genes[cluster_id].append(gene_id)

                                # Extract comprehensive annotation information
                                annotations = []
                                domains_found = {}  # Track domains to consolidate info

                                # Extract SMCOGs
                                if "gene_functions" in feature_quals:
                                    for gene_func in feature_quals["gene_functions"]:
                                        # Clean up multi-line qualifiers by removing newlines and extra whitespace
                                        gene_func_clean = " ".join(gene_func.split())

                                        if (
                                            "smcogs" in gene_func_clean.lower()
                                            and "SMCOG" in gene_func_clean
                                        ):
                                            # Extract SMCOG ID and description
                                            smcog_match = gene_func_clean.split(
                                                "SMCOG"
                                            )[1]
                                            if ":" in smcog_match:
                                                smcog_id = (
                                                    "SMCOG" + smcog_match.split(":")[0]
                                                )
                                                smcog_desc = smcog_match.split(":", 1)[
                                                    1
                                                ].strip()
                                                annotations.append(
                                                    f"note\t{smcog_id}: {smcog_desc}"
                                                )

                                        # Extract biosynthetic domain names for consolidation
                                        elif "biosynthetic" in gene_func_clean.lower():
                                            if "rule-based-clusters" in gene_func_clean:
                                                # Extract domain name
                                                parts = gene_func_clean.split(")")
                                                if len(parts) >= 2:
                                                    domain_name = (
                                                        parts[1].split(":")[-1].strip()
                                                    )
                                                    domains_found[domain_name] = {
                                                        "source": "biosynthetic"
                                                    }

                                # Get detailed domain info from sec_met_domain and consolidate
                                if "sec_met_domain" in feature_quals:
                                    for domain in feature_quals["sec_met_domain"]:
                                        domain_clean = " ".join(domain.split())
                                        if "(" in domain_clean:
                                            domain_name = domain_clean.split("(")[
                                                0
                                            ].strip()
                                            domain_info = domain_clean.split("(")[
                                                1
                                            ].split(")")[0]
                                            # Extract E-value for summary
                                            evalue = "unknown"
                                            if "E-value:" in domain_info:
                                                evalue = (
                                                    domain_info.split("E-value:")[1]
                                                    .split(",")[0]
                                                    .strip()
                                                )

                                            # Store or update domain info
                                            domains_found[domain_name] = {
                                                "source": "sec_met",
                                                "evalue": evalue,
                                            }

                                # Get NRPS/PKS domain information
                                if "NRPS_PKS" in feature_quals:
                                    for nrps_pks in feature_quals["NRPS_PKS"]:
                                        nrps_pks_clean = " ".join(nrps_pks.split())
                                        if "Domain:" in nrps_pks_clean:
                                            # Extract domain name: "Domain: PKS_KR (26-118). E-value: 7.5e-08..."
                                            domain_part = nrps_pks_clean.split(
                                                "Domain:"
                                            )[1].strip()
                                            domain_name = domain_part.split("(")[
                                                0
                                            ].strip()
                                            # Extract E-value
                                            evalue = "unknown"
                                            if "E-value:" in nrps_pks_clean:
                                                evalue = (
                                                    nrps_pks_clean.split("E-value:")[1]
                                                    .split(".")[0]
                                                    .strip()
                                                )

                                            domains_found[domain_name] = {
                                                "source": "NRPS_PKS",
                                                "evalue": evalue,
                                            }

                                # Create consolidated domain annotations
                                for domain_name, info in domains_found.items():
                                    if "evalue" in info and info["evalue"] != "unknown":
                                        annotations.append(
                                            f"note\tantiSMASH domain: {domain_name} (E-value: {info['evalue']})"
                                        )
                                    else:
                                        annotations.append(
                                            f"note\tantiSMASH domain: {domain_name}"
                                        )

                                # Store all annotations for this gene
                                if annotations:
                                    if gene_id not in backbone_domains:
                                        backbone_domains[gene_id] = []
                                    backbone_domains[gene_id].extend(annotations)

                                # Check if this is a backbone enzyme (biosynthetic gene)
                                if "gene_functions" in feature_quals:
                                    for gene_func in feature_quals["gene_functions"]:
                                        # Clean up multi-line qualifiers
                                        gene_func_clean = " ".join(gene_func.split())
                                        if (
                                            "biosynthetic" in gene_func_clean.lower()
                                            and "rule-based-clusters" in gene_func_clean
                                        ):
                                            if cluster_id not in backbone_enzymes:
                                                backbone_enzymes[cluster_id] = []
                                            if (
                                                gene_id
                                                not in backbone_enzymes[cluster_id]
                                            ):
                                                backbone_enzymes[cluster_id].append(
                                                    gene_id
                                                )
                counter += 1

        # Write annotations file if specified
        if annotations_file:
            with open(annotations_file, "w") as f:
                # Write header
                f.write("#gene_id\tannotation_type\tannotation_value\n")

                # Write cluster membership annotations
                for cluster_id, genes in cluster_genes.items():
                    for gene_id in genes:
                        cluster_type = backbone_subtype.get(cluster_id, ["unknown"])
                        f.write(
                            f"{gene_id}\tnote\tantiSMASH: {cluster_id} [{','.join(cluster_type)}]\n"
                        )

                # Write detailed annotations (SMCOGs, domains, etc.)
                for gene_id, annotations in backbone_domains.items():
                    if isinstance(annotations, list):
                        for annotation in annotations:
                            if "\t" in annotation:
                                ann_type, ann_value = annotation.split("\t", 1)
                                f.write(f"{gene_id}\t{ann_type}\t{ann_value}\n")
                    else:
                        # Handle old format (string) - shouldn't happen with new code
                        f.write(f"{gene_id}\tsec_met_domain\t{annotations}\n")

        # Write clusters file if specified
        if clusters_file:
            with open(clusters_file, "w") as f:
                for cluster_id, genes in cluster_genes.items():
                    cluster_type = backbone_subtype.get(cluster_id, "unknown")
                    backbone = ", ".join(backbone_enzymes.get(cluster_id, []))
                    f.write(f"#Cluster: {cluster_id}\n")
                    f.write(f"#Type: {cluster_type}\n")
                    f.write(f"#Backbone enzymes: {backbone}\n")
                    f.write(f"#Genes: {', '.join(genes)}\n\n")

        return backbone_domains, backbone_subtype, backbone_enzymes, cluster_genes

    except Exception as e:
        logger.error(f"Error parsing antiSMASH GenBank file: {e}")
        return None, None, None, None


def extract_clusters(input_file, output_dir):
    """
    Extract individual cluster GenBank files from antiSMASH output.

    Args:
        input_file (str): Path to antiSMASH GenBank output file
        output_dir (str): Directory to store extracted cluster files

    Returns:
        list: Paths to extracted cluster files
    """
    # Check if input file exists
    if not os.path.isfile(input_file):
        logger.error(f"Input file not found: {input_file}")
        return []

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Parse GenBank file
    cluster_files = []

    try:
        # Parse GenBank file using gb_io instead of BioPython
        with open(input_file, "rb") as f:
            records = list(gb_io.load(f))

        for record_idx, record in enumerate(records):
            # Get regions as the primary cluster units (matches HTML interface)
            # antiSMASH v8.0.0 organizes clusters into regions
            clusters = [f for f in record.features if f.kind == "region"]

            for cluster in clusters:
                # Get cluster information
                # Convert qualifiers to dictionary format
                cluster_quals = _qualifiers_to_dict(cluster.qualifiers)

                # Extract region ID (now focusing on regions as primary clusters)
                if "region_number" not in cluster_quals:
                    continue

                # Create globally unique cluster ID across all records
                region_num = cluster_quals["region_number"][0]
                cluster_id = f"chr{record_idx + 1}_region_{region_num}"

                # Get cluster type
                if "product" in cluster_quals:
                    cluster_type = cluster_quals["product"][0]
                    cluster_id = f"{cluster_id}_{cluster_type.replace(' ', '_')}"

                # Create a new record for the cluster
                cluster_record = record[cluster.location.start : cluster.location.end]

                # Write cluster to file
                cluster_file = os.path.join(output_dir, f"{cluster_id}.gbk")
                with open(cluster_file, "wb") as f:
                    gb_io.write(cluster_record, f)
                cluster_files.append(cluster_file)

        return cluster_files

    except Exception as e:
        logger.error(f"Error extracting clusters from antiSMASH GenBank file: {e}")
        return []


def antismash_to_json(input_file, output_file):
    """
    Convert antiSMASH cluster information to JSON format.

    Args:
        input_file (str): Path to antiSMASH GenBank output file
        output_file (str): Path to output JSON file

    Returns:
        bool: True if successful, False otherwise
    """
    backbone_domains, backbone_subtype, backbone_enzymes, cluster_genes = (
        parse_antismash_gbk(input_file)
    )

    if backbone_domains is None:
        return False

    # Create JSON structure
    data = {
        "clusters": {},
        "backbone_domains": backbone_domains,
        "backbone_enzymes": backbone_enzymes,
    }

    for cluster_id, genes in cluster_genes.items():
        data["clusters"][cluster_id] = {
            "type": backbone_subtype.get(cluster_id, "unknown"),
            "genes": genes,
            "backbone_enzymes": backbone_enzymes.get(cluster_id, []),
        }

    # Write to JSON file
    try:
        with open(output_file, "w") as f:
            json.dump(data, f, indent=2)
        return True
    except Exception as e:
        logger.error(f"Error writing antiSMASH JSON file: {e}")
        return False


def antismash_subparser(subparsers):
    """
    Add argparse subparser for antiSMASH functionality.

    Args:
        subparsers: argparse subparsers object
    """
    parser = subparsers.add_parser(
        "antismash",
        description="Run antiSMASH and parse results",
        help="Run antiSMASH and parse results",
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
        help="Path to GenBank file (.gbk or .gbff) (alternative to --input)",
    )
    input_group.add_argument(
        "--parse",
        help="Path to pre-computed antiSMASH GenBank output file to parse (skips running antiSMASH)",
    )
    input_group.add_argument(
        "-o", "--output", help="Output directory (default: input_dir/annotate_misc)"
    )
    input_group.add_argument(
        "--antismash_path", default="antismash", help="Path to antismash executable"
    )

    # antiSMASH options
    antismash_group = parser.add_argument_group("antiSMASH options")
    antismash_group.add_argument(
        "--cpus", type=int, default=1, help="Number of CPU cores to use"
    )
    antismash_group.add_argument(
        "--taxon",
        default="fungi",
        choices=["bacteria", "fungi"],
        help="Taxonomy of the input genome (default: fungi)",
    )

    # Additional analysis options (disabled by default to match antiSMASH v8.0.0)
    additional_group = parser.add_argument_group("Additional analysis")
    additional_group.add_argument(
        "--fullhmmer",
        action="store_true",
        help="Run a whole-genome HMMer analysis using Pfam profiles",
    )
    additional_group.add_argument(
        "--cassis",
        action="store_true",
        help="Motif based prediction of SM gene cluster regions",
    )
    additional_group.add_argument(
        "--clusterhmmer",
        action="store_true",
        help="Run a cluster-limited HMMer analysis using Pfam profiles",
    )
    additional_group.add_argument(
        "--tigrfam",
        action="store_true",
        help="Annotate clusters using TIGRFam profiles",
    )
    additional_group.add_argument(
        "--asf",
        action="store_true",
        help="Run active site finder analysis",
    )
    additional_group.add_argument(
        "--cc-mibig",
        action="store_true",
        help="Run a comparison against the MIBiG dataset",
    )
    additional_group.add_argument(
        "--cb-general",
        action="store_true",
        help="Compare identified clusters against a database of antiSMASH-predicted clusters",
    )
    additional_group.add_argument(
        "--cb-subclusters",
        action="store_true",
        help="Compare identified clusters against known subclusters",
    )
    additional_group.add_argument(
        "--cb-knownclusters",
        action="store_true",
        help="Compare identified clusters against known gene clusters from the MIBiG database",
    )
    additional_group.add_argument(
        "--pfam2go",
        action="store_true",
        help="Run Pfam to Gene Ontology mapping module",
    )
    additional_group.add_argument(
        "--rre",
        action="store_true",
        help="Run RREFinder precision mode on all RiPP gene clusters",
    )
    additional_group.add_argument(
        "--smcog-trees",
        action="store_true",
        dest="smcog_trees",
        help="Generate phylogenetic trees of sec. met. cluster orthologous groups",
    )
    additional_group.add_argument(
        "--tfbs",
        action="store_true",
        help="Run TFBS finder on all gene clusters",
    )

    # HTML output options (HTML is always enabled in antiSMASH v8.0.0)
    html_group = parser.add_argument_group("HTML output")
    html_group.add_argument(
        "--html-description",
        default="antiSMASH results",
        help="Description for HTML output (HTML is always enabled)",
    )

    # Output format options
    format_group = parser.add_argument_group("Output format")
    format_group.add_argument(
        "--json", action="store_true", help="Output cluster information in JSON format"
    )
    format_group.add_argument(
        "--extract-clusters",
        action="store_true",
        help="Extract individual cluster GenBank files",
    )

    # Set the function to call
    parser.set_defaults(func=run_antismash_cli)


def run_antismash_cli(args):
    """
    Command-line interface for running antiSMASH.

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

        logger.info(f"Parsing pre-computed antiSMASH results: {args.parse}")

        # Parse antiSMASH results
        annotations_file = os.path.join(output_dir, "antismash.annotations.txt")
        clusters_file = os.path.join(output_dir, "antismash.clusters.txt")

        result = parse_antismash_gbk(
            input_file=args.parse,
            output_dir=output_dir,
            annotations_file=annotations_file,
            clusters_file=clusters_file,
        )
        backbone_domains, _, _, cluster_genes = result

        if backbone_domains is not None:
            logger.info(f"Parsed {len(cluster_genes)} clusters from {args.parse}")
            logger.info(f"Wrote annotations to {annotations_file}")
            logger.info(f"Wrote cluster information to {clusters_file}")

            # Extract individual cluster files if requested
            if args.extract_clusters:
                clusters_dir = os.path.join(output_dir, "clusters")
                cluster_files = extract_clusters(args.parse, clusters_dir)
                if cluster_files:
                    logger.info(
                        f"Extracted {len(cluster_files)} cluster files to {clusters_dir}"
                    )

            # Convert to JSON if requested
            if args.json:
                json_file = os.path.join(output_dir, "antismash.json")
                if antismash_to_json(args.parse, json_file):
                    logger.info(f"Wrote JSON cluster information to {json_file}")
                else:
                    logger.error(
                        f"Failed to write JSON cluster information to {json_file}"
                    )
        else:
            logger.error(f"Failed to parse clusters from {args.parse}")

        return

    # Check if at least one input option is provided
    if not args.input and not args.file:
        logger.error("Either --input, --file, or --parse must be specified")
        return

    # Get input file
    input_file = get_input_file(args, "genbank") if args.input else args.file

    if not input_file:
        logger.error("No GenBank file found")
        return

    # Log which GenBank file is being used for transparency
    if args.input:
        logger.info(
            f"Using GenBank file from funannotate2 predict directory: {input_file}"
        )
    else:
        logger.info(f"Using GenBank file: {input_file}")

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

    # Validate that the input file is a GenBank file
    if not input_file.endswith((".gbk", ".gbff", ".gb")):
        logger.error(
            f"Input file must be a GenBank file (.gbk, .gbff, or .gb): {input_file}"
        )
        return

    # Run antiSMASH
    gbk_file = run_antismash(
        input_file=input_file,
        output_dir=output_dir,
        antismash_path=args.antismash_path,
        cpu=args.cpus,
        taxon=args.taxon,
        fullhmmer=args.fullhmmer,
        cassis=args.cassis,
        clusterhmmer=args.clusterhmmer,
        tigrfam=args.tigrfam,
        asf=args.asf,
        cc_mibig=args.cc_mibig,
        cb_general=args.cb_general,
        cb_subclusters=args.cb_subclusters,
        cb_knownclusters=args.cb_knownclusters,
        pfam2go=args.pfam2go,
        rre=args.rre,
        smcog_trees=args.smcog_trees,
        tfbs=args.tfbs,
        html_description=args.html_description,
    )

    if gbk_file:
        # Parse antiSMASH results
        annotations_file = os.path.join(output_dir, "antismash.annotations.txt")
        clusters_file = os.path.join(output_dir, "antismash.clusters.txt")

        result = parse_antismash_gbk(
            input_file=gbk_file,
            output_dir=output_dir,
            annotations_file=annotations_file,
            clusters_file=clusters_file,
        )
        backbone_domains, _, _, cluster_genes = result

        if backbone_domains is not None:
            logger.info(f"Parsed {len(cluster_genes)} clusters from {gbk_file}")
            logger.info(f"Wrote annotations to {annotations_file}")
            logger.info(f"Wrote cluster information to {clusters_file}")

            # Extract individual cluster files if requested
            if args.extract_clusters:
                clusters_dir = os.path.join(output_dir, "clusters")
                cluster_files = extract_clusters(gbk_file, clusters_dir)
                if cluster_files:
                    logger.info(
                        f"Extracted {len(cluster_files)} cluster files to {clusters_dir}"
                    )

            # Convert to JSON if requested
            if args.json:
                json_file = os.path.join(output_dir, "antismash.json")
                if antismash_to_json(gbk_file, json_file):
                    logger.info(f"Wrote JSON cluster information to {json_file}")
                else:
                    logger.error(
                        f"Failed to write JSON cluster information to {json_file}"
                    )
        else:
            logger.error(f"Failed to parse clusters from {gbk_file}")
    else:
        logger.error("Failed to run antiSMASH")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run antiSMASH and parse results")
    subparsers = parser.add_subparsers(dest="command")
    antismash_subparser(subparsers)
    args = parser.parse_args()

    if args.command == "antismash":
        run_antismash_cli(args)
    else:
        parser.print_help()
