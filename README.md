# funannotate2-addons

Addon scripts for enhancing functional annotation in funannotate2 genome annotation pipeline.

## Standardized Output Format

All modules in funannotate2-addons produce a standardized 3-column TSV file that can be directly used with funannotate2 for functional annotation:

```
#gene_id	annotation_type	annotation_value
```

Where:
- **gene_id**: The identifier for the protein or gene
- **annotation_type**: The type of annotation (e.g., "name", "note", "go_term", "EC_number")
- **annotation_value**: The actual annotation content

This standardized format makes it easy to integrate the results with funannotate2's annotation system.

## Overview

funannotate2-addons provides tools for running and parsing results from external annotation tools:

- **EggNOG Mapper**: Functional annotation using orthology
- **InterProScan**: Protein domain and family annotation
- **antiSMASH**: Secondary metabolite gene cluster annotation
- **SignalP**: Signal peptide prediction

## Installation

```bash
pip install funannotate2-addons
```

## Usage

```bash
# Run
f2a <command> [options]
```

### Standardized Command-Line Interface

All funannotate2-addons commands use a consistent command-line interface:

```bash
# Basic usage pattern
f2a <command> -i /path/to/funannotate2_predict_folder [options]
```

Where:
- `-i, --input`: Path to a funannotate2 predict output directory
- `-o, --output`: (Optional) Output directory. If not specified, output will be written to `<input_dir>/annotate_misc`

This standardized interface makes it easy to run multiple annotation tools on the same funannotate2 predict output:

### EggNOG Mapper

EggNOG Mapper is a tool for fast functional annotation of novel sequences. This module provides a wrapper for running EggNOG Mapper and parsing its results.

```bash
# Basic usage
f2a emapper -i /path/to/funannotate2_predict_folder --cpu 8

# Specify output directory
f2a emapper -i /path/to/funannotate2_predict_folder -o /path/to/output_dir --cpu 8

# Specify database
f2a emapper -i /path/to/funannotate2_predict_folder --cpu 8 --database eggnog_proteins

# Export annotations to JSON
f2a emapper -i /path/to/funannotate2_predict_folder --cpu 8 --json
```

The module will generate:
1. A standard eggnog-mapper output file
2. An annotations file with gene names, descriptions, and EC numbers in the standardized 3-column format
3. Optional JSON file with annotation information

### InterProScan

InterProScan is a tool that scans protein sequences for matches against the InterPro protein signature databases. This module provides a wrapper for running InterProScan and parsing its results.

```bash
# Basic usage
f2a iprscan -i /path/to/funannotate2_predict_folder --cpu 8

# Specify output directory
f2a iprscan -i /path/to/funannotate2_predict_folder -o /path/to/output_dir --cpu 8

# Specify applications to run
f2a iprscan -i /path/to/funannotate2_predict_folder --cpu 8 --applications "Pfam,SMART,CDD"

# Disable GO term lookup
f2a iprscan -i /path/to/funannotate2_predict_folder --cpu 8 --no-goterms
```

### antiSMASH

antiSMASH is a tool for identifying and analyzing secondary metabolite biosynthesis gene clusters in bacterial and fungal genomes. This module provides a wrapper for running antiSMASH and parsing its results.

```bash
# Basic usage
f2a antismash -i /path/to/funannotate2_predict_folder --cpu 8

# Specify output directory
f2a antismash -i /path/to/funannotate2_predict_folder -o /path/to/output_dir --cpu 8

# Specify taxonomy
f2a antismash -i /path/to/funannotate2_predict_folder --cpu 8 --taxon fungi

# Extract individual cluster files
f2a antismash -i /path/to/funannotate2_predict_folder --cpu 8 --extract-clusters

# Export cluster information to JSON
f2a antismash -i /path/to/funannotate2_predict_folder --cpu 8 --json
```

The module will generate:
1. A GenBank file with annotated clusters
2. An annotations file with cluster and domain information in the standardized 3-column format
3. A clusters file with detailed information about each cluster
4. Optional individual GenBank files for each cluster
5. Optional JSON file with cluster information

**Note**: This module only accepts GenBank files from the funannotate2 predict output. Gene finding is disabled to ensure we use the gene models from the GenBank file rather than having antiSMASH attempt to predict genes, which doesn't work well for eukaryotic organisms.

### SignalP

SignalP is a tool for predicting the presence and location of signal peptides in proteins. This module provides a wrapper for running SignalP 6.0 and parsing its results.

```bash
# Basic usage
f2a signalp -i /path/to/funannotate2_predict_folder

# Specify output directory
f2a signalp -i /path/to/funannotate2_predict_folder -o /path/to/output_dir

# Specify organism group
f2a signalp -i /path/to/funannotate2_predict_folder --organism euk

# Generate plots
f2a signalp -i /path/to/funannotate2_predict_folder --plot

# Export predictions to JSON
f2a signalp -i /path/to/funannotate2_predict_folder --json
```

The module will generate:
1. A summary file with signal peptide predictions
2. An annotations file with signal peptide information in the standardized 3-column format
3. Optional plots for visualizing signal peptide predictions
4. Optional JSON file with prediction information

**Note**: This module is specifically designed for SignalP 6.0, which has a different command-line interface and output format compared to previous versions.

## Requirements

- Python 3.7+
- External tools must be installed separately:
  - EggNOG Mapper: https://github.com/eggnogdb/eggnog-mapper
  - InterProScan: https://interproscan-docs.readthedocs.io/en/latest/
  - antiSMASH: https://docs.antismash.secondarymetabolites.org/
  - SignalP: https://services.healthtech.dtu.dk/services/SignalP-6.0/

## Dependencies

funannotate2-addons has minimal dependencies:
- packaging
- requests
- gfftk
- gb-io (for GenBank file parsing)

Note: The antiSMASH module uses gb-io instead of BioPython to reduce dependencies.

## Development

### Testing

funannotate2-addons uses pytest for testing. To run the tests:

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run all tests
pytest

# Run tests with coverage report
pytest --cov=funannotate2_addons
```

Test files are located in the `tests/` directory.
