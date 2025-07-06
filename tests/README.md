# Tests for funannotate2-addons

This directory contains tests for the funannotate2-addons package.

## Running Tests

You can run the tests using pytest:

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run all tests
pytest

# Run tests with coverage report
pytest --cov=funannotate2_addons

# Run specific test file
pytest tests/unit/test_emapper.py
```

## Test Structure

- `unit/`: Unit tests for individual modules
  - `test_emapper.py`: Tests for the EggNOG Mapper module
  - `test_iprscan.py`: Tests for the InterProScan module
  - `test_antismash.py`: Tests for the antiSMASH module
  - `test_signalp.py`: Tests for the SignalP module

## Adding New Tests

When adding new functionality, please also add corresponding tests. Follow the existing test structure and naming conventions.
