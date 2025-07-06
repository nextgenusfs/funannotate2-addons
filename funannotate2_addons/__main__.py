#!/usr/bin/env python3

import sys
import os
import argparse
from .__init__ import __version__
from .help_formatter import MyParser, MyHelpFormatter
from .emapper import emapper_subparser, run_emapper_cli
from .iprscan import iprscan_subparser, run_iprscan_cli
from .antismash import antismash_subparser, run_antismash_cli
from .signalp6 import signalp_subparser, run_signalp_cli


def main():
    args = parse_args(sys.argv[1:])

    if args.subparser_name == "emapper":
        run_emapper_cli(args)
    elif args.subparser_name == "iprscan":
        run_iprscan_cli(args)
    elif args.subparser_name == "antismash":
        run_antismash_cli(args)
    elif args.subparser_name == "signalp":
        run_signalp_cli(args)
    # elif args.subparser_name == 'update':
    #    update(args)


def parse_args(args):
    description = "Funannotate2-addons: eukaryotic genome annotation scripts"
    parser = MyParser(
        description=description, formatter_class=MyHelpFormatter, add_help=False
    )
    subparsers = parser.add_subparsers(title="Commands", dest="subparser_name")
    # add subtools here
    emapper_subparser(subparsers)
    iprscan_subparser(subparsers)
    antismash_subparser(subparsers)
    signalp_subparser(subparsers)

    # update_subparser(subparsers)

    help_args = parser.add_argument_group("Help")
    help_args.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit",
    )
    help_args.add_argument(
        "--version",
        action="version",
        version="{} v{}".format(
            os.path.basename(os.path.dirname(os.path.realpath(__file__))), __version__
        ),
        help="show program's version number and exit",
    )

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        raise SystemExit(0)

    return parser.parse_args(args)


if __name__ == "__main__":
    main()
