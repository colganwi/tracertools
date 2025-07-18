import argparse

import tracertools.scripts as scripts


def main():
    """Main entry point for the command line interface."""
    parser = argparse.ArgumentParser(prog="tracertools", description="Tracertools command line interface")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Add subcommands by importing the parsers from each script
    subparsers.add_parser(
        "alleles-from-bam", parents=[scripts.alleles_from_bam_script.get_parser()], help="Call alleles from bam file"
    )
    subparsers.add_parser(
        "count-integrations", parents=[scripts.count_integrations_script.get_parser()], help="Count integrations from allele counts files"
    )
    subparsers.add_parser(
        "edit-fractions", parents=[scripts.edit_fractions_script.get_parser()], help="Calculate edit fractions from allele counts files"
    )

    # Parse arguments and dispatch the function
    args = parser.parse_args()
    func = args.func
    kwargs = vars(args)
    func(**kwargs)


if __name__ == "__main__":
    main()
