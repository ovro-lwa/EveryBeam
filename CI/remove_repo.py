#!/usr/bin/env python3

"""
Script to un-publish and remove a Debian repository from a Debian package repository server.
"""

import argparse
import logging
import sys
from aptly_api import Client
from aptly_api.base import AptlyAPIException

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger()


def parse_arguments():
    """
    Parse command-line arguments.
    :return: argparse.Namespace
    """
    # Do not add help automatically to prevent that two sections of optional arguments appear
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description="Un-publish and remove a Debian repository "
                                                 "from a Debian package repository server.",
                                     add_help=False)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("reponame",
                          help="name of the repository to remove")
    required.add_argument("-d", "--distribution", required=True,
                          help="distribution of the repository (e.g. bionic)")
    optional.add_argument("-h", "--help", action="help",
                          help="show this help message and exit")
    optional.add_argument("-p", "--prefix", default=None,
                          help="prefix of the repository (default: <reponame>")
    optional.add_argument("-u", "--url", default="http://packages.astron.nl:8080",
                          help="URL of repository server (default: %(default)s)")
    args = parser.parse_args()
    args.prefix = args.reponame if args.prefix is None else args.prefix
    return args


def main(args):
    """
    Main function. It creates a new repository on the repository server, based on the command line
    arguments that were supplied.
    :param args: Arguments as parsed by `parse_arguments()`.
    :return: Error status: 0: success; 1: error; 2: warning
    """
    # Open connection to the server
    aptly = Client(args.url)
    status = 0

    # Un-publish the repository
    try:
        aptly.publish.drop(prefix=args.prefix,
                        distribution=args.distribution)
        logger.info("Un-published repository %s", args.reponame)
    except AptlyAPIException as excp:
        logger.warning(excp)
        status = 2

    # Remove the repository from the server
    try:
        aptly.repos.delete(args.reponame)
        logger.info("Deleted repository %s", args.reponame)
    except AptlyAPIException as excp:
        logger.error(excp)
        status = 1

    return status


if __name__ == "__main__":
    sys.exit(main(parse_arguments()))
