#!/usr/bin/env python3

"""
Script to create and publish a Debian repository on a Debian package repository server.
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
                                     description="Create and publish a Debian repository "
                                                 "on a Debian package repository server.",
                                     add_help=False)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("reponame",
                          help="name of the repository to create")
    required.add_argument("-a", "--architectures", nargs="+", required=True,
                          help="white-space separated list of supported architectures (e.g. amd64)")
    required.add_argument("-d", "--distribution", required=True,
                          help="default distribution when publishing (e.g. bionic)")
    optional.add_argument("-c", "--component",
                          help="default component when publishing (e.g. main)")
    # Add help manually to optional group.
    optional.add_argument("-h", "--help", action="help",
                          help="show this help message and exit")
    optional.add_argument("-p", "--prefix", default=None,
                          help="publishing prefix (default: <reponame>)")
    optional.add_argument("-r", "--rollback", action="store_true",
                          help="roll back all changes upon failure")
    optional.add_argument("-s", "--skip-signing", action="store_true",
                          help="skip signing when publishing repository (not recommended)")
    optional.add_argument("-u", "--url", default="http://packages.astron.nl:8080",
                          help="URL of repository server (default: %(default)s)")
    optional.add_argument("--gpg-key",
                          help="GPG key ID to use when signing the release")
    optional.add_argument("--gpg-passphrase",
                          help="GPG passphrase to unlock private key (possibly insecure)")
    args = parser.parse_args()
    args.prefix = args.reponame if args.prefix is None else args.prefix
    return args


def main(args):
    """
    Main function. It creates a new repository on the repository server, based on the command line
    arguments that were supplied.
    :param args: Arguments as parsed by `parse_arguments()`.
    :return: Exit status: 0: success; 1: error
    """
    # Open connection to the server
    aptly = Client(args.url)

    # Create repository on the server
    try:
        result = aptly.repos.create(args.reponame,
                                    default_component=args.component,
                                    default_distribution=args.distribution)
        logger.info("Created repository %s: %s", args.reponame, result)
    except AptlyAPIException as excp:
        logger.error("Failed to create repository %s: %s", args.reponame, excp)
        return 1

    # Publish the repository on the server
    try:
        if args.skip_signing:
            result = aptly.publish.publish(sources=[{'Name': args.reponame}],
                                           architectures=args.architectures,
                                           distribution=args.distribution,
                                           prefix=args.prefix,
                                           sign_skip=True)
        else:
            result = aptly.publish.publish(sources=[{'Name': args.reponame}],
                                           architectures=args.architectures,
                                           distribution=args.distribution,
                                           prefix=args.prefix,
                                           sign_gpgkey=args.gpg_key,
                                           sign_passphrase=args.gpg_passphrase,
                                           sign_batch=True)
        logger.info("Published repository %s: %s", args.reponame, result)
    except AptlyAPIException as excp:
        logger.error("Failed to publish repository %s: %s", args.reponame, excp)
        # Roll back any changes if requested by the user
        if args.rollback:
            try:
                aptly.repos.delete(args.reponame)
                logger.info("Rolled back changes")
            except AptlyAPIException as excp:
                logger.error("Failed to rollback changes: %s", excp)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main(parse_arguments()))
