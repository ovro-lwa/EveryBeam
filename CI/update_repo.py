#!/usr/bin/env python3

"""
Script to update a published Debian repository on a Debian package repository server.
"""

import argparse
import logging
import sys
from aptly_api import Client
from aptly_api.base import AptlyAPIException

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger()

def parse_arguments():
    """
    Parse command-line arguments.
    :return: argparse.Namespace
    """
    # Do not add help automatically to prevent that two sections of optional arguments appear
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description="Update a published Debian repository "
                                                 "on a Debian package repository server.",
                                     add_help=False)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("reponame",
                          help="name of the repository to update")
    required.add_argument("files", nargs="+",
                          help="white-space separated list of files to upload and publish")
    required.add_argument("-d", "--distribution", required=True,
                          help="name of the distribution (e.g. bionic)")
    # Add help manually to optional group.
    optional.add_argument("-h", "--help", action="help",
                          help="show this help message and exit")
    optional.add_argument("-p", "--prefix", default=None,
                          help="publishing prefix (default: <reponame>)")
    optional.add_argument("-s", "--skip-signing", action="store_true",
                          help="skip signing when publishing repository (not recommended)")
    optional.add_argument("-u", "--url", default="http://packages.astron.nl:8080",
                          help="URL of repository server (default: %(default)s)")
    optional.add_argument("--cleanup", action="store_true",
                          help="cleanup upload directory on the server")
    optional.add_argument("--gpg-key",
                          help="GPG key ID to use when signing the release")
    optional.add_argument("--gpg-passphrase",
                          help="GPG passphrase to unlock private key (possibly insecure)")
    optional.add_argument("--upload-dir", default=None,
                          help="upload directory on the server to temporarily store files "
                               "(default: <reponame>)")
    args = parser.parse_args()
    args.prefix = args.reponame if args.prefix is None else args.prefix
    args.upload_dir = args.reponame if args.upload_dir is None else args.upload_dir
    return args


def main(args):
    """
    Main function. It creates a new repository on the repository server, based on the command line
    arguments that were supplied.
    :param args: Arguments as parsed by `parse_arguments()`.
    """
    # Open connection to the server
    aptly = Client(args.url)
    status = 0

    try:
        # Upload file(s) to the server
        result = aptly.files.upload(args.upload_dir,
                                    *args.files)
        logger.info("Uploaded file(s): %s", " ".join(result))

        # Add uploaded files to the local repository; this will automatically remove them from the
        # upload directory
        result = aptly.repos.add_uploaded_file(reponame=args.reponame,
                                               force_replace=True,
                                               dir=args.upload_dir)
        added_files = result.report['Added']
        removed_files = result.report['Removed']

        if added_files or removed_files:
            added = "added: %s" % ", ".join(added_files) if added_files else ""
            removed = "removed: %s" % ", ".join(removed_files) if removed_files else ""
            logger.info("Updated local repository: %s", "; ".join((added, removed)))

        if result.failed_files:
            for warning in result.report['Warnings']:
                logger.warning(warning)
            logger.error("Updating file(s): %s", ", ".join(result.failed_files))
            status = 2

        # Update the published repository with the latest local content
        if args.skip_signing:
            result = aptly.publish.update(distribution=args.distribution,
                                        prefix=args.prefix,
                                        sign_skip=True)
        else:
            result = aptly.publish.update(distribution=args.distribution,
                                        prefix=args.prefix,
                                        sign_gpgkey=args.gpg_key,
                                        sign_passphrase=args.gpg_passphrase)
        logger.info("Published %s: %s", args.reponame, result)

    except AptlyAPIException as excp:
        logger.error(excp)
        status = 1

    finally:
        if args.cleanup:
            aptly.files.delete(args.upload_dir)

    return status


if __name__ == "__main__":
    sys.exit(main(parse_arguments()))
