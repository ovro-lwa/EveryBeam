# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import pytest
from subprocess import call


@pytest.fixture(scope="session", autouse=True)
def download_msets():
    """
    Download LBA, HBA and OSKAR measurement sets if they are not yet found.
    Fixture is run only once at start of session.
    """
    download_scripts = [
        "download_lofar_lba_ms.sh",
        "download_lofar_hba_ms.sh",
        "download_oskar_ms.sh",
        "download_ska_mid_ms.sh",
    ]
    for script in download_scripts:
        call(
            f"sh {os.path.join(os.environ['SCRIPTS_DIR'], script)}", shell=True
        )
