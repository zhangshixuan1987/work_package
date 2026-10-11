"""Fixed input locations for observations and land/grid support files.

Observations live in the shared archive ``/compyfs/zhan391/acme_init/Observations``
(one copy per product; the former ``v3_dart_cda_scratch/reference`` copy was
merged into it). Land masks and soil-layer thicknesses, and regrid grids/maps,
live next to it in ``acme_init``. Each location can be overridden with an
environment variable.
"""

import os
from pathlib import Path

OBSERVATION_ROOT = Path(os.environ.get("E3SM_DART_OBS_ROOT", "/compyfs/zhan391/acme_init/Observations"))
LAND_SUPPORT_DIR = Path(os.environ.get("E3SM_DART_LAND_REFERENCE_DIR", "/compyfs/zhan391/acme_init/lnd_sea_mask"))
REGRID_MAP_DIR = Path(os.environ.get("E3SM_DART_REGRID_MAP_DIR", "/compyfs/zhan391/acme_init/map_file"))

# Products stored below a sub-directory of their archive folder.
# ESA CCI holds several products; processed soil moisture (H2OSOI daily and
# monthly on 1 degree) is in ESA_CCI/H20SOI.
PRODUCT_DIRS = {"ESA_CCI": "ESA_CCI/H20SOI"}


def product_dir(product, root=OBSERVATION_ROOT):
    """Archive directory of an observation product, e.g. ``product_dir("ESA_CCI") / "daily"``."""
    return Path(root) / PRODUCT_DIRS.get(product, product)
