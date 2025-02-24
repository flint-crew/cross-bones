from __future__ import annotations

import logging
import sys

handler = logging.StreamHandler(sys.stdout)
logging.getLogger().addHandler(handler)
logger = logging.getLogger("askapmetry")
logger.setLevel(logging.INFO)
