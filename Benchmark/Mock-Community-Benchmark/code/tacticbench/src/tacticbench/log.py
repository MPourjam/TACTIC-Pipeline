from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

_LOGGER: Optional[logging.Logger] = None


def get_logger(log_path: Optional[Path] = None) -> logging.Logger:
    """
    Create (once) and return a module-level logger.

    - Always logs INFO to stderr.
    - If log_path is provided, also logs INFO to that file.
    """
    global _LOGGER
    if _LOGGER is not None:
        return _LOGGER

    logger = logging.getLogger("tacticbench")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    # console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    logger.addHandler(ch)

    # optional file handler
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_path, encoding="utf-8")
        fh.setLevel(logging.INFO)
        fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        logger.addHandler(fh)

    _LOGGER = logger
    return logger
