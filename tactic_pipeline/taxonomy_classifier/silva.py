"""
This module holds SILVA taxonomy classifier related code.
It contains functions and classes necessary to classify sequences based on
the SILVA taxonomy database.
"""
from __future__ import annotations
from pathlib import Path
from tactic_pipeline.taxonomy_classifier.classifier import TaxonomyClassifier
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union


class SILVATaxonomyClassifier(TaxonomyClassifier):
    """
    A concrete implementation of TaxonomyClassifier for the SILVA database.
    """

    db_dir: str = "/database/SILVA/LATEST/"
    db_download_base_url: str = "https://ftp.arb-silva.de/"

    def __init__(self, db_dir: str = None, db_version: str = "138_2"):
        """
        Initializes the SILVATaxonomyClassifier with the given database path.

        Args:
            db_dir (str, optional): Path to the SILVA taxonomy database directory.
                                    Defaults to the class-level db_dir.
        """
        super().__init__(db_dir)

    def __get_silva_versions(self) -> List[str]:
        """
        Returns a list of available SILVA database versions.

        Returns:
            List[str]: List of available SILVA database versions.
        """
        #TODO
        return ["138_2", "132", "128", "123"]

    @property
    def db_versions(self) -> str:
        """
        Returns the SILVA database versions.

        Returns:
            str: The SILVA database versions.
        """
        return self.__get_silva_versions()

    @property
    def db_latest_version_tag(self) -> str:
        """
        Returns the latest SILVA database version tag.

        Returns:
            str: The latest SILVA database version tag.
        """
        return self.db_versions[0]

    def download_db_files(self, url: str = None):
        """
        Downloads the SILVA taxonomy database files from the specified URL.

        Args:
            url (str, optional): The URL to download the database files from.
                                 Defaults to the class-level db_download_base_url.
        """
        #TODO
        pass

    def check_db_fasta_integrity(self) -> bool:
        """
        Checks the integrity of the SILVA taxonomy database FASTA file.

        Returns:
            bool: True if the FASTA file is intact, False otherwise.
        """
        # Implement SILVA-specific integrity check logic here
        fasta_path = Path(self.db_dir) / f"SILVA_{self.db_version}_SSURef_NR99_tax_silva.fasta"
        if not fasta_path.exists():
            return False
        # Additional integrity checks can be added here
        return True

    def check_db_arb_integrity(self) -> bool:
        """
        Checks the integrity of the SILVA taxonomy database ARB file.

        Returns:
            bool: True if the ARB file is intact, False otherwise.
        """
        # Implement SILVA-specific integrity check logic here
        arb_path = Path(self.db_dir) / "SILVA_138_SSURef_NR99_tax_silva.arb"
        if not arb_path.exists():
            return False
        # Additional integrity checks can be added here
        return True

    def check_db_idx_integrity(self) -> bool:
        """
        Checks the integrity of the SILVA taxonomy database index file.

        Returns:
            bool: True if the index file is intact, False otherwise.
        """
        # Implement SILVA-specific integrity check logic here
        idx_path = Path(self.db_dir) / "SILVA_138_SSURef_NR99_tax_silva.idx"
        if not idx_path.exists():
            return False
        # Additional integrity checks can be added here
        return True

    def check_db_integrity(self) -> bool:
        """
        Checks the integrity of the SILVA taxonomy database.

        Returns:
            bool: True if the database is intact, False otherwise.
        """
        # Implement SILVA-specific integrity check logic here
        expected_fiels = [
            "SILVA_138_SSURef_NR99_tax_silva.fasta",


