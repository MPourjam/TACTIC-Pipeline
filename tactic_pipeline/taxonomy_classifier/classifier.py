from __future__ import annotations
import abc
import pickle
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union
import numpy as np


class TaxonomyClassifier(abc.ABC):

    db_dir: str = "/database/MyDatabase/LATEST/"
    db_download_base_url: str = ""

    def __init__(self, db_dir: str = None):
        """
        Initializes the TaxonomyClassifier with the given database path.

        Args:
            db_dir (str): Path to the taxonomy database directory.
        """
        if db_dir is None:
            db_dir = self.db_dir

    @abc.abstractmethod
    def check_db_integrity(self) -> bool:
        """
        Checks the integrity of the taxonomy database.

        Returns:
            bool: True if the database is intact, False otherwise.
        """
        # Placeholder for integrity check logic
        raise NotImplementedError("Database integrity check not implemented")
    
    @abc.abstractmethod
    def download_db_files(self, url: str):
        """
        Downloads the taxonomy database files from the specified URL.

        Args:
            url (str): The URL to download the database files from.
        """
        # Placeholder for download logic
        print(f"Downloading database files from {url}")
        raise NotImplementedError("Download database files not implemented")
    
    @abc.abstractmethod
    def index_db(self):
        """
        Indexes the taxonomy database for efficient querying.
        """
        # Placeholder for indexing logic
        print("Indexing database")
        raise NotImplementedError("Indexing database not implemented")

    @abc.abstractmethod
    def load_db(self, db_file: Union[str, Path] = None):
        """
        Loads the taxonomy database from the specified path.
        """
        # Placeholder for loading logic
        print(f"Loading database from {self.db_dir}")
        raise NotImplementedError("Loading database not implemented")

    @abc.abstractmethod
    def classify_sequences(self, sequences: List[str]) -> List[Tuple[str, str]]:
        """
        Classifies the given sequences using the taxonomy database.

        Args:
            sequences (List[str]): List of sequences to classify.

        Returns:
            List[Tuple[str, str]]: List of tuples containing sequence and its classification.
        """
        # Placeholder for classification logic
        print("Classifying sequences")
        raise NotImplementedError("Classifying sequences not implemented")

    @abc.abstractmethod
    def update_otu_table(self, otu_table_path: str, classifications: List[Tuple[str, str]], output_path: str):
        """
        Updates the OTU table with the given classifications.

        Args:
            otu_table_path (str): Path to the input OTU table.
            classifications (List[Tuple[str, str]]): List of sequence classifications.
            output_path (str): Path to write the updated OTU table.
        """
        # Placeholder for updating OTU table logic
        raise NotImplementedError("Updating OTU table not implemented")


    @abc.abstractmethod
    def update_otu_fasta(self, fasta_path: str, classifications: List[Tuple[str, str]], output_path: str):
        """
        Updates the OTU fasta file with the given classifications.

        Args:
            fasta_path (str): Path to the input fasta file.
            classifications (List[Tuple[str, str]]): List of sequence classifications.
            output_path (str): Path to write the updated fasta file.
        """
        # Placeholder for updating OTU fasta logic
        print(f"Updating OTU fasta at {fasta_path} and writing to {output_path}")
        raise NotImplementedError("Updating OTU fasta not implemented")
