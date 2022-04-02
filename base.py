import dataclasses
import collections
import operator
import numpy as np
from typing import List, Any


# Modification object
# name: name of the modification
# mass: mass of the modification
# residue: target residue
# position: Position of the modification in peptide sequence
Modification = collections.namedtuple(
    "Modification", ("name", "mass", "residue", "position"),
    defaults=["None"] * 4)

# XIC
XIC = collections.namedtuple("XIC", ["mz", "tol", "xic"])


@dataclasses.dataclass
class Assignment():
    """
    Assignment for each mass spectrum after database search
    """
    seq: str = None
    modifications: List = dataclasses.field(default_factory=list)
    charge: int = None
    delta: float = None
    missed_cleavage: int = None
    protein_accessions: List = dataclasses.field(default_factory=list)
    rank: int = None
    pep_type: str = None

    def get(self, attr):
        return self.__dict__.get(attr)


@dataclasses.dataclass
class SearchResult():
    """
    Results output by database search engine for each mass spectrum
    """
    data_id: str = None
    spectrum_id: Any = None
    mz: float = None
    nmatch: int = None
    matches: List = dataclasses.field(default_factory=list)

    def get(self, attr):
        return self.__dict__.get(attr)


@dataclasses.dataclass
class MassSpectrum():
    """
    Mass Spectrum object
    """
    scan: int = None
    spectrum: np.ndarray = None
    rt: float = None
    id: str = None
    precursormz: float = None
    charge: int = None
    collision: str = None
    energy: float = None
    ms_level: int = None

    def get(self, attr):
        return self.__dict__.get(attr)


@dataclasses.dataclass
class Localization():
    """
    Localization
    """
    score: float = None
    site_prob: float = None
    site: int = None
    mod: str = None
    target: str = None
