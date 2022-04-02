import re
import collections
import dataclasses
from typing import List, Any

from .base import Modification, Assignment
from .base import SearchResult as _SearchResult

try:
    from lxml import etree as et
except ImportError:
    try:
        import xml.etree.cElementTree as et
    except ImportError:
        import xml.etree.ElementTree as et


ITRAQ = collections.namedtuple("ITRAQ", ["tag", "area", "error"])

_Tags = ["SPECTRUM", "MATCH", "MOD_FEATURE", "TERM_MOD_FEATURE",
         "ITRAQPEAKS", "STATISTIC", "SEQUENCE", "PROTEIN_CONTEXT"]

# because possible special characters in number, will check
Non_decimal_match = re.compile(r"[^\d.]+")


@dataclasses.dataclass
class Match(Assignment):
    score: float = None
    evalue: float = None
    confidence: float = None
    mod_prob: float = None


@dataclasses.dataclass
class SearchResult(_SearchResult):
    rt: float = None
    itraq: Any = None


class ProteinPilotReader():
    """
    Read ProteinPilot search results
    """

    _Match_fields = [("score", "score", float),
                     ("evalue", "eval", float),
                     ("confidence", "confidence", float),
                     ("seq", "seq", lambda x: x),
                     ("charge", "charge", int),
                     ("delta", "da_delta", float),
                     ("pep_type", "type", int),
                     ("mod_prob", "mod_prob", float)]

    def __init__(self):
        pass

    def get_results(self, filename, tags=_Tags,
                    replace=dataclasses.replace):
        """
        Get search results
        """
        results, accessions = [], collections.defaultdict()
        matches, num_hypo, mods, accx = [], [], [], []
        for _, _elem in et.iterparse(filename, recover=True, encoding="iso-8859-1"):
            _tag = _elem.tag
            if _tag in tags:

                if _tag == "MOD_FEATURE" or _tag == "TERM_MOD_FEATURE":
                    mod = self._parse_modification(_elem)
                    if mod is not None:
                        mods.append(mod)

                elif _tag == "MATCH":
                    match = self._parse_match(_elem)
                    _t = "normal" if match.pep_type == 0 else "decoy"
                    if mods:
                        mods = [m._replace(residue=match.seq[m.position-1])
                                if isinstance(m.position, int) else m
                                for m in mods]
                    match = replace(match, modifications=mods, pep_type=_t,
                                    rank=len(matches) + 1)
                    matches.append(match)

                    # initialize for next match
                    mods = []

                elif _tag == "ITRAQPEAKS":
                    itraq = self._parse_peaks(_elem)

                elif _tag == "STATISTIC":
                    num_hypo.append(self._parse_statistic(_elem))

                elif _tag == "SPECTRUM":
                    if matches:
                        spectrum = self._parse_spectrum_head(_elem)
                        spectrum = replace(spectrum, itraq=itraq,
                                           nmatch=num_hypo, matches=matches)
                        results.append(spectrum)

                        # initialize for next spectrum
                        matches, num_hypo = [], []

                elif _tag == "SEQUENCE":
                    # get peptide sequence for protein accession
                    # assignments
                    seq = self._parse_sequence(_elem)
                    accessions[seq] = accx

                    # initialize for next sequence
                    accx = []

                elif _tag == "PROTEIN_CONTEXT":
                    # get protein accessions
                    accx.append(self._parse_protein_context(_elem))

            _elem.clear()

        return self._merge_results(results, accessions)

    def _parse_spectrum_head(self, element):
        """
        Parse the spectrum head
        """
        _id = element.keys()[-2]
        spectrum = {"rt": float(element.get("elution")),
                    "mz": float(element.get("precursormass")),
                    "spectrum_id": ".".join(element.get(_id).split(".")[:-1])}

        return SearchResult(**spectrum)

    def _parse_match(self, element):
        """
        Parse the peptide spectrum matches
        """
        return Match(**{_field: _convert(element.get(_key))
                        for _field, _key, _convert in self._Match_fields})

    def _parse_modification(self, element):
        """
        Parse the modifications
        """
        name, pos = element.get("mod"), element.get("pos")
        if name.startswith("No ") or pos is None:
            return None

        if "&gt;" in name:
            name = name.replace("&gt;", ">")

        if name.startswith("Protein N-Terminal"):
            name = name.lstrip("Protein N-Terminal")

        if name.startswith("Protein C-Terminal"):
            name = name.lstrip("Protein C-Terminal")

        pos = int(pos)
        if element.tag == "TERM_MOD_FEATURE":
            pos = "nterm" if pos == 1 else "cterm"

        return Modification(name=name, mass=None, residue=None, position=pos)

    def _parse_peaks(self, element, ITRAQ=ITRAQ, sub=Non_decimal_match.sub):
        """
        Parse MS peak
        """
        if element.text is None:
            return None

        try:
            itraqs = [ITRAQ(*(float(x) for x in peak.split())) for peak in
                      element.text.splitlines() if peak]
        except ValueError:
            itraqs = [ITRAQ(*(float(sub("", x)) for x in peak.split()))
                      for peak in element.text.splitlines() if peak]

        return itraqs

    def _parse_statistic(self, element):
        """
        Parse Statistc to get number of hypotheses
        """
        return int(element.get("numberofhypotheses"))

    def _parse_sequence(self, element):
        """
        Parse Sequence to get protein accessions the peptide belongs to
        """
        return element.get(element.keys()[0])

    def _parse_protein_context(self, element):
        """
        Protein accessions
        """
        return element.get("protein")

    def _merge_results(self, results, accessions, replace=dataclasses.replace):
        """
        Merge the identifications and protein accessions
        """
        if not accessions:
            return results

        for i, spectrum in enumerate(results):
            matches = [m if accessions.get(m.seq) is None else
                       replace(m, protein_accessions=accessions.get(m.seq))
                       for m in spectrum.matches]
            results[i] = replace(spectrum, matches=matches)
        return results
