import csv
import dataclasses
import collections
from typing import Any

try:
    import xml.etree.cElementTree as et
except ImportError:
    import xml.etree.ElementTree as et

from .base import Modification, Assignment, Localization
from .base import SearchResult as _SearchResult


_SearchResultLocal = dataclasses.make_dataclass(
    "_SearchResultLocal",
    [(_fd, _fd_prop.type, _fd_prop.default)
     for _fd, _fd_prop in _SearchResult.__dataclass_fields__.items()
     if _fd != "matches"])

_TABLE_HEADS = [("seq", "Sequence", lambda x: x),
                ("charge", "Charge", lambda x: x if x is None else int(x)),
                ("delta", "Mass error [Da]",
                 lambda x: x if x is None else float(x)),
                ("missed_cleavage", "Missed cleavages", int),
                ("protein_accessions", "Proteins", lambda x: x.split(";")),
                ("pep_type", "Reverse",
                 lambda x: "normal" if not x else "decoy"),
                ("data_id", "Raw file", lambda x: x),
                ("spectrum_id", "Scan number", int),
                ("mz", "m/z", float),
                ("rt", "Retention time", float),
                ("pep", "PEP", float),
                ("score", "Score", float),
                ("delta_score", "Delta score", float)]


@dataclasses.dataclass
class SearchResult(Assignment, _SearchResultLocal):
    localization: Any = None
    rt: float = None
    pep: float = None
    score: float = None
    delta_score: float = None


def _get_bracket_index(string):
    """
    Parse the seq with brackets to split the string and get
    the indices of brackets in the sequence.
    """
    if "(" not in string:
        return None, None

    br_index, stack = [], []
    for i, s in enumerate(string):
        if s == "(":
            if not stack:
                br_index.append(i)
            stack.append(i)
        elif s == ")":
            _ = stack.pop()
            if not stack:
                br_index.append(i)

    splitted_str = [
        string[i+1:j] for i, j in zip(br_index, br_index[1:] + [None])
    ]
    if br_index[0] != 0:
        splitted_str.insert(0, string[:br_index[0]])

    # parse the modified peptides
    n = len(splitted_str)
    # starting index of sequence
    i0 = 0 if br_index[0] > 0 else 1
    # starting index of modification
    r0 = (i0 + 1) % 2
    sites = []
    for i in range(r0, n, 2):
        seq = "".join(splitted_str[i0:i:2])
        sites.append((splitted_str[i], len(seq)))

    return sites, "".join(splitted_str[i0::2])


class AndromedaParser():
    """
    Parse Andromeda search results, i.e., msms.txt file
    """
    def __init__(self):
        pass

    def get_results(self, search_results, param_file=None,
                    reader=csv.DictReader, header=_TABLE_HEADS):
        """
        Get search results
        """
        fixed_mods = self._get_param_fixed_mod(param_file)

        # get search results
        res = []
        with open(search_results, "r") as f:
            f_table = reader(f, delimiter="\t")
            # get headers
            self._headers = {fd.lower(): fd for fd in f_table.fieldnames}
            for cell in f_table:
                # parse modification and get the site probabilities
                mods_cp = self._split_mod(cell.get("Modifications"))
                mods, locs = [], []
                if mods_cp is not None:
                    mod_pep = cell.get("Modified sequence")
                    mods = self._parse_mod_peptide(mod_pep, mods_cp)

                    # localization
                    locs = self._get_localizations(cell, mods_cp)
                mods = self._update_fix_mod(mods, fixed_mods,
                                             cell.get("Sequence"))

                res.append(
                    SearchResult(modifications=mods, localization=locs,
                                 **{fd: convert(cell.get(head))
                                    for fd, head, convert in header}))

        return res

    def _get_param_fixed_mod(self, param_file):
        """
        Parse parameter file.
        """
        if param_file is None:
            return [("Carbamidomethyl", "C", None)]

        fixed_mods = []
        for _, el in et.iterparse(param_file):
            if el.tag.endswith("fixedModifications"):
                for ch in el:
                    fixed_mods.append(self._parse_mod_type(ch.text))

        if not fixed_mods:
            fixed_mods = None

        return fixed_mods

    def _split_mod(self, string):
        """
        Split modification string to get modification composition
        """
        if string == "Unmodified":
            return None

        mods = collections.defaultdict()
        for mod_type in string.split(","):
            if "(" not in mod_type:
                mod_head = mod_type.split()[-1]
            else:
                mod_head = " ".join(mod_type.split()[-2:])
            mods[mod_head] = self._parse_mod_type(mod_head)

        return mods

    def _parse_mod_peptide(self, string, mods):
        """
        Parse modified peptides
        """
        # get nested brackets from the modified peptide
        mods_sites, seq = _get_bracket_index(string[1:-1])
        # get modifications
        if mods_sites is None:
            return []

        pep_mods = []
        if mods_sites is not None:
            for name, i in mods_sites:
                r = seq[i - 1]
                if name in mods:
                    mod, _, _ = mods[name]
                else:
                    mod = [m for _, (m, mr, _) in mods.items() if r in mr][0]
                pep_mods.append(Modification(
                    name=mod, mass=None, residue=seq[i-1], position=i
                ))

        return pep_mods

    def _update_fix_mod(self, mods, fixed_mods, seq):
        """
        Update modifications if fixed modifications are set.
        """
        # get fixed modification
        if fixed_mods is None:
            return mods

        for name, res, t in fixed_mods:
            if t is None and res is not None:
                mods += [Modification(name=name, mass=None,
                                      residue=s, position=i+1)
                         for i, s in enumerate(seq) if s in res]
            elif t is None:
                mods.append(Modification(name=name, mass=None,
                                         residue=None, position=t))
        return mods

    def _get_localizations(self, cell, mods,
                           replace=dataclasses.replace):
        """
        Get the localizations
        """
        # localizations
        locates = []
        # get site localizations
        for m, (name, res, t) in mods.items():
            probs, _ = _get_bracket_index(
                cell.get(self._headers.get(f"{m.lower()} probabilities")))
            scores, _ = _get_bracket_index(
                cell.get(self._headers.get(f"{m.lower()} score diffs")))

            # store them
            locates += [Localization(score=float(s), site_prob=float(p),
                                      site=i0, mod=name, target=res)
                         for (i0, p), (_, s) in zip(probs, scores)]
        return locates

    def _parse_mod_type(self, mod_type):
        """
        Parse type of modification
        """
        if "(" not in mod_type:
            return mod_type, None, None

        name, res = tuple(mod_type.split())
        return ((name, None, "nterm") if res.endswith("N-term") else
                (name, None, "cterm") if res.endswith("C-term") else
                (name, res, None))
