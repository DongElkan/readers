"""
Parse Mascot search results (.dat)
"""
import re
import collections
import urllib
import dataclasses

from .base import Modification, Assignment
from .base import SearchResult as _SearchResult


_param_fields = ["ftol", "ftol_u", "ptol", "ptol_u", "decoy", "database",
                 "cleavage", "err_tol", "charges", "version", "quantitation",
                 "total_queries"]
Parameters = collections.namedtuple("Parameters", _param_fields,
                                    defaults=(None,) * len(_param_fields))

Mod_param_pattern = re.compile("([A-Za-z]+)([0-9]+)=(-?[0-9.]+),(.*) \((.*)\)")


@dataclasses.dataclass
class Match(Assignment):
    ionscores: float = None


@dataclasses.dataclass
class SearchResult(_SearchResult):
    max_mz: float = None
    min_mz: float = None
    scan: int = None
    query: int = None


class MascotdatParser:
    """
    Get results from Mascot dat files
    """
    _ParamNamePairs = {"ITOL": "ftol", "ITOLU": "ftol_u", "TOL": "ptol",
                       "TOLU": "ptol_u", "ERRORTOLERANT": "err_tol",
                       "DB": "database", "CLE": "cleavage", "decoy": "DECOY",
                       "CHARGE": "charges", "FORMVER": "version",
                       "quantitation": "QUANTITATION"}

    _MatchNames = ("seq", "ionscores", "delta", "missed_cleavage", "charge",
                   "modifications", "rank", "pep_type", "protein_accessions")

    _ConvertSpectrumItem = {"title": ("spectrum_id", urllib.parse.unquote),
                            "scans": ("scan", int),
                            "mass_min": ("min_mz", float),
                            "mass_max": ("max_mz", float)}

    _AttributeforClear = ["res", "_charges", "_alt_charges", "_curr_match",
                          "_curr_spectrum", "fixed_mods", "variable_mods"]

    def __init__(self):
        # set up parameters
        self.fixed_mods = collections.defaultdict()
        self.variable_mods = collections.defaultdict()
        self.search_parameters = Parameters()
        self._parser = {"parameters": self._parse_params,
                        "masses": self._parse_mod_param,
                        "header": self._parse_header,
                        "summary": self._parse_summary,
                        "et_summary": self._parse_summary,
                        "decoy_summary": self._parse_summary,
                        "peptides": self._parse_peptide,
                        "et_peptides": self._parse_peptide,
                        "decoy_peptides": self._parse_peptide,
                        "query": self._parse_spectrum}

        # initialization
        self.res = collections.defaultdict()
        self._charges = []
        self._alt_charges = []
        self._curr_match = collections.defaultdict(list)
        self._curr_spectrum = collections.defaultdict(list)

    def _clear(self):
        """Clear the contents."""
        for _attribute in self._AttributeforClear:
            self.__dict__.get(_attribute).clear()

    def get_results(self, result_file):
        """Get Mascot search results."""
        skip_read = True
        with open(result_file, "r") as f:
            for line in f:
                tline = line.rstrip()
                if "=" not in tline:
                    continue

                # headers
                if tline.startswith('Content-Type'):
                    content = re.findall('"(.*)"', tline)
                    content = content[0] if content else None
                    if (not content or not (content in self._parser or
                                            content.startswith("query"))):
                        skip_read = True
                        continue
                    else:
                        skip_read = False

                    parser = (self._parse_spectrum
                              if content.startswith("query")
                              else self._parser[content])

                    if (content.endswith("peptides")
                            or content.endswith("summary")):
                        self._pep_type = (
                            "err_tol" if content.startswith("et") else
                            "decoy" if content.startswith("decoy") else
                            "normal"
                        )

                    self.content = content
                else:
                    if skip_read or not tline:
                        continue
                    _ = parser(tline)

        # the last records
        self._parse_last()

        res = [psm for psm in self.res.values() if psm.matches]
        # release memories
        self._clear()

        return res

    def _parse_last(self):
        """ Parse last record """
        # the last spectrum record
        if self._curr_spectrum:
            self._parse_spectrum_info(self._curr_spectrum)
        if self._curr_match:
            _ = self._parse_peptide_records(self._curr_match)

    def _parse_params(self, string):
        """Parse parameters."""
        param_name, param_val = tuple(string.split("="))
        if param_name in self._ParamNamePairs:
            field = self._ParamNamePairs[param_name]
            if not param_val:
                param_val = None
            self.search_parameters = self.search_parameters._replace(
                **{field: param_val}
            )

        return self

    def _parse_mod_param(self, string):
        """Parse modification string."""
        if string.startswith("delta") or string.startswith("FixedMod"):
            p = Mod_param_pattern.match(string)
            if p is not None:
                t, no, m, s, r = p.groups()
                param_vals = {"mass": float(m), "name": s}

                if r.endswith("N-term"):
                    param_vals.update({"residues": None, "position": "nterm"})
                elif r.endswith("C-term"):
                    param_vals.update({"residues": None, "position": "cterm"})
                else:
                    param_vals.update({"residues": r, "position": None})

                if t == "delta":
                    self.variable_mods[no] = param_vals
                elif t == "FixedMod":
                    self.fixed_mods[no] = param_vals

        return self

    def _parse_header(self, string):
        """Get the parameters in header subsection."""
        if string.startswith("queries"):
            total_queries = int(string.split("=")[1])
            self.search_parameters = self.search_parameters._replace(
                total_queries=total_queries
            )

        return self

    def _parse_summary(self, string, SearchResult=SearchResult):
        """Get m/z, charges of mass spectra and number of matches."""
        if self._pep_type == "normal":
            field, val_str = tuple(string.split("=")[:2])
            parsed = re.match("([a-z]+)([0-9]+)", field)
            if parsed is None:
                return self
            s, n = parsed.groups()
            n = int(n)

            if s == "qmass":
                # initialization of new mass spectrum search
                self.res[n] = SearchResult(query=n)
            elif s == "qexp":
                mz, c = re.match("(.*),([0-9]+)", val_str).groups()
                self.res[n].mz = float(mz)
                self._charges.append(int(c))
            elif s == "qmatch":
                self.res[n].nmatch = int(val_str)

        elif string.startswith("qexp"):
            # charge state of alternative identification types,
            # i.e., decoy or error tolerant search
            _, c = re.match("(.*),([0-9]+)", string.split("=")[1]).groups()
            self._alt_charges.append(int(c))

        return self

    def _parse_peptide(self, string):
        """Get peptide matches."""
        if string.endswith("-1"):
            return self

        id_info, _ = tuple(string.split("="))

        # query number and rank
        n, r = re.match("q([0-9]+)_p([0-9]+)", id_info).groups()
        if self._curr_match and n != next(iter(self._curr_match)):
            _ = self._parse_peptide_records(self._curr_match)
            self._curr_match.clear()

        self._curr_match[n].append((int(r), id_info.endswith(r), string))
        return self

    def _parse_peptide_records(self, records):
        """Parse records for the matches to current mass spectrum."""
        n, records = next(iter(records.items()))
        n = int(n) - 1
        _type = self._pep_type
        c = self._charges[n] if _type == "normal" else self._alt_charges[n]

        # matches
        matches = []
        for r, is_pep, record in records:
            id_info, match_info = tuple(record.split("="))

            if is_pep:
                pep_info, pro_info = tuple(match_info.split(";"))
                # get peptide matches and protein accessions
                (nmc, _, delta, _, seq, _, modstr, ionscores, _) = tuple(
                    pep_info.split(",")[:9]
                )

                accessions = tuple(re.findall('"([^"]*)"', pro_info))
                # parse modifications
                mods = self._parse_modstr(modstr, seq)
                # to assignments
                match = zip(self._MatchNames,
                            [seq, float(ionscores), float(delta),
                             int(nmc), c, mods, r, _type, accessions])
                matches.append(Match(**dict(match)))

            elif _type == "err_tol" and id_info.endswith("et_mods"):
                # update error tolerant search modifications
                m, _, name, _ = self._parse_et_mods(match_info)
                new_mods = [mod._replace(name=name, mass=m)
                            if mod.name == "X" else mod
                            for mod in matches[r - 1].modifications]
                matches[r - 1].modifications = new_mods

        self.res[n + 1].matches += matches

        return self

    def _parse_spectrum(self, string):
        """Get title."""
        if self._curr_spectrum and self.content not in self._curr_spectrum:
            self._parse_spectrum_info(self._curr_spectrum)
            self._curr_spectrum.clear()

        self._curr_spectrum[self.content].append(string)
        return self

    def _parse_spectrum_info(self, spectrum_info):
        """Parse spectrum info."""
        content, strings = next(iter(spectrum_info.items()))
        n = int(content.lstrip("query"))
        spectrum_info = []
        for string in strings:
            name, content = tuple(string.split("="))
            if name in self._ConvertSpectrumItem:
                item_name, convert = self._ConvertSpectrumItem.get(name)
                spectrum_info.append((item_name, convert(content)))

        match = self.res[n]
        self.res[n] = dataclasses.replace(match, **dict(spectrum_info))

    def _parse_modstr(self, string, seq, Modification=Modification):
        """Parse modification string."""
        mods = []
        if self.fixed_mods is not None:
            # fixed modifications
            for m in self.fixed_mods.values():
                name, mass = m.get("name"), m.get("mass")
                pos, res = m.get("position"), m.get("residues")
                if pos is not None:
                    mods.append(Modification(
                        name=name, mass=mass, residue=None, position=pos
                    ))
                else:
                    mods += [Modification(name=name, mass=mass,
                                          residue=r, position=i + 1)
                             for i, r in enumerate(seq) if r in res]

        # variable modifications
        n = len(string) - 1
        mod_tags = [(i, x) for i, x in enumerate(string) if x != "0"]
        if not mod_tags:
            return mods

        for i, _s in mod_tags:
            pos = i if n > i > 0 else "nterm" if i == 0 else "cterm"
            res = None if pos in ["nterm", "cterm"] else seq[i - 1]
            if _s == "X":
                name, mass = "X", None
            else:
                name = self.variable_mods[_s].get("name")
                mass = self.variable_mods[_s].get("mass")
            mods.append(Modification(
                name=name, mass=mass, residue=res, position=pos
            ))

        return mods

    def _parse_et_mods(self, string):
        """Get error tolerance modifications."""
        return re.match(r"([-+]?[0-9.]+),([-+]?[0-9.]+),(.*) \((.*)\)",
                        string).groups()
