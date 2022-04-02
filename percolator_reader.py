import csv
import os
import collections
import warnings
import re
import dataclasses

from .comet_reader import Match, SearchResult


@dataclasses.dataclass
class Percolator:
    seq: str = None
    mod_seq: str = None
    charge: int = None
    mz: float = None
    qvalue: float = None
    pl_score: float = None
    pep: float = None


_SearchResult = dataclasses.make_dataclass(
    "_SearchResult", 
    [(_fd, _fd_prop.type, _fd_prop.default)
     for _fd, _fd_prop in SearchResult.__dataclass_fields__.items()
     if _fd != "matches"])


@dataclasses.dataclass
class Result(Match, Percolator, _SearchResult):
    pass


PercolatorFiles = collections.namedtuple(
    "PercolatorFiles", ("target", "decoy", "log", "path")
)

Pattern = re.compile("([A-Z]+)")

Pick_up_fields = [("charge", "charge", int),
                  ("spectrum precursor m/z", "mz", float),
                  ("percolator score", "pl_score", float),
                  ("percolator q-value", "qvalue", float),
                  ("percolator PEP", "pep", float),
                  ("sequence", "seq", lambda s: "".join(Pattern.findall(s))),
                  ("sequence", "mod_seq", lambda s: s)]


class PercolatorReader():
    """
    Read post-processing results output by Percolator
    """
    def __init__(self, result_file):
        self.result = result_file
        self._collect_processing_results()
        self._file_index = self._parse_file_index()

    def get_results(self):
        """
        Get the results from "results_file"
        """
        validation_results = self._read_results(
            os.path.join(self._percolator_files.path,
                         self._percolator_files.target)
        )

        if self._percolator_files.decoy is not None:
            val_res_decoy = self._read_results(
                os.path.join(self._percolator_files.path,
                             self._percolator_files.decoy)
            )
            # merge targets and decoys together
            for file_ in val_res_decoy:
                if file_ in validation_results:
                    validation_results[file_].update(val_res_decoy[file_])
                else:
                    validation_results[file_] = val_res_decoy[file_]

        return validation_results

    def merge_results(self, validation_results, search_results):
        """
        Merge the validation results to search results
        """
        self._is_compatible()

        search_file = None
        results = []
        for spectrum in search_results:
            if search_file is None:
                search_file = spectrum.spectrum_id.split(".")[0]
            pl_res = validation_results[search_file].get(spectrum.scan)

            if pl_res is None:
                continue

            seq_split = Pattern.findall(pl_res.mod_seq)
            if len(seq_split) == 1:
                psm = [m for m in spectrum.matches if m.seq == pl_res.seq]

            else:
                pos, i = set(), 0
                for _s in seq_split[:-1]:
                    i += len(_s)
                    pos.add(i)
                # TODO: should identify whether it is a
                # terminal modification
                if pl_res.mod_seq.endswith("]"):
                    pos.add(i + len(seq_split[-1]))
                psm = [m for m in spectrum.matches if m.seq == pl_res.seq
                       and pos.issubset(mod.position
                                        for mod in m.modifications)]

            # TODO: raise error if they are not matched
            if not psm:
                continue

            # merge the two sets of results
            res = self._merge_result(spectrum, psm[0], pl_res)
            res["data_id"] = search_file
            results.append(Result(**res))

        return results

    def _collect_processing_results(self):
        """
        Collect all outputs after processing, including
        the log file, targets, and decoys. If log or decoys
        is missing, warning is raised, and
        [the "target" file is required, otherwise an exception will
        be raised].
        """
        path, file_ = os.path.split(self.result)
        if not path:
            path = os.getcwd()

        files = os.listdir(path)
        pl_files = PercolatorFiles(
            target="percolator.target.psms.txt",
            decoy="percolator.decoy.psms.txt",
            log="percolator.log.txt",
            path=path
        )

        # define warnings and exceptions if some files are missing
        if pl_files.log not in files:
            warnings.warn("Log file is missing, will use file index instead",
                          ResourceWarning)
            pl_files = pl_files._replace(log=None)

        if pl_files.decoy not in files:
            warnings.warn("Decoy processing file is missing", ResourceWarning)
            pl_files = pl_files._replace(decoy=None)

        if pl_files.target not in files:
            raise ValueError("The target processing file is not found.")

        self._percolator_files = pl_files

    def _parse_file_index(self):
        """
        Parse log file to get file index
        """
        if self._percolator_files.log is None:
            return None

        # get the indices of files
        file_index = {}
        with open(os.path.join(self._percolator_files.path,
                               self._percolator_files.log)) as f:
            for line in f:
                if "Assigning index" in line:
                    idx, full_path = re.match(
                        "INFO: Assigning index ([0-9]+) to (.*).",
                        line).groups()
                    res_file_c = os.path.basename(full_path).split(".")
                    file_index[idx] = (
                        ".".join(res_file_c[:-2]) if "pep" in res_file_c
                        else ".".join(res_file_c[:-1])
                    )

                # finish parsing when goes to statistics of matches
                elif line.startswith("INFO: There are"):
                    break

        return file_index if file_index else None

    def _read_results(self, result_file):
        """
        Read the results
        """
        # read results
        res = collections.defaultdict(list)
        with open(result_file) as f:
            f_table = csv.DictReader(f, delimiter="\t")
            for r in f_table:
                idx = r["file_idx"]
                res[idx].append((
                    int(r["scan"]),
                    Percolator(**{field: convert(r[ky])
                                  for ky, field, convert in Pick_up_fields})
                ))

        # convert to dictionary and return
        if self._file_index is not None:
            return {self._file_index.get(idx): dict(items)
                    for idx, items in res.items()}
        return {idx: dict(items) for idx, items in res.items()}

    def _is_compatible(self):
        """
        Check compatibility of the validation results and
        search results
        """
        if self._file_index is None:
            raise ValueError("The file indices are not identified, "
                             "Can't map raw files to Percolator indices.")

        if self._percolator_files.decoy is None:
            warnings.warn("Decoy results are not obtained, only target"
                          "results are returned.")

    @staticmethod
    def _merge_result(spectrum, match, validation):
        """
        Merge the sets of identifications and validations into a
        single search result item.
        """
        merged = {**match.__dict__,
                  **{field: val for field, val
                     in spectrum.__dict__.items() if field != "matches"}}
        vals = []
        for field, val in validation.__dict__.items():
            if field in merged and val is not None:
                merged[field] = val
            else:
                vals.append((field, val))
        merged.update(vals)

        return merged
