import dataclasses

from .base import Assignment
from .base import SearchResult as _SearchResult

from .pepxml_reader import PepXMLReader


@dataclasses.dataclass
class Match(Assignment):
    xcorr: float = None
    dcn: float = None
    dcnstar: float = None
    sp: float = None
    evalue: float = None


@dataclasses.dataclass
class SearchResult(_SearchResult):
    mass: float = None
    scan: int = None
    spectrum_native_id: str = None
    retention_time: float = None


class CometReader(PepXMLReader):
    """
    Read Comet search result from pep.xml file
    """
    def __init__(self, search_result):
        super().__init__(search_result)

    def get_result(self):
        """ Read search results """
        # get prefix of the tag in name spaces
        res = self.parse_matches()
        return self._parse_results(res)

    def _parse_results(self, xml_results):
        """ Parse search results from """
        spectrum_matches = []
        for spectrum_query in xml_results:
            query_info = spectrum_query.get("query")
            query_matches = spectrum_query.get("matches")
            # parse the matches
            c = int(query_info.get("assumed_charge"))
            nmatch = None
            matches = []
            for match_info in query_matches:
                match = self._parse_spectrum_match(match_info)
                nmatch = match_info["info"].get("num_matched_peptides")
                if nmatch is not None:
                    nmatch = int(nmatch)
                matches.append(dataclasses.replace(match, charge=c))

            # get search results
            spec_match = SearchResult(
                spectrum_id=query_info["spectrum"],
                mass=float(query_info["precursor_neutral_mass"]),
                nmatch=nmatch,
                matches=matches,
                scan=int(query_info["start_scan"]),
                spectrum_native_id=query_info.get("spectrumNativeID"),
                retention_time=float(query_info["retention_time_sec"]) / 60
            )
            spectrum_matches.append(spec_match)

        xml_results.clear()
        return spectrum_matches

    def _parse_spectrum_match(self, hit):
        """ Construct search result. """
        hit_info = hit.get("info")
        scores = hit.get("scores")
        mods = hit.get("modifications")
        proteins = hit.get("alter_protein")

        # type of peptides
        protein = hit_info.get("protein")
        proteins.append(protein)
        pep_type = "decoy" if protein.startswith(self.df) else "normal"

        return Match(seq=hit_info["peptide"],
                     modifications=mods,
                     delta=float(hit_info["massdiff"]),
                     missed_cleavage=int(hit_info["num_missed_cleavages"]),
                     protein_accessions=proteins,
                     rank=int(hit_info["hit_rank"]),
                     pep_type=pep_type,
                     xcorr=float(scores["xcorr"]),
                     dcn=float(scores["deltacn"]),
                     sp=float(scores["spscore"]),
                     dcnstar=(float(scores["deltacnstar"])
                              if "deltacnstar" in scores else None),
                     evalue=float(scores["expect"]))