import dataclasses

from .base import Assignment
from .base import SearchResult as _SearchResult

from .pepxml_reader import PepXMLReader


@dataclasses.dataclass
class ProphetMatch(Assignment):
    # Sequest/Comet scores
    xcorr: float = None
    dcnstar: float = None
    dcn: float = None
    sp: float = None
    # X!Tandem scores
    hyperscore: float = None
    nextscore: float = None
    yscore: float = None
    bscore: float = None
    # OMSSA scores
    pvalue: float = None
    # Expect scores
    evalue: float = None
    # PeptideProphet probability
    prob: float = None
    params: dict = None
    # iProphet probability
    iprob: float = None
    iparams: dict = None


@dataclasses.dataclass
class SearchResult(_SearchResult):
    mass: float = None
    scan: int = None
    spectrum_native_id: str = None
    retention_time: float = None


class ProphetReader(PepXMLReader):
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
                retention_time=(float(query_info["retention_time_sec"]) / 60
                                if "retention_time_sec" in query_info
                                else None)
            )
            spectrum_matches.append(spec_match)

        xml_results.clear()
        return spectrum_matches

    def _parse_spectrum_match(self, hit):
        """ Construct search result. """
        hit_info = hit.get("info")
        mods = hit.get("modifications")
        proteins = hit.get("alter_protein")
        post_anal = hit.get("post_analysis")

        # scores of post-analysis by PeptideProphet or iProphet
        if post_anal is None:
            raise ValueError("The post-analysis using PeptideProphet"
                             "/InterProphet is not performed.")
        post_anal_scores = self._parse_post_analysis_probs(post_anal)

        # search engine scores
        engine_scores = self._parse_scores(hit.get("scores"))

        # proteins
        proteins.append(hit_info.get("protein"))

        return ProphetMatch(
            seq=hit_info["peptide"],
            modifications=mods,
            delta=float(hit_info["massdiff"]),
            missed_cleavage=int(hit_info["num_missed_cleavages"]),
            protein_accessions=proteins,
            rank=int(hit_info["hit_rank"]),
            pep_type="normal",
            **post_anal_scores, **engine_scores
        )

    def _parse_scores(self, scores_info):
        """ Scores based on search engines. """
        def _comet_scores(scores):
            scores_ = {"xcorr": float(scores["xcorr"]),
                       "dcn": float(scores["deltacn"]),
                       "sp": float(scores["spscore"])}
            if "deltacnstar" in scores:
                scores_["dcnstar"] = float(scores["deltacnstar"])
            return scores_

        def _tandem_scores(scores):
            scores_ = {"hyperscore": float(scores["hyperscore"]),
                       "nextscore": float(scores["nextscore"])}
            if "yscore" in scores:
                scores_["yscore"] = float(scores["yscore"])
                scores_["bscore"] = float(scores["bscore"])
            return scores_

        def _omssa_scores(scores):
            scores_ = {"pvalue": float(scores["pvalue"])}
            return scores_

        scores = None
        for engine_score in [_comet_scores, _tandem_scores, _omssa_scores]:
            try:
                scores = engine_score(scores_info)
            except KeyError:
                pass

        if scores is None:
            raise ValueError("Unrecognizable engine scores, Comet/Sequest, "
                             "X! Tandem and OMSSA are currengly implemented.")

        if "expect" in scores_info:
            scores["evalue"] = float(scores_info["expect"])
        return scores

    def _parse_post_analysis_probs(self, scores_info) -> dict:
        """ Parse scores from post analysis. """
        scores = []
        if "peptideprophet" in scores_info:
            prophet_score = scores_info["peptideprophet"]
            scores.append(("prob", prophet_score["prob"]))
            scores.append(("params", prophet_score["param"]))

        if "interprophet" in scores_info:
            prophet_score = scores_info["interprophet"]
            scores.append(("iprob", prophet_score["prob"]))
            scores.append(("iparams", prophet_score["param"]))

        return dict(scores)
