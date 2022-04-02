import collections

from .base import Modification

try:
    import xml.etree.cElementTree as et
except ImportError:
    import xml.etree.ElementTree as et


class PepXMLReader():
    """ Read search result pep.xml file. """

    def __init__(self, search_result):
        self._ns_xml = self._get_name_space(search_result)
        self.df = self._get_decoy_prefix(search_result)
        self.result_file = search_result

    def parse_matches(self):
        """ Read search results """
        ns_xml = f"{{{self._ns_xml}}}" if self._ns_xml is not None else ""
        # read search results
        res = []
        read_spectrum = False
        for event, elem in et.iterparse(self.result_file,
                                        events=("start", "end")):
            if event == "start" and elem.tag == f"{ns_xml}spectrum_query":
                query_info = dict(elem.items())
                hit_info = []
                read_spectrum = True

            elif event == "end" and elem.tag == f"{ns_xml}spectrum_query":
                matches = self._parse_query(hit_info)
                res.append({"query": query_info, "matches": matches})
                read_spectrum = False

            elif read_spectrum:
                hit_info.append((event, elem.tag, dict(elem.items())))

            elem.clear()

        return res

    def _parse_query(self, match_list):
        """ Parse query to get spectrum info. """
        # hit sections
        hit_section_index = [i for i, (_, tag, _) in enumerate(match_list)
                             if tag.endswith("search_hit")]
        # get matches
        matches = []
        for i0, i1 in zip(hit_section_index[::2], hit_section_index[1::2]):
            hit_section = match_list[i0: i1+1]
            # check whether it is post-processed by iProphet or PeptideProphet
            post_anal_section_index = [
                i for i, (_, tag, _) in enumerate(hit_section)
                if tag.endswith("analysis_result")
            ]

            # hits with scores
            _, _, match_info = hit_section[0]
            mods, alter_protein, scores = [], [], collections.defaultdict()
            post_anal = None
            # parse post-analysis results
            j0 = i1
            if post_anal_section_index:
                post_anal = self._parse_post_analysis(
                    hit_section, post_anal_section_index
                )
                j0 = post_anal_section_index[0] + 1

            # get the match
            for event, tag, info in hit_section[1:j0]:
                # terminal modifications
                if event == "start":
                    if tag.endswith("modification_info"):
                        mod = self._parse_term_mod(info)
                        if mod is not None:
                            mods.append(mod)
                    # variable and fixed modifications
                    elif tag.endswith("mod_aminoacid_mass"):
                        mods.append(self._parse_mod(info))
                    elif tag.endswith("alternative_protein"):
                        alter_protein.append(self._parse_alter_proteins(info))
                    # search scores
                    elif tag.endswith("search_score"):
                        scores.update(self._parse_score(info))

            match = {"scores": scores, "alter_protein": alter_protein,
                     "modifications": mods, "post_analysis": post_anal,
                     "info": match_info}

            matches.append(match)

        return matches

    def _parse_post_analysis(self, hit, index):
        """ Parse post-analysis results. """
        post_anal = collections.defaultdict()
        for i0, i1 in zip(index[::2], index[1::2]):
            sub_sec = hit[i0: i1 + 1]
            _, _, anal_method = sub_sec[0]
            _, _, prob_info = sub_sec[1]
            # get parameters
            params = collections.defaultdict()
            for event, tag, info in sub_sec[2:]:
                if event=="start" and tag.endswith("parameter"):
                    params.update(self._parse_score(info))

            post_anal[anal_method.get("analysis")] = {
                "prob": float(prob_info.get("probability")),
                "param": params
            }
        return post_anal

    def _parse_term_mod(self, info):
        """ Parse terminal modifications. """
        if "mod_nterm_mass" in info:
            return Modification(mass=float(info.get("mod_nterm_mass")),
                                position="nterm")

        if "mod_cterm_mass" in info:
            return Modification(mass=float(info.get("mod_cterm_mass")),
                                position="cterm")

        return None

    def _parse_mod(self, info):
        """ Parse variable and fixed modifications. """
        return Modification(mass=float(info.get("mass")),
                            position=int(info.get("position")))

    def _parse_score(self, info):
        return {info.get("name"): float(info.get("value"))}

    def _parse_alter_proteins(self, element):
        """ Alternative protein. """
        return element.get("protein")

    def _get_name_space(self, search_result):
        """ Get prefix of the tag in name spaces. """
        for _, elem in et.iterparse(search_result, events=["start-ns"]):
            if elem[0] == "":
                return elem[1] if elem[1] else None

    def _get_decoy_prefix(self, search_result):
        """
        Get prefix for decoy proteins
        """
        param = "parameter"
        if self._ns_xml is not None:
            param = f"{{{self._ns_xml}}}parameter"

        for _, elem in et.iterparse(search_result):
            if elem.tag == param and elem.get("name") == "decoy_prefix":
                decoy_prefix = elem.get("value")
                elem.clear()
                return decoy_prefix
            elem.clear()
