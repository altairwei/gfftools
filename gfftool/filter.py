from typing import Dict, Tuple, Iterator

from HTSeq import GenomicFeature

from gfftool.reader import GFF_Reader


class GFF_Filter(GFF_Reader):

    def __init__(self, gff_file: str, filter_params: Dict, end_included=True):
        GFF_Reader.__init__(self, gff_file, end_included)
        self.gff_file = gff_file
        self.params = filter_params

    def __iter__(self) -> Iterator[Tuple[GenomicFeature, str]]:
        for feature, raw_line in GFF_Reader.__iter__(self):
            if ("seqid" in self.params and self.params["seqid"] and
                feature.iv.chrom not in self.params["seqid"]):
                continue
            if ("type" in self.params and self.params["type"] and
                feature.type not in self.params["type"]):
                continue
            if ("source" in self.params and self.params["source"] and
                feature.source not in self.params["source"]):
                continue
            if "attributes" in self.params:
                matched = []
                for attr_keyval in self.params["attributes"]:
                    key, val = attr_keyval.split("=")
                    if key in feature.attr and val == feature.attr[key]:
                        matched.append(True)
                    else:
                        matched.append(False)
                if not all(matched):
                    continue

            yield (feature, raw_line)
