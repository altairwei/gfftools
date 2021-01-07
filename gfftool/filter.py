import sys
from typing import Dict, Tuple, Iterator, MutableSequence

from HTSeq import GenomicFeature

from gfftool.reader import GFF_Reader


class GFF_Filter(GFF_Reader):

    def __init__(self, gff_file: str, filter_params: Dict, end_included=True, show_progress=False):
        GFF_Reader.__init__(self, gff_file, end_included, show_progress=show_progress)
        self.gff_file = gff_file
        self.check_seqid = False
        self.check_type = False
        self.check_source = False
        self.check_strand = False
        self.check_attributes = False
        self.eval_expression = False
        self.params = {}
        if ("seqid" in filter_params and filter_params["seqid"] and
            isinstance(filter_params["seqid"], MutableSequence)):
            self.params["seqid"] = filter_params["seqid"]
            self.check_seqid = True
        if ("type" in filter_params and filter_params["type"] and
            isinstance(filter_params["type"], MutableSequence)):
            self.params["type"] = filter_params["type"]
            self.check_type = True
        if ("source" in filter_params and filter_params["source"] and
            isinstance(filter_params["source"], MutableSequence)):
            self.params["source"] = filter_params["source"]
            self.check_source = True
        if ("strand" in filter_params and filter_params["strand"] and
            isinstance(filter_params["strand"], MutableSequence)):
            self.params["strand"] = filter_params["strand"]
            self.check_strand = True
        if ("attributes" in filter_params and filter_params["attributes"] and
            isinstance(filter_params["attributes"], MutableSequence)):
            self.params["attributes"] = []
            for attr_keyval in filter_params["attributes"]:
                self.params["attributes"].append(attr_keyval.split("="))
            self.check_attributes = True
        if ("expression" in filter_params and filter_params["expression"] and
            isinstance(filter_params["expression"], str)):
            self.params["expression"] = filter_params["expression"]
            self.eval_expression = True

    def __iter__(self) -> Iterator[Tuple[GenomicFeature, str]]:
        for feature, raw_line in GFF_Reader.__iter__(self):
            if (self.check_seqid and
                feature.iv.chrom not in self.params["seqid"]):
                continue
            if (self.check_type and
                feature.type not in self.params["type"]):
                continue
            if (self.check_source and
                feature.source not in self.params["source"]):
                continue
            if (self.check_strand and
                feature.iv.strand not in self.params["strand"]):
                continue
            if self.check_attributes:
                matched = []
                for key, val in self.params["attributes"]:
                    if key in feature.attr and val == feature.attr[key]:
                        matched.append(True)
                    else:
                        matched.append(False)
                if not all(matched):
                    continue
            if self.eval_expression:
                env = {
                    "seqid": feature.iv.chrom,
                    "source": feature.source,
                    "type": feature.type,
                    "start": feature.iv.start + 1,
                    "end": feature.iv.end,
                    "score": feature.score,
                    "strand": feature.iv.strand,
                    "phase": str(feature.frame),
                    "attributes": feature.attr
                }

                cond = bool(eval(self.params["expression"], env))

                if not cond:
                    continue

            yield (feature, raw_line)
