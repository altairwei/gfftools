import warnings
from typing import (
    Dict, List, Tuple, Iterator, Sequence
)
from abc import ABC, abstractmethod
from collections import namedtuple

from HTSeq import GenomicFeature

from pygff.reader import GFF_Reader


class FilterError(Exception):
    pass


class Filter(ABC):
    """A class represents a GFF filter."""

    @abstractmethod
    def validate(self, feature: GenomicFeature) -> bool:
        """Determine if the feature can pass this filter."""


class SimpleValueUnionFilter(Filter):
    """Pass the filter if value meets any one of candidates."""

    def __init__(self, param):
        self.valid_values = []
        # Empty string is not allowed
        if isinstance(param, str) and param:
            self.valid_values.append(param)
        elif isinstance(param, Sequence):
            self.valid_values.extend(param)

    def is_valid(self, value):
        # Do not check if filter is not set.
        if not self.valid_values:
            return True
        if value in self.valid_values:
            return True
        else:
            return False


class SeqIdFilter(SimpleValueUnionFilter):
    def validate(self, feature: GenomicFeature) -> bool:
        return self.is_valid(feature.iv.chrom)


class TypeFilter(SimpleValueUnionFilter):
    def validate(self, feature: GenomicFeature) -> bool:
        return self.is_valid(feature.type)


class SourceFilter(SimpleValueUnionFilter):
    def validate(self, feature: GenomicFeature) -> bool:
        return self.is_valid(feature.source)


class StrandFilter(SimpleValueUnionFilter):
    def validate(self, feature: GenomicFeature) -> bool:
        return self.is_valid(feature.iv.strand)


class AttributesFilter(Filter):
    def __init__(self, param):
        self.attr_pairs = []
        # Empty string is not allowed
        if isinstance(param, str) and param:
            self.attr_pairs.append(param.split("="))
        elif isinstance(param, Sequence):
            for attr_keyval in param:
                self.attr_pairs.append(attr_keyval.split("="))

    def validate(self, feature: GenomicFeature) -> bool:
        # Skip if filter isn't set
        if not self.attr_pairs:
            return True
        matched = []
        for key, val in self.attr_pairs:
            if key in feature.attr and val == feature.attr[key]:
                matched.append(True)
            else:
                matched.append(False)
        if not all(matched):
            return False
        else:
            return True


class ExpressionFilter(Filter):
    def __init__(self, param):
        self.expression = None
        # Emtpry string is not allowed
        if isinstance(param, str) and param:
            self.expression = param

    def validate(self, feature: GenomicFeature) -> bool:
        if not self.expression:
            return True
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

        return bool(eval(self.expression, env))


class RegionsFilter(Filter):
    Region = namedtuple('Region', ['seqname', 'start', 'end'])

    def __init__(self, param):
        self.regions = []
        # Empty string is not allowed
        if isinstance(param, str) and param:
            self.regions.append(self.parse_region(param))
        elif isinstance(param, Sequence):
            for region_string in param:
                self.regions.append(self.parse_region(region_string))

    def parse_region(self, region_string: str):
        seqname, start, end = None, None, None
        parts_1st = region_string.split(":", 1)
        if len(parts_1st) > 0 and parts_1st[0]:
            seqname = parts_1st[0]
        if len(parts_1st) > 1:
            parts_2nd = parts_1st[1].split("-", 1)
            if len(parts_2nd) > 0 and parts_2nd[0]:
                try:
                    # Convert 1-based indexing to 0-based half-open interval
                    start = int(parts_2nd[0]) - 1
                    if start < 0:
                        raise FilterError("Region start must larger than 0")
                except ValueError:
                    warnings.warn(
                        "Start position of region filter is invalid"
                        " and will be ignored", category=RuntimeWarning)

            if len(parts_2nd) > 1 and parts_2nd[1]:
                try:
                    # In the user's view this is closed interval
                    end = int(parts_2nd[1])
                    if end < 1:
                        raise FilterError("Region end must larger than 0")
                except ValueError:
                    warnings.warn(
                        "End position of region filter is invalid"
                        " and will be ignored.", category=RuntimeWarning)

        if start and end and start > end:
            raise FilterError("Region start must be less than end")

        return RegionsFilter.Region(seqname, start, end)

    def validate(self, feature: GenomicFeature) -> bool:
        """
        A feature is retained only if it starts and ends in the region.
        """
        if not self.regions:
            return True

        for region in self.regions:
            if feature.iv.chrom != region.seqname:
                continue
            if region.start and feature.iv.start < region.start:
                continue
            if region.end and feature.iv.end > region.end:
                continue
            return True

        return False


FILTER_NAME_MAP = {
    "seqid": SeqIdFilter,
    "type": TypeFilter,
    "source": SourceFilter,
    "strand": StrandFilter,
    "attributes": AttributesFilter,
    "expression": ExpressionFilter,
    "region": RegionsFilter
}


def make_filters(filter_params: Dict) -> List[Filter]:
    filters = []
    for key, val in filter_params.items():
        if key in FILTER_NAME_MAP.keys():
            filters.append(FILTER_NAME_MAP[key](val))
    return filters


class FilterChain:
    def __init__(self, filter_params):
        self.filters = make_filters(filter_params)

    def add_filter(self, filter: Filter):
        self.filters.append(filter)

    def validate(self, feature):
        for filter in self.filters:
            if not filter.validate(feature):
                return False
        return True


class GFF_Filter(GFF_Reader):

    def __init__(
            self, gff_file: str, filter_params: Dict,
            end_included=True, show_progress=False):
        GFF_Reader.__init__(
            self, gff_file, end_included, show_progress=show_progress)
        self.gff_file = gff_file
        self.filter_chain = FilterChain(filter_params)

    def __iter__(self) -> Iterator[Tuple[GenomicFeature, str]]:
        for feature, raw_line in GFF_Reader.__iter__(self):
            if self.filter_chain.validate(feature):
                yield (feature, raw_line)
            else:
                continue
