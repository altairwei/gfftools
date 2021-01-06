import re
from typing import Dict, List, Tuple, Iterator

import HTSeq
from HTSeq import (
    GenomicFeature,
    GenomicInterval,
    parse_GFF_attribute_string,
    FileOrSequence
)


class GFF_Reader(FileOrSequence):
    """Parse a GFF file (Modified from HTSeq.GFF_Reader)

    Pass the constructor either a file name or an iterator of lines of a
    GFF files. If a file name is specified, it may refer to a gzip compressed
    file.

    Iterating over the object then yields GenomicFeature objects.
    """

    def __init__(self, filename_or_sequence, end_included=True):
        FileOrSequence.__init__(self, filename_or_sequence)
        self.end_included = end_included
        self.metadata = {}

    def __iter__(self) -> Iterator[Tuple[GenomicFeature, str]]:
        for line in FileOrSequence.__iter__(self):
            if isinstance(line, bytes):
                line = line.decode()
            if line == "\n":
                continue
            if line.startswith('#'):
                if line.startswith("##"):
                    mo = re.compile(r"##\s*(\S+)\s+(\S*)").match(line)
                    if mo:
                        self.metadata[mo.group(1)] = mo.group(2)
                continue
            (seqname, source, feature, start, end, score,
             strand, frame, attributeStr) = line.split("\t", 8)
            (attr, name) = parse_GFF_attribute_string(attributeStr, True)
            if self.end_included:
                iv = GenomicInterval(
                        seqname,
                        int(start) - 1, int(end),
                        strand)
            else:
                iv = GenomicInterval(
                        seqname,
                        int(start) - 1, int(end) - 1,
                        strand)
            f = GenomicFeature(name, feature, iv)
            if score != ".":
                score = float(score)
            if frame != ".":
                frame = int(frame)
            f.source = source
            f.score = score
            f.frame = frame
            f.attr = attr
            yield (f, line)