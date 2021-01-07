import re
import os
from typing import Dict, List, Tuple, Iterator

import HTSeq
from HTSeq import (
    GenomicFeature,
    GenomicInterval,
    parse_GFF_attribute_string
)

from gfftool.utils import ProgressBar

class TextFile(object):
    def __init__(self, filename, show_progress=False):
        self.show_progress = show_progress
        self.filename = filename
        self.filesize = os.stat(self.filename).st_size
        self.line_no = None

    def __iter__(self):
        self.line_no = 1
        lines = open(self.filename, encoding="UTF-8")
        try:
            if self.show_progress:
                with ProgressBar(self.filesize, "Processing: ", "bytes") as bar:
                    for line in lines:
                        line_size = len(line.encode("UTF-8"))
                        bar.update(line_size)
                        yield line
                        self.line_no += 1
            else:
                for line in lines:
                    yield line
                    self.line_no += 1
        finally:
            if isinstance(self.filename, str):
                lines.close()
        self.line_no = None

    def __repr__(self):
        if isinstance(self.filename, str):
            return "<%s object, connected to file name '%s'>" % (
                self.__class__.__name__, self.filename)
        else:
            return "<%s object, connected to %s >" % (
                self.__class__.__name__, repr(self.filename))

    def get_line_number_string(self):
        if self.line_no is None:
            if isinstance(self.filename, str):
                return "file %s closed" % self.filename
            else:
                return "file closed"
        if isinstance(self.filename, str):
            return "line %d of file %s" % (self.line_no, self.filename)
        else:
            return "line %d" % self.line_no


class GFF_Reader(TextFile):
    """Parse a GFF file (Modified from HTSeq.GFF_Reader)

    Pass the constructor either a file name or an iterator of lines of a
    GFF files. If a file name is specified, it may refer to a gzip compressed
    file.

    Iterating over the object then yields GenomicFeature objects.
    """

    def __init__(self, filename_or_sequence, end_included=True, show_progress=False):
        TextFile.__init__(self, filename_or_sequence, show_progress)
        self.end_included = end_included
        self.metadata = {}

    def __iter__(self) -> Iterator[Tuple[GenomicFeature, str]]:
        for line in TextFile.__iter__(self):
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