#!/usr/bin/env python

import argparse
import sys
import json
import re
from argparse import ArgumentError, Namespace
from typing import Dict, List, Tuple, Iterator
from collections import defaultdict

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


def attr_to_string(attrs: Dict):
    attr_list = []
    for id_attr in ("gene_id", "transcript_id"):
        if id_attr in attrs:
            attr_list.append('{} "{}"'.format(id_attr, attrs[id_attr]))
            del attrs[id_attr]
    attr_list.extend(
        '{} "{}"'.format(str(key), str(val)) for (key, val) in attrs.items()
    )
    return "; ".join(attr_list)


def split_prefix(x: str, delimiter: str) -> Tuple[str, str]:
    try:
        i = x.index(delimiter)
        return (x[0:i], x[i+1:])
    except ValueError:
        return ("", x)


def get_gtf_line(
    feature: GenomicFeature,
    transcript_parent: Dict[str, str],
    id_prefix: List[str],
    type_mapping: Dict[str, str],
    type_delimiter: str
):
    attr_dict: Dict[str, str] = {}
    attr_dict.update(feature.attr)

    if feature.type in ("mRNA", "tRNA", "rRNA"):
        if "transcript_id" not in attr_dict:
            attr_dict["transcript_id"] = split_prefix(attr_dict["ID"], type_delimiter)[1]

    # Fill necessary attributes with information extracted from `Parent`
    if "Parent" in attr_dict:
        parent_type, parent_id = split_prefix(attr_dict["Parent"], type_delimiter)
        # Remove ID prefix if necessary
        if parent_type in id_prefix:
            full_id = attr_dict["Parent"]
        else:
            full_id = parent_id
        # Fill in the required attributes for GTF format
        if parent_type in ("transcript", "rna"):
            attr_dict["transcript_id"] = full_id
        elif parent_type == "gene":
            attr_dict["gene_id"] = full_id
            # Record gene-transcript relation
            if "transcript_id" in attr_dict:
                transcript_parent[attr_dict["transcript_id"]] = full_id

    if "gene_id" not in attr_dict:
        if feature.type == "gene":
            attr_dict["gene_id"] = split_prefix(attr_dict["ID"], type_delimiter)[1]
        else:
            if "transcript_id" in attr_dict:
                attr_dict["gene_id"] = transcript_parent[attr_dict["transcript_id"]]
   

    # Reserve or replace ID prefix
    for prefix in id_prefix:
        attr_key = prefix + "_id"
        if attr_key in attr_dict:
            if not attr_dict[attr_key].startswith(prefix):
                attr_dict[attr_key] = prefix + type_delimiter + attr_dict[attr_key]

    # Change feature type
    if feature.type in type_mapping:
        feature.type = type_mapping[feature.type]

    if feature.type == "exon":
        # Check for gene_id and transcript_id exists
        if "gene_id" not in attr_dict or "transcript_id" not in attr_dict:
            raise Exception("Exon must contain both 'gene_id' and 'transcript_id'")

    return (
        "\t".join(
            [
                feature.iv.chrom,
                feature.source,
                feature.type,
                str(
                    feature.iv.start + 1
                ),  # See https://htseq.readthedocs.io/en/master/genomic.html#HTSeq.GenomicInterval
                str(
                    feature.iv.end
                ),  # See https://htseq.readthedocs.io/en/master/genomic.html#HTSeq.GenomicInterval
                feature.score,
                feature.iv.strand,
                str(feature.frame),
                attr_to_string(attr_dict),
            ]
        )
        + "\n"
    )


def stats_action(options: Namespace) -> None:
    #TODO: Guess end included or not. See: https://htseq.readthedocs.io/en/master/features.html 
    summary = {
        "closed_intervals": defaultdict(int),
        "half_open_intervals": defaultdict(int),
        "types": defaultdict(int)
    }
    for _, line in GFF_Reader(options.gff_file):
        (seqname, source, feature_type, start, end, score,
            strand, frame, attributeStr) = line.split("\t", 8)
        if feature_type in ("exon", "CDS", "start_codon", "stop_codon"):
            # The GFF specification is unclear on whether the end coordinate marks
            # the last base-pair of the feature (closed intervals, end_included=True)
            # or the one after (half-open intervals, end_included=False).
            is_end_included = ((int(end) - int(start) + 1) % 3) == 0
            if is_end_included:
                summary["closed_intervals"][feature_type] += 1
            else:
                summary["half_open_intervals"][feature_type] += 1
        summary["types"][feature_type] += 1
    json.dump(summary, sys.stdout, indent=2)
    print("", file=sys.stdout)


def convert_action(options: Namespace) -> None:
    type_mapping = {}
    for type_aes in options.type_mapping:
        old_type, new_type = type_aes.split(":")
        type_mapping[old_type] = new_type
    gff3 = GFF_Reader(options.gff_file)
    i = 0
    transcript_parent = {}
    for feature, _ in gff3:
        line = get_gtf_line(
            feature, transcript_parent, options.id_prefix,
            type_mapping, options.type_delimiter
        )
        sys.stdout.write(line)
        i += 1
        if i % 100000 == 0:
            print("%d GFF lines processed." % i, file=sys.stderr)


def filter_action(options: Namespace) -> None:
    for feature, raw_line in GFF_Reader(options.gff_file):
        if options.seqid and feature.iv.chrom not in options.seqid:
            continue
        if options.type and feature.type not in options.type:
            continue
        if options.source and feature.source not in options.source:
            continue
        if options.attributes:
            matched = []
            for attr_keyval in options.attributes:
                key, val = attr_keyval.split("=")
                if key in feature.attr and val == feature.attr[key]:
                    matched.append(True)
                else:
                    matched.append(False)
            if not all(matched):
                continue

        # Print out selected fields
        if options.print_field == "all":
            sys.stdout.write(raw_line)
        elif options.print_field == "attributes":
            (*_, attributeStr) = raw_line.split("\t", 8)
            sys.stdout.write(attributeStr)
        else:
            if options.print_field in feature.attr:
                sys.stdout.write(
                    feature.attr[options.print_field] + "\n")


def seq_action(options: Namespace) -> None:
    pass


if __name__ == "__main__":
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("gff_file", help="GFF3 file obtained from Ensembl.", type=str, metavar="GFF_FILE")

    parser = argparse.ArgumentParser(description="GFF tool.")
    subparsers = parser.add_subparsers()

    stats_cmd = subparsers.add_parser(
        "stats", help="Print overview stats of GFF file.", parents=[parent_parser]
    )
    stats_cmd.set_defaults(func=stats_action)

    convert_cmd = subparsers.add_parser(
        "conv", help="Converts Ensembl's favored GFF3 to GTF.", parents=[parent_parser]
    )
    convert_cmd.set_defaults(func=convert_action)
    convert_cmd.add_argument(
        "-p",
        "--reserve-id-prefix",
        dest="id_prefix",
        action="append",
        default=[],
        help="The ID prefix to be reserved. Example `-p transcript -p gene`"
        " will produce 'transcript:G9200.1' and 'gene:G9200' as IDs.",
    )
    convert_cmd.add_argument(
        "-t",
        "--transform-feature-type",
        dest="type_mapping",
        action="append",
        default=[],
        help="Rewrite the type of a feature to another one. Example `-t mRNA:transcript`"
        " will change the type of `mRNA` feature to `transcript` type.",
    )
    convert_cmd.add_argument(
        "-d",
        "--type-delimiter",
        dest="type_delimiter",
        default=":",
        help="Value of `Parent` or `ID` attribute generally contains an ID prefixed with a type, "
        "%(prog)s will use this prefix to determine the type of `Parent`. Therefore you "
        "need to specify the delimiter used to separate the prefix from the real ID. (default: %(default)s)",
    )
    convert_cmd.add_argument(
        "-e",
        "--end-included",
        dest="end_included",
        action="store_true",
        default=True,
        help="Specifies whether the end coordinate of marks the last base-pair"
        " in output GFF file. (default: %(default)s)",
    )

    filter_cmd = subparsers.add_parser(
        "filter", help="Filter records in GFF files based on specified parameters.", parents=[parent_parser]
    )
    filter_cmd.set_defaults(func=filter_action)
    filter_cmd.add_argument(
        "-i",
        "--seqid",
        dest="seqid",
        action="append",
        default=[],
        help="Filter records with given seqid (aka. chromosome name), such as `-i chr1A -i chr3B`.",
    )
    filter_cmd.add_argument(
        "-s",
        "--source",
        dest="source",
        action="append",
        default=[],
        help="Filter records with given source, such as `-s IWGSC -s Genbank`.",
    )
    filter_cmd.add_argument(
        "-t",
        "--type",
        dest="type",
        action="append",
        default=[],
        help="Filter records with given feature types, such as `-t exon -t CDS`.",
    )
    filter_cmd.add_argument(
        "-a",
        "--attributes",
        dest="attributes",
        action="append",
        default=[],
        help="Filter records with given key value pairs, such as `-a ID=GENE0545 -a Name=nad2`."
        "Note that this option behaves differently from the other filtering options in that "
        "a record will only pass the filter if all of the specified attributes match.",
    )
    filter_cmd.add_argument(
        "-p",
        "--print-field",
        dest="print_field",
        default="all",
        help="Specify which field to print after filtering. Available fields: "
        "seqid, source, type, start, end, score, strand, phase and attributes. "
        "Any other value will be treated as a key of attributes and the value "
        "of that key will be printed out. (default: %(default)s)",
    )

    seq_cmd = subparsers.add_parser(
        "seq", help="Extract sequences from FASTA files based on GFF annotation.", parents=[parent_parser]
    )
    seq_cmd.set_defaults(func=seq_action)
    seq_cmd.add_argument(
        "-g",
        "--genome-file",
        dest="genome",
        help="Full path to a multi-fasta file with the genomic sequences.",
    )

    options = parser.parse_args()
    try:
        options.func(options)
    except BrokenPipeError:
        pass

