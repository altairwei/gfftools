from typing import Dict

from pyfaidx import Fasta


def genome_extract(
    fasta_file: str, chromosome: str = None,
    start: int = None, end: int = None, strand = "+"
) -> str:
    genome = Fasta(fasta_file)
    if not chromosome:
        raise Exception("Chromosome name must be provided.")
    if not start or not end:
        raise Exception("Position start and end must be provided.")
    # Note: pyfaidx uses 0-based indexing
    seq_obj = genome[chromosome][slice(start, end)]
    if strand == "+":
        return seq_obj.seq
    if strand == "-":
        return seq_obj.reverse.complement
