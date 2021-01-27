# PyGFF - Useful tools for working with GTF/GFF3 format files.

## How to install gfftools

Install PyGFF from git repo via pip:

```shell
pip install git+https://github.com/altairwei/PyGFF.git
```

## How to use`gfftools

`gfftools` contains a series of subcommands, includes `stats`, `conv`, `filter` and `seq`.

### Filter GFF features with given conditions

#### Match specified fields

We can specify multiple fields of the same type, and as long as the feature meets one of them, it will be preserved.

For example:

```shell
gfftools filter --seqid 1 --seqid 3 --source ensembl --type CDS Homo_sapiens.GRCh38.99.gtf > chr1_or_chr3_CDS_from_ensembl.gtf
```

#### All attributes need to be matched

We can specify multiple attributes at the same time, and only features that satisfy all of them will be preserved.

For example:

```shell
gfftools filter --attributes gene_version=5 --attributes gene_biotype=lncRNA > lncRNA.gtf
```

#### Evaluate python condition expression

The filter will execute a user-specified python conditional expression for each feature, and the feature will be preserved if the expression is true. The environment in which the expression is executed contains 9 predefined variables.

```python
seqid: str # chromosome name
source: str # Where the feature is generated
type: str # Feature type, CDS, exon and so on
start: int # Start position, 1-based indexing
end: int # End position, 1-based indexing
score: str # The score of the feature
strand: str # strand of DNA
phase: str # aka. frame
attributes: dict # Feature attributes are wraped in a dict.
```

For example:

```shell
gfftools filter -e "(end - start + 1) % 3 == 0" Homo_sapiens.GRCh38.99.gtf > triad.gtf
```

### Extract gene sequences from the genome based on GFF file

In general, we want to extract gene sequences from genome with a small number of features, so the `gfftools seq` command supports the same GFF filter as `gfftools filter`.

```shell
gfftools seq -a gene_id=ENSG00000177133 -g Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.99.gtf > ENSG00000177133.fa
```

The file `ENSG00000177133.fa` in above example would be a multi-fasta file with the sequences corresponding to features of filtered GFF content. The sequence names of genome fasta file (in this case, it means `Homo_sapiens.GRCh38.dna.primary_assembly.fa`) should match what in 1st column of the input GFF file.

### Convert GFF3 to GTF

Work in progress.


### Count features' properties of GFF file

```shell
gfftools stats Homo_sapiens.GRCh38.99.gtf
```