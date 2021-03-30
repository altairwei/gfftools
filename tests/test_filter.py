import unittest
import tempfile
import os
from contextlib import contextmanager

from pygff.filter import GFF_Filter, FilterError


@contextmanager
def tempinput(data):
    temp = tempfile.NamedTemporaryFile(delete=False)
    temp.write(data.encode("UTF-8"))
    temp.close()
    try:
        yield temp.name
    finally:
        os.unlink(temp.name)


def filter_gff(gff_file, filter_params):
    output_lines = []
    for _, raw_line in GFF_Filter(gff_file, filter_params):
        output_lines.append(raw_line)
    return "".join(output_lines)


GTF_CONTENT = """140	Twinscan	inter	5141	8522	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	inter_CNS	8523	9711	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	inter	9712	13182	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	3UTR	65149	65487	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	3UTR	66823	66992	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	stop_codon	66993	66995	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	66996	66999	.	-	1	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	intron_CNS	70103	70151	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	70207	70294	.	-	2	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	71696	71807	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	start_codon	71805	71806	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	start_codon	73222	73222	.	-	2	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	73222	73222	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	5UTR	73223	73504	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
381	Twinscan	exon	150	200	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	300	401	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	380	401	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	501	650	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	501	650	.	+	2	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	700	800	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	700	707	.	+	2	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	900	1000	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	start_codon	380	382	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	stop_codon	708	710	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
"""


class GffFilterTestCase(unittest.TestCase):

    def test_seqid_filter(self):
        with tempinput(GTF_CONTENT) as gff_file:
            # Empty filter should be ignored
            for param in [[], ""]:
                with self.subTest(seqid=param):
                    self.assertEqual(filter_gff(gff_file, {"seqid": param}), GTF_CONTENT)

            for param in [["381"], "381"]:
                with self.subTest(seqid=param):
                    self.assertEqual(filter_gff(gff_file, {"seqid": param}),
"""381	Twinscan	exon	150	200	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	300	401	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	380	401	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	501	650	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	501	650	.	+	2	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	700	800	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	700	707	.	+	2	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	900	1000	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	start_codon	380	382	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	stop_codon	708	710	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
"""
                    )

    def test_source_filter(self):
        gff_content = """381	Twinscan	CDS	380	401	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	havana	exon	700	800	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	insdc	start_codon	380	382	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	insdc	stop_codon	708	710	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
"""
        with tempinput(gff_content) as gff_file:
            for param in [[], ""]:
                with self.subTest(seqid=param):
                    self.assertEqual(filter_gff(gff_file, {"source": param}), gff_content)
            for param in [["insdc"], "insdc"]:
                with self.subTest(seqid=param):
                    self.assertEqual(filter_gff(gff_file, {"source": param}),
"""381	insdc	start_codon	380	382	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	insdc	stop_codon	708	710	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
"""
                    )
            self.assertEqual(filter_gff(gff_file, {"source": ["Twinscan", "havana"]}),
"""381	Twinscan	CDS	380	401	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	havana	exon	700	800	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
"""
            )

    def test_type_filter(self):
        with tempinput(GTF_CONTENT) as gff_file:
            for param in [[], ""]:
                with self.subTest(seqid=param):
                    self.assertEqual(filter_gff(gff_file, {"type": param}), GTF_CONTENT)
            for param in [["CDS"], "CDS"]:
                with self.subTest(seqid=param):
                    self.assertEqual(filter_gff(gff_file, {"type": param}),
"""140	Twinscan	CDS	66996	66999	.	-	1	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	70207	70294	.	-	2	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	71696	71807	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	73222	73222	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
381	Twinscan	CDS	380	401	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	501	650	.	+	2	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	700	707	.	+	2	gene_id "381.000"; transcript_id "381.000.1";
"""
                    )

            self.assertEqual(filter_gff(gff_file, {"type": ["CDS", "stop_codon"]}),
"""140	Twinscan	stop_codon	66993	66995	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	66996	66999	.	-	1	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	70207	70294	.	-	2	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	71696	71807	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	73222	73222	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
381	Twinscan	CDS	380	401	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	501	650	.	+	2	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	700	707	.	+	2	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	stop_codon	708	710	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
"""
            )

            self.assertEqual(filter_gff(gff_file, {"type": ["start_codon", "stop_codon"]}),
"""140	Twinscan	stop_codon	66993	66995	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	start_codon	71805	71806	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	start_codon	73222	73222	.	-	2	gene_id "140.000"; transcript_id "140.000.1";
381	Twinscan	start_codon	380	382	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	stop_codon	708	710	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
"""
            )

    def test_strand_filter(self):
        with tempinput(GTF_CONTENT) as gff_file:
            for param in [[], ""]:
                with self.subTest(seqid=param):
                    self.assertEqual(filter_gff(gff_file, {"strand": param}), GTF_CONTENT)
            self.assertEqual(filter_gff(gff_file, {"strand": ["+", "-"]}), GTF_CONTENT)
            for param in [["+"], "+"]:
                with self.subTest(seqid=param):
                    self.assertEqual(filter_gff(gff_file, {"strand": param}),
"""381	Twinscan	exon	150	200	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	300	401	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	380	401	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	501	650	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	501	650	.	+	2	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	700	800	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	700	707	.	+	2	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	900	1000	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	start_codon	380	382	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	stop_codon	708	710	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
"""
                    )
            for param in [["-"], "-"]:
                with self.subTest(seqid=param):
                    self.assertEqual(filter_gff(gff_file, {"strand": param}),
"""140	Twinscan	inter	5141	8522	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	inter_CNS	8523	9711	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	inter	9712	13182	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	3UTR	65149	65487	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	3UTR	66823	66992	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	stop_codon	66993	66995	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	66996	66999	.	-	1	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	intron_CNS	70103	70151	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	70207	70294	.	-	2	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	71696	71807	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	start_codon	71805	71806	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	start_codon	73222	73222	.	-	2	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	73222	73222	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	5UTR	73223	73504	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
"""
                    )

    def test_attributes_filter(self):
        gtf_content = """140	Twinscan	CDS	66996	66999	.	-	1	gene_id "140.000"; transcript_id "140.000.1"; ensembl_phase "2"; rank "11"
140	Twinscan	CDS	73222	73222	.	-	0	gene_id "140.000"; transcript_id "140.000.1"; ensembl_phase "0"; rank "2"
140	Twinscan	5UTR	73223	73504	.	-	.	gene_id "140.000"; transcript_id "140.000.1"
"""

        gff3_content = """140	Twinscan	CDS	66996	66999	.	-	1	gene_id=140.000;transcript_id=140.000.1;ensembl_phase=2;rank=11
140	Twinscan	CDS	73222	73222	.	-	0	gene_id=140.000;transcript_id=140.000.1;ensembl_phase=0;rank=2
140	Twinscan	5UTR	73223	73504	.	-	.	gene_id=140.000;transcript_id=140.000.1
"""

        with tempinput(gtf_content) as gtf_file:
            with tempinput(gff3_content) as gff3_file:
                for param in [[], ""]:
                    with self.subTest(seqid=param):
                        self.assertEqual(filter_gff(gtf_file, {"attributes": param}), gtf_content)
                        self.assertEqual(filter_gff(gff3_file, {"attributes": param}), gff3_content)
                for param in [["gene_id=140.000"], "gene_id=140.000"]:
                    with self.subTest(seqid=param):
                        self.assertEqual(filter_gff(gtf_file, {"attributes": param}), gtf_content)
                        self.assertEqual(filter_gff(gff3_file, {"attributes": param}), gff3_content)

                self.assertEqual(filter_gff(gtf_file, {"attributes": ["gene_id=140.000", "ensembl_phase=2", "rank=11"]}), 
"""140	Twinscan	CDS	66996	66999	.	-	1	gene_id "140.000"; transcript_id "140.000.1"; ensembl_phase "2"; rank "11"
"""
                    )
                self.assertEqual(filter_gff(gff3_file, {"attributes": ["gene_id=140.000", "ensembl_phase=2", "rank=11"]}), 
"""140	Twinscan	CDS	66996	66999	.	-	1	gene_id=140.000;transcript_id=140.000.1;ensembl_phase=2;rank=11
"""
                    )

    def test_expression_filter(self):
        with tempinput(GTF_CONTENT) as gff_file:
            self.assertEqual(filter_gff(gff_file, {"expression": []}), GTF_CONTENT)
            self.assertEqual(filter_gff(gff_file, {
                "expression": "seqid == '140' and type == 'inter' and start == 5141"
            }), """140	Twinscan	inter	5141	8522	.	-	.	gene_id ""; transcript_id "";\n""")
            self.assertEqual(filter_gff(gff_file, {
                "expression": "(end - start) == 3381"
            }), """140	Twinscan	inter	5141	8522	.	-	.	gene_id ""; transcript_id "";\n""")

    def test_region_filter(self):
        with tempinput(GTF_CONTENT) as gff_file:
            for param in [[], ""]:
                with self.subTest(seqid=param):
                    self.assertEqual(filter_gff(gff_file, {"region": param}), GTF_CONTENT)
            for param in [["140:8523-9711"], "140:8523-9711"]:
                with self.subTest(seqid=param):
                    self.assertEqual(filter_gff(gff_file, {"region": param}),
"""140	Twinscan	inter_CNS	8523	9711	.	-	.	gene_id ""; transcript_id "";
"""
                    )
            # Test multiple regions
            self.assertEqual(filter_gff(gff_file, {"region": ["140:8523-9711", "140:65149-65487"]}),
"""140	Twinscan	inter_CNS	8523	9711	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	3UTR	65149	65487	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
"""
            )
            # Only seqname
            self.assertEqual(filter_gff(gff_file, {"region": ["140"]}),
"""140	Twinscan	inter	5141	8522	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	inter_CNS	8523	9711	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	inter	9712	13182	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	3UTR	65149	65487	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	3UTR	66823	66992	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	stop_codon	66993	66995	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	66996	66999	.	-	1	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	intron_CNS	70103	70151	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	70207	70294	.	-	2	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	71696	71807	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	start_codon	71805	71806	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	start_codon	73222	73222	.	-	2	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	73222	73222	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	5UTR	73223	73504	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
"""
            )
            # Only start position
            self.assertEqual(filter_gff(gff_file, {"region": ["140:71805", "381:700"]}),
"""140	Twinscan	start_codon	71805	71806	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	start_codon	73222	73222	.	-	2	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	CDS	73222	73222	.	-	0	gene_id "140.000"; transcript_id "140.000.1";
140	Twinscan	5UTR	73223	73504	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
381	Twinscan	exon	700	800	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	700	707	.	+	2	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	900	1000	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	stop_codon	708	710	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
"""
            )
            # Only end position
            self.assertEqual(filter_gff(gff_file, {"region": ["140:-13182", "381:-401"]}),
"""140	Twinscan	inter	5141	8522	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	inter_CNS	8523	9711	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	inter	9712	13182	.	-	.	gene_id ""; transcript_id "";
381	Twinscan	exon	150	200	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	exon	300	401	.	+	.	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	CDS	380	401	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
381	Twinscan	start_codon	380	382	.	+	0	gene_id "381.000"; transcript_id "381.000.1";
"""
            )

        # Test exactly region match
        # TODO: test `end_included` option
        region_exactly_test = """140	Twinscan	inter_CNS	8523	9711	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	inter_CNS	8522	9711	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	inter_CNS	8523	9712	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	3UTR	65149	65487	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
"""
        with tempinput(region_exactly_test) as gff_file:
            self.assertEqual(
                filter_gff(gff_file, {"region": ["140:8523-9711", "140:65149-65487"]}),
"""140	Twinscan	inter_CNS	8523	9711	.	-	.	gene_id ""; transcript_id "";
140	Twinscan	3UTR	65149	65487	.	-	.	gene_id "140.000"; transcript_id "140.000.1";
"""
            )

        # Test wrong params
        with tempinput(GTF_CONTENT) as gff_file:
            with self.assertRaisesRegex(FilterError, "Region start must larger than 0"):
                filter_gff(gff_file, {"region": "140:0-490"})
            with self.assertRaisesRegex(FilterError, "Region end must larger than 0"):
                filter_gff(gff_file, {"region": "140:1-0"})
            with self.assertRaisesRegex(FilterError, "Region start must be less than end"):
                filter_gff(gff_file, {"region": "140:170-12"})
            with self.assertWarnsRegex(RuntimeWarning, "Start position of region filter is invalid"):
                filter_gff(gff_file, {"region": "140:a-132"})
            with self.assertWarnsRegex(RuntimeWarning, "End position of region filter is invalid"):
                filter_gff(gff_file, {"region": "140:1-b"})