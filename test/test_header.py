
import tempfile
import unittest
from pathlib import Path

from pysam import VariantFile
from pybcf import BcfReader

from .helpers import record, write_bcf

class TestHeader(unittest.TestCase):
    ''' class to make sure the header of a BCF is parsed correctly
    '''

    def test_header_access(self):
        ''' check this package matches pysam for the header fields
        '''
        path = Path(__file__).parent / 'data' / 'hapmap_3.3.hg38.shrunk.bcf'
        vcf_pysam = VariantFile(path)
        vcf_pybcf = BcfReader(path)
        
        self.assertEqual(list(vcf_pysam.header.contigs), vcf_pybcf.header.contigs)
        self.assertEqual(list(vcf_pysam.header.info), vcf_pybcf.header.info)
        self.assertEqual(list(vcf_pysam.header.filters), vcf_pybcf.header.filters)
        self.assertEqual(list(vcf_pysam.header.formats), vcf_pybcf.header.formats)
        self.assertEqual(list(vcf_pysam.header.samples), vcf_pybcf.header.samples)
        
        path = Path(__file__).parent / 'data' / '1000G.shrunk.bcf'
        vcf_pysam = VariantFile(path)
        vcf_pybcf = BcfReader(path)
        
        self.assertEqual(list(vcf_pysam.header.contigs), vcf_pybcf.header.contigs)
        self.assertEqual(list(vcf_pysam.header.info), vcf_pybcf.header.info)
        self.assertEqual(list(vcf_pysam.header.filters), vcf_pybcf.header.filters)
        self.assertEqual(list(vcf_pysam.header.formats), vcf_pybcf.header.formats)
        self.assertEqual(list(vcf_pysam.header.samples), vcf_pybcf.header.samples)

    def test_header_quoted_commas_and_blank_lines(self):
        ''' check header parsing handles blank lines, and quoted commas/equals/backslashes
        '''
        header = (
            '##fileformat=VCFv4.2\n'
            '\n'
            '#CHROMOSOME_NOTE=ignored single-hash comment\n'
            '##FILTER=<ID=LowQual,Description="Low quality",IDX=1>\n'
            '##FILTER=<ID=PASS,Description="All filters passed",IDX=5>\n'
            '##contig=<ID=chr1,length=1000>\n'
            '##INFO=<ID=REAL_ID,Number=1,Type=Integer,Description="has comma, ID=FAKE_ID, and slash \\\\",IDX=2>  \t\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'hdr.bcf'
            write_bcf(path, header)

            bcf = BcfReader(path)
            self.assertEqual(bcf.header.contigs, ['chr1'])
            self.assertEqual(bcf.header.info, ['REAL_ID'])
            self.assertEqual(bcf.header.filters, ['LowQual', 'PASS'])
            self.assertEqual(list(bcf), [])

    def test_header_without_idx_shares_one_dictionary(self):
        ''' check a header without IDX numbers FILTER/INFO/FORMAT names together, ignoring repeats
        '''
        header = (
            '##fileformat=VCFv4.2\n'
            '##contig=<ID=chr1,length=1000>\n'
            '##contig=<ID=chr1,length=1000>\n'
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="info dp">\n'
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="dup info dp">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="fmt dp">\n'
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="gq">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        )
        # DP is number 1 as both an INFO and a FORMAT field, GT is 2 and GQ is 3
        rec = record(n_info=1, n_fmt=3, n_sample=1, info=b'\x11\x01\x11\x63',
                     indiv=b'\x11\x02\x21\x02\x04\x11\x01\x11\x0f\x11\x03\x11\x2a')
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'shared_dict.bcf'
            write_bcf(path, header, [rec])

            bcf = BcfReader(path)
            self.assertEqual(bcf.header.contigs, ['chr1'])
            self.assertEqual(bcf.header.formats, ['DP', 'GT', 'GQ'])
            var = next(bcf)
            self.assertEqual(var.info['DP'], 99)
            self.assertEqual(float(var.samples['DP'][0]), 15.0)
            self.assertEqual(float(var.samples['GQ'][0]), 42.0)

    def test_invalid_headers_raise(self):
        ''' check headers with two #CHROM lines, or a bad ID or IDX, are rejected
        '''
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'bad_hdr.bcf'
            for lines in (
                '#CHROM\tPOS\n',
                '##contig=<ID=chr1,IDX=-1>\n',
                '##contig=<ID=chr1,IDX=99999999999999999999>\n',
                '##INFO=<Number=1,Type=Integer>\n',
                # IDX has to be on every FILTER/INFO/FORMAT line, or none
                '##INFO=<ID=A,Number=1,Type=Integer>\n##INFO=<ID=B,Number=1,Type=Integer,IDX=1>\n',
                '##INFO=<ID=A,Number=1,Type=Integer,IDX=1>\n##INFO=<ID=B,Number=1,Type=Integer>\n',
            ):
                write_bcf(path, '##fileformat=VCFv4.2\n' + lines + '#CHROM\tPOS\n')
                with self.assertRaises(ValueError):
                    BcfReader(path)

    def test_header_idx_limits(self):
        ''' check IDX values too large to be real are rejected, with a message saying so
        '''
        lines = '##fileformat=VCFv4.2\n{}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'idx.bcf'
            for line in (
                '##INFO=<ID=A,Number=1,Type=Integer,IDX=65537>',
                '##FORMAT=<ID=A,Number=1,Type=Integer,IDX=65537>',
                '##FILTER=<ID=A,IDX=65537>',
                '##contig=<ID=chr1,IDX=4194305>',
                '##contig=<ID=chr1,IDX=99999999999999999999>',
            ):
                write_bcf(path, lines.format(line))
                with self.assertRaisesRegex(ValueError, 'IDX field out of range'):
                    BcfReader(path)

            # contigs can be far more numerous than INFO/FORMAT/FILTER fields
            write_bcf(path, lines.format('##contig=<ID=chr1,IDX=70000>'))
            self.assertEqual(BcfReader(path).header.contigs, ['chr1'])
            write_bcf(path, lines.format('##INFO=<ID=A,Number=1,Type=Integer,IDX=65536>'))
            self.assertEqual(BcfReader(path).header.info, ['A'])
