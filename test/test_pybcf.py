
import ast
import gc
import gzip
import os
import resource
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from pysam import VariantFile
from pybcf import BcfReader
import pybcf
import pybcf.reader as pybcf_reader

from .helpers import record, write_bcf

def get_pysam_probs(var, alt_idx):
    ''' get unphased biallelic genotype data
    '''
    geno = np.zeros((len(var.samples), 3), dtype=np.float64)
    
    for i, sample in enumerate(var.samples.itervalues()):
        alleles = sample.allele_indices
        if alleles[0] is None:
            # check for missing data first
            geno[i, :] = float('nan')
        else:
            geno[i, alleles.count(alt_idx)] = 1.0
    
    return geno

def get_pybcf_probs(var, alt_indices):
    ''' get unphased biallelic genotype data
    '''
    geno = var.samples['GT']
    
    # count alts
    is_ref = (geno == 0).all(axis=1)
    is_nan = np.isnan(geno[:, 0])
    hom_ref_idx = np.where(is_ref)
    is_nan_idx = np.where(is_nan)
    
    has_alt = ~is_ref & ~is_nan
    has_alt_idx = np.where(has_alt)[0]
    
    probs = np.zeros((len(geno), 3), dtype=np.float64)
    probs[hom_ref_idx, 0] = 1.0
    probs[is_nan_idx, :] = float('nan')
    has_alt = geno[has_alt_idx]
    
    for alt_idx in alt_indices:
        if alt_idx > 1:
            probs[has_alt_idx] = 0  # clear from previous alts
        n_alts = (has_alt == alt_idx).sum(axis=1)
        probs[has_alt_idx, n_alts] = 1.0
        yield alt_idx, probs

def open_sample_bcf(path):
    ''' start writing a BCF with GT, integer and float FORMAT fields

    It has 37 samples, to cover the SIMD loops (which take 32 bytes or 16 samples
    at a time) as well as the plain loop that decodes whatever is left over.
    '''
    out = VariantFile(path, 'wb')
    out.header.contigs.add('chr1', length=10000)
    out.header.formats.add('GT', 1, 'String', 'Genotype')
    out.header.formats.add('I8', '.', 'Integer', '8-bit signed/missing ints')
    out.header.formats.add('I16', '.', 'Integer', '16-bit signed/missing ints')
    out.header.formats.add('I32', '.', 'Integer', '32-bit signed/missing ints')
    out.header.formats.add('FL', '.', 'Float', 'Float values with missing/EOV')
    for i in range(37):
        out.header.add_sample(f'S{i}')
    return out

# header with every type of INFO and FORMAT field, for records written by hand. Each
# line has an IDX, which is the number the records use for that field.
TYPED_HEADER = (
    '##fileformat=VCFv4.2\n'
    '##FILTER=<ID=PASS,Description="All filters passed",IDX=0>\n'
    '##contig=<ID=chr1,length=1000,IDX=0>\n'
    '##INFO=<ID=FLG,Number=0,Type=Flag,Description="f",IDX=1>\n'
    '##INFO=<ID=IS,Number=1,Type=Integer,Description="is",IDX=2>\n'
    '##INFO=<ID=IV,Number=.,Type=Integer,Description="iv",IDX=3>\n'
    '##INFO=<ID=FS,Number=1,Type=Float,Description="fs",IDX=4>\n'
    '##INFO=<ID=FV,Number=.,Type=Float,Description="fv",IDX=5>\n'
    '##INFO=<ID=SS,Number=1,Type=String,Description="ss",IDX=6>\n'
    '##INFO=<ID=SV,Number=.,Type=String,Description="sv",IDX=7>\n'
    '##INFO=<ID=IV_EOV,Number=.,Type=Integer,Description="iv eov",IDX=8>\n'
    '##INFO=<ID=FV_EOV,Number=.,Type=Float,Description="fv eov",IDX=9>\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp",IDX=10>\n'
    '##FORMAT=<ID=FF,Number=1,Type=Float,Description="ff",IDX=11>\n'
    '##FORMAT=<ID=FT,Number=1,Type=String,Description="ft",IDX=12>\n'
    '##INFO=<ID=FLG_V21,Number=0,Type=Flag,Description="fv21",IDX=13>\n'
    '##FORMAT=<ID=FZ,Number=1,Type=Integer,Description="fz",IDX=14>\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n'
)

class TestBcfReader(unittest.TestCase):
    ''' class to make sure BcfReader works correctly
    '''
    
    def test_without_sample_data(self):
        ''' check this package matches pysam for BCF without sample data
        '''
        path = Path(__file__).parent / 'data' / 'hapmap_3.3.hg38.shrunk.bcf'
        vcf_pysam = VariantFile(path)
        vcf_pybcf = BcfReader(path)
        
        for var_pysam, var_pybcf in zip(vcf_pysam, vcf_pybcf):
            self.assertEqual(var_pysam.chrom, var_pybcf.chrom)
            self.assertEqual(var_pysam.pos, var_pybcf.pos)
            self.assertEqual(var_pysam.ref, var_pybcf.ref)
            self.assertEqual(var_pysam.alts, var_pybcf.alts)
            self.assertEqual(var_pysam.qual, var_pybcf.qual)
            self.assertEqual(list(var_pysam.filter), var_pybcf.filter)
            self.assertEqual(var_pysam.id, var_pybcf.id)
            self.assertEqual(list(var_pybcf.info), list(var_pysam.info))
            
            # check all the info fields match
            for field in var_pysam.info:
                self.assertEqual(var_pysam.info[field], var_pybcf.info[field],)
    
    def test_with_sample_data(self):
        ''' check this package matches pysam for BCF with per sample data
        '''
        path = Path(__file__).parent / 'data' / '1000G.shrunk.bcf'
        vcf_pysam = VariantFile(path)
        vcf_pybcf = BcfReader(path)
        for var_pysam, var_pybcf in zip(vcf_pysam, vcf_pybcf):
            self.assertEqual(var_pysam.chrom, var_pybcf.chrom)
            self.assertEqual(var_pysam.pos, var_pybcf.pos)
            self.assertEqual(var_pysam.ref, var_pybcf.ref)
            self.assertEqual(var_pysam.alts, var_pybcf.alts)
            self.assertEqual(var_pysam.qual, var_pybcf.qual)
            self.assertEqual(list(var_pysam.filter), var_pybcf.filter)
            self.assertEqual(var_pysam.id, var_pybcf.id)
            self.assertEqual(list(var_pybcf.info), list(var_pysam.info))
            
            # check all the info fields match
            for field in var_pysam.info:
                self.assertEqual(var_pysam.info[field], var_pybcf.info[field],)
            
            alt_indices = np.arange(1, len(var_pybcf.alts) + 1)
            for alt_idx, geno_pybcf in get_pybcf_probs(var_pybcf, alt_indices):
                geno_pysam = get_pysam_probs(var_pysam, alt_idx)
                
                is_nan_pysam = np.isnan(geno_pysam).any(axis=1)
                is_nan_pybcf = np.isnan(geno_pybcf).any(axis=1)
                
                self.assertTrue((is_nan_pysam == is_nan_pybcf).all())
                self.assertTrue((geno_pybcf[~is_nan_pysam] == geno_pysam[~is_nan_pysam]).all())
            
            self.assertEqual(list(var_pysam.format), list(var_pybcf.samples))
            for key in var_pybcf.samples:
                if key == 'GT':
                    continue
                data_pybcf = var_pybcf.samples[key]
                data_pysam = [x[key] for x in var_pysam.samples.itervalues()]
                
                if isinstance(data_pybcf, list):
                    for x, y in zip(data_pysam, data_pybcf):
                        if isinstance(x, str) and isinstance(y, tuple):
                            assert len(y) == 1
                            y = y[0]
                        self.assertEqual(x, y)
                
                # pysam and pybcf represent the data differently. pysam returns
                # a ragged list, whereas pybcf inserts nan values in the gaps
                elif len(data_pybcf.shape) > 1:
                    sums_pybcf = (~np.isnan(data_pybcf)).sum(axis=1)
                    sums_pysam = np.array([len(x) for x in data_pysam])
                    
                    # finding values to skip 
                    skip = sums_pybcf == 0
                    if skip.sum() > 0:
                        # the places where sums_bcf is zero correspond to pysam
                        # gave a (None, ) tuple
                        self.assertTrue((sums_pysam[skip] == 1).all())
                        self.assertTrue([x for x, keep in zip(data_pysam, skip) if keep] == [(None, )] * skip.sum())
                        
                        sums_pybcf = sums_pybcf[~skip]
                        sums_pysam = sums_pysam[~skip]
                        
                        data_pybcf = data_pybcf[~skip]
                        data_pysam = [x for x, drop in zip(data_pysam, skip) if not drop]
                    
                    self.assertTrue((sums_pybcf == sums_pysam).all())
                    data_pybcf = [tuple(x[~np.isnan(x)].tolist()) for x in data_pybcf]
                    self.assertEqual(data_pybcf, data_pysam)
                else:
                    data_pysam = np.array(data_pysam)
                    if len(data_pysam.shape) > 1 and max(len(x) for x in data_pysam):
                        data_pysam = data_pysam.T
                    
                    not_nan = ~np.isnan(data_pybcf)
                    data_pybcf = data_pybcf[not_nan]
                    data_pysam = data_pysam[not_nan]
                    
                    for i, (x, y) in enumerate(zip(data_pybcf, data_pysam)):
                        if x != y:
                            print(i, x, y)
                    
                    self.assertTrue((data_pybcf == data_pysam).all())
    
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

    def test_uncompressed_and_plain_gzip_bcf(self):
        ''' check a BCF without compression, or with plain gzip, reads the same as with bgzip
        '''
        path = Path(__file__).parent / 'data' / 'bcfv2.2.bcf'
        with tempfile.TemporaryDirectory() as tmpdir:
            raw_path = Path(tmpdir) / 'raw.bcf'
            raw_path.write_bytes(gzip.open(path).read())
            self.assertEqual(raw_path.read_bytes()[:5], b'BCF\x02\x02')
            gzip_path = Path(tmpdir) / 'gzip.bcf'
            with gzip.open(gzip_path, 'wb') as gz:
                gz.write(raw_path.read_bytes())

            variants = list(BcfReader(path))
            for other_path in (raw_path, gzip_path):
                others = list(BcfReader(other_path))
                self.assertEqual(len(variants), len(others))
                for var, other in zip(variants, others):
                    self.assertEqual((var.chrom, var.pos, var.ref, var.alts),
                                     (other.chrom, other.pos, other.ref, other.alts))
                    self.assertEqual(repr(var.info.items()), repr(other.info.items()))
                    self.assertEqual(list(var.samples), list(other.samples))
                    for key in var.samples:
                        if isinstance(var.samples[key], list):
                            self.assertEqual(var.samples[key], other.samples[key])
                        else:
                            np.testing.assert_array_equal(var.samples[key], other.samples[key])

    @unittest.skipUnless(os.path.exists('/dev/stdin'), 'needs /dev/stdin')
    def test_read_from_pipe(self):
        ''' check a BCF can be read from a pipe, with or without compression
        '''
        path = Path(__file__).parent / 'data' / 'hapmap_3.3.hg38.shrunk.bcf'
        sites = [f'{x.chrom} {x.pos}' for x in BcfReader(path)]
        code = ('from pybcf import BcfReader\n'
                'for x in BcfReader("/dev/stdin"):\n'
                '    print(x.chrom, x.pos)')
        for data in (path.read_bytes(), gzip.open(path).read()):
            piped = subprocess.run([sys.executable, '-c', code], input=data,
                                   capture_output=True, check=True)
            self.assertEqual(piped.stdout.decode().splitlines(), sites)

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

    def test_negative_info_ints_and_empty_info(self):
        ''' check negative int8/int16/int32 INFO values and zero-INFO variants match pysam
        '''
        with tempfile.TemporaryDirectory() as tmpdir:
            bcf_path = Path(tmpdir) / 'neg_info.bcf'
            with VariantFile(bcf_path, 'wb') as out:
                out.header.contigs.add('chr1', length=10000)
                out.header.info.add('I_SCALAR', 1, 'Integer', 'Scalar signed int')
                out.header.info.add('I_VEC', '.', 'Integer', 'Vector signed ints')
                out.header.info.add('F_SCALAR', 1, 'Float', 'Scalar float')

                # Variant 1: int8 negative values (-1, -126)
                r1 = out.new_record(contig='chr1', start=99, alleles=('A', 'G'))
                r1.info['I_SCALAR'] = -5
                r1.info['I_VEC'] = (-126, -1, 0, 127)
                r1.info['F_SCALAR'] = -2.5
                out.write(r1)

                # Variant 2: int16 negative values (-300, -32766)
                r2 = out.new_record(contig='chr1', start=199, alleles=('C', 'T'))
                r2.info['I_SCALAR'] = -300
                r2.info['I_VEC'] = (-32766, -129, 32767)
                out.write(r2)

                # Variant 3: int32 negative values (-100000)
                r3 = out.new_record(contig='chr1', start=299, alleles=('G', 'A'))
                r3.info['I_SCALAR'] = -100000
                r3.info['I_VEC'] = (-2000000, -32769, 100000)
                out.write(r3)

                # Variant 4: zero INFO fields
                r4 = out.new_record(contig='chr1', start=399, alleles=('T', 'C'))
                out.write(r4)

            vcf_pysam = VariantFile(bcf_path)
            vcf_pybcf = BcfReader(bcf_path)
            records = list(zip(vcf_pysam, vcf_pybcf))
            self.assertEqual(len(records), 4)
            for var_pysam, var_pybcf in records:
                self.assertEqual(set(var_pysam.info), set(var_pybcf.info))
                for field in var_pysam.info:
                    self.assertEqual(var_pysam.info[field], var_pybcf.info[field])

    def check_samples_match_pysam(self, path):
        ''' check GT, phasing, and the int and float FORMAT fields of a one-variant BCF match pysam
        '''
        pairs = list(zip(VariantFile(path), BcfReader(path)))
        self.assertEqual(len(pairs), 1)
        for var_pysam, var_pybcf in pairs:
            if 'GT' in var_pysam.format:
                expected_phased = np.array([smp.phased for smp in var_pysam.samples.values()], dtype=bool)
                np.testing.assert_array_equal(var_pybcf.samples.phased, expected_phased)

                gt = var_pybcf.samples['GT']
                for idx, smp in enumerate(var_pysam.samples.values()):
                    alleles = smp.allele_indices
                    for j, a in enumerate(alleles):
                        if a is None:
                            self.assertTrue(np.isnan(gt[idx, j]))
                        else:
                            self.assertEqual(gt[idx, j], float(a))
                    for j in range(len(alleles), gt.shape[1]):
                        self.assertTrue(np.isnan(gt[idx, j]))
            else:
                self.assertEqual(var_pybcf.samples.phased.dtype, np.bool_)
                self.assertEqual(len(var_pybcf.samples.phased), 0)

            for key in ('I8', 'I16', 'I32', 'FL'):
                if key not in var_pysam.format:
                    continue
                arr = var_pybcf.samples[key]
                for idx, smp in enumerate(var_pysam.samples.values()):
                    vals = smp[key]
                    for j, expected in enumerate(vals):
                        if expected is None:
                            self.assertTrue(np.isnan(arr[idx, j]))
                        else:
                            self.assertAlmostEqual(float(arr[idx, j]), float(expected), places=5)
                    for j in range(len(vals), arr.shape[1]):
                        self.assertTrue(np.isnan(arr[idx, j]))

    def test_sample_ints_floats_and_phasing(self):
        ''' check FORMAT parsing of 8-bit GT, phasing, and int8/int16/int32/float vectors
        '''
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'samples.bcf'
            with open_sample_bcf(path) as out:
                rec = out.new_record(contig='chr1', start=100, alleles=('A', 'C', 'G'))
                for i, s in enumerate(out.header.samples):
                    # diploid, missing and haploid genotypes, phased and unphased
                    mod = i % 6
                    if mod == 0:
                        rec.samples[s]['GT'] = (0, 1)
                        rec.samples[s].phased = False
                    elif mod == 1:
                        rec.samples[s]['GT'] = (1, 2)
                        rec.samples[s].phased = True
                    elif mod == 2:
                        rec.samples[s]['GT'] = (None, None)
                        rec.samples[s].phased = False
                    elif mod == 3:
                        rec.samples[s]['GT'] = (None, None)
                        rec.samples[s].phased = True
                    elif mod == 4:
                        rec.samples[s]['GT'] = (1,)
                    else:
                        rec.samples[s]['GT'] = (None,)

                    # missing values, single values, and pairs of values
                    if i % 5 == 0:
                        rec.samples[s]['I8'] = (None,)
                        rec.samples[s]['I16'] = (None,)
                        rec.samples[s]['I32'] = (None,)
                        rec.samples[s]['FL'] = (None,)
                    elif i % 5 == 1:
                        rec.samples[s]['I8'] = (-5,)
                        rec.samples[s]['I16'] = (-300,)
                        rec.samples[s]['I32'] = (1,)  # low 16 bits == 0x0001 regression test
                        rec.samples[s]['FL'] = (-1.25,)
                    else:
                        rec.samples[s]['I8'] = (-120 + i, 50)
                        rec.samples[s]['I16'] = (-2000 + i, 2000 + i)
                        rec.samples[s]['I32'] = (-100000 - i, 65537 + i)
                        rec.samples[s]['FL'] = (float(i) * 0.5, None)
                out.write(rec)
            self.check_samples_match_pysam(path)

    def test_sample_gt_with_many_alleles(self):
        ''' check GT with more than 63 alleles, which BCF has to store as 16-bit values
        '''
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'many_alleles.bcf'
            with open_sample_bcf(path) as out:
                many_alleles = ('A',) + tuple(f'A{"C" * k}' for k in range(1, 70))
                rec = out.new_record(contig='chr1', start=200, alleles=many_alleles)
                for i, s in enumerate(out.header.samples):
                    if i % 3 == 0:
                        rec.samples[s]['GT'] = (0, 68)
                        rec.samples[s].phased = True
                    elif i % 3 == 1:
                        rec.samples[s]['GT'] = (65,)
                    else:
                        rec.samples[s]['GT'] = (None, None)
                        rec.samples[s].phased = False
                out.write(rec)
            self.check_samples_match_pysam(path)

    def test_samples_without_gt(self):
        ''' check a variant with samples but no GT field has empty phasing
        '''
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'no_gt.bcf'
            with open_sample_bcf(path) as out:
                rec = out.new_record(contig='chr1', start=300, alleles=('A', 'T'))
                for i, s in enumerate(out.header.samples):
                    rec.samples[s]['I8'] = (i % 10, 1)
                out.write(rec)
            self.check_samples_match_pysam(path)

    sample_bcf = Path(__file__).parent / 'data' / '1000G.shrunk.bcf'

    def test_wrappers_keep_parents_alive(self):
        ''' check header, info and sample objects work when nothing else holds their reader or variant
        '''
        # the readers and variants here are temporaries, which each result has to keep alive
        contigs = BcfReader(self.sample_bcf).header.contigs
        self.assertGreater(len(contigs), 0)
        info_keys = list(next(BcfReader(self.sample_bcf)).info)
        self.assertGreater(len(info_keys), 0)
        gt = next(BcfReader(self.sample_bcf)).samples['GT']
        self.assertEqual(gt.ndim, 2)

        # variants and the header stay usable after the `with` block closes the reader
        with BcfReader(self.sample_bcf) as bcf:
            var = next(bcf)
        self.assertGreater(len(var.info), 0)
        self.assertEqual(len(var.info), len(list(var.info)))
        self.assertEqual(dict(var.info.items()), {k: var.info[k] for k in var.info})
        self.assertEqual(len(var.samples), len(list(var.samples)))
        sample_items = list(var.samples.items())
        self.assertEqual([k for k, _ in sample_items], list(var.samples))
        np.testing.assert_array_equal(var.samples.get('GT'), gt)
        self.assertEqual(bcf.header.contigs, contigs)

    def test_closed_reader_raises(self):
        ''' check a closed reader raises ValueError when asked for more variants
        '''
        with BcfReader(self.sample_bcf) as bcf:
            chrom = next(bcf).chrom
        with self.assertRaises(ValueError):
            next(bcf)

        # close() can be called more than once
        bcf = BcfReader(self.sample_bcf)
        bcf.close()
        bcf.close()
        with self.assertRaises(ValueError):
            next(bcf)
        with self.assertRaises(ValueError):
            bcf.fetch(chrom)

    def test_wrappers_cannot_be_created_directly(self):
        ''' check the header, variant, info and sample classes raise TypeError if created by hand
        '''
        bcf = BcfReader(self.sample_bcf)
        var = next(bcf)
        for cls in (type(bcf.header), type(var.info), type(var.samples)):
            with self.assertRaises(TypeError):
                cls()
        with self.assertRaises(TypeError):
            type(var)(bcf)
        # `__new__` skips `__init__`, so the wrappers must also refuse to exist without
        # the object they read from
        for cls in (type(var.info), type(var.samples), type(var)):
            for args in ((), (None,)):
                with self.assertRaises(TypeError):
                    cls.__new__(cls, *args)

        # A variant made with `__new__` is empty, and has no contig to name if the
        # header lacks contigs
        with tempfile.TemporaryDirectory() as tmpdir:
            no_contigs = Path(tmpdir) / 'no_contigs.bcf'
            write_bcf(no_contigs, '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            empty = type(var).__new__(type(var), BcfReader(no_contigs))
            self.assertEqual((empty.pos, empty.ref, empty.alts, empty.filter), (0, '', None, []))
            self.assertEqual((len(empty.info), len(empty.samples)), (0, 0))
            with self.assertRaises(IndexError):
                empty.chrom

    def test_phased_without_samples(self):
        ''' check phasing is an empty bool array for a BCF without samples
        '''
        sites_bcf = Path(__file__).parent / 'data' / 'hapmap_3.3.hg38.shrunk.bcf'
        var = next(BcfReader(sites_bcf))
        self.assertEqual(var.samples.phased.dtype, np.bool_)
        self.assertEqual(len(var.samples.phased), 0)

    @unittest.skipUnless(Path('/proc/self/fd').exists(), 'needs /proc to count open files')
    def test_discarded_reader_closes_file(self):
        ''' check a reader with no references left closes its file at once

        This fails if the reader and its header or variants refer to each other in
        a cycle, as the file then stays open until the garbage collector runs.
        '''
        gc.collect()
        gc.disable()
        try:
            before = len(os.listdir('/proc/self/fd'))
            for _ in range(25):
                r = BcfReader(self.sample_bcf)
                r.header.contigs
                r.samples
                var = next(r)
                list(var.info)
                var.samples['GT']
                del r, var
            after = len(os.listdir('/proc/self/fd'))
        finally:
            gc.enable()
        self.assertEqual(before, after)

    def test_type_stubs_match_extension(self):
        ''' check reader.pyi and py.typed exist, parse cleanly, and match Cython classes
        '''
        pkg_dir = Path(pybcf.__file__).parent
        self.assertTrue((pkg_dir / 'py.typed').is_file())
        stub_path = pkg_dir / 'reader.pyi'
        self.assertTrue(stub_path.is_file())

        # The stub is only parsed, never run: type checkers accept `int | None` in stubs
        # on every Python version, but running it fails before Python 3.10.
        source = stub_path.read_text(encoding='utf8')
        tree = ast.parse(source, filename=str(stub_path))

        stub_classes = {
            node.name: {
                item.name for item in node.body
                if isinstance(item, ast.FunctionDef) and not item.name.startswith('_')
            }
            for node in tree.body if isinstance(node, ast.ClassDef)
        }
        runtime_classes = {
            name: {attr for attr in dir(cls) if not attr.startswith('_')}
            for name, cls in vars(pybcf_reader).items()
            if isinstance(cls, type) and not name.startswith('_')
        }
        self.assertEqual(stub_classes, runtime_classes)

    def test_missing_info_ints_and_floats(self):
        ''' check missing scalar/vector INFO ints and floats return None like pysam
        '''
        vcf_txt = (
            '##fileformat=VCFv4.2\n'
            '##contig=<ID=chr1,length=1000>\n'
            '##INFO=<ID=IS,Number=1,Type=Integer,Description="scalar int">\n'
            '##INFO=<ID=IV,Number=.,Type=Integer,Description="vec int">\n'
            '##INFO=<ID=FS,Number=1,Type=Float,Description="scalar float">\n'
            '##INFO=<ID=FV,Number=.,Type=Float,Description="vec float">\n'
            '##INFO=<ID=SS,Number=1,Type=String,Description="scalar str">\n'
            '##INFO=<ID=SV,Number=.,Type=String,Description="vec str">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
            'chr1\t10\t.\tA\tG\t.\t.\tIS=.;IV=1,.,3;FS=.;FV=1.5,.,2.5;SS=hello;SV=a,b,c\n'
            'chr1\t20\t.\tC\tT\t.\t.\tIV=5;FV=3.5;SS=.;SV=only_one\n'
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            vcf_path = Path(tmpdir) / 'missing_info.vcf'
            bcf_path = Path(tmpdir) / 'missing_info.bcf'
            vcf_path.write_text(vcf_txt)
            with VariantFile(vcf_path) as vin, VariantFile(bcf_path, 'wb', header=vin.header) as vout:
                for rec in vin:
                    vout.write(rec)

            pairs = list(zip(VariantFile(bcf_path), BcfReader(bcf_path)))
            self.assertEqual(len(pairs), 2)
            for var_pysam, var_pybcf in pairs:
                for field in var_pysam.info:
                    self.assertEqual(var_pysam.info[field], var_pybcf.info[field])

    def test_zero_length_values_and_end_of_vector(self):
        ''' check values with no elements, and INFO vectors cut short by an end-of-vector value
        '''
        info = (
            b'\x11\x01\x00'  # FLG (flag, n_vals=0) -> True
            b'\x11\x0d\x11\x01'  # FLG_V21 (BCFv2.1 int8=1 flag) -> True
            b'\x11\x02\x01'  # IS (int8, n_vals=0) -> None
            b'\x11\x03\x01'  # IV (int8, n_vals=0) -> ()
            b'\x11\x04\x05'  # FS (float, n_vals=0) -> None
            b'\x11\x05\x05'  # FV (float, n_vals=0) -> ()
            b'\x11\x06\x07'  # SS (char, n_vals=0) -> None
            b'\x11\x07\x07'  # SV (char, n_vals=0) -> ()
            b'\x11\x08\x31\x05\x80\x81'  # IV_EOV = [5, missing, EOV] -> (5, None)
            + b'\x11\x09\x35' + struct.pack('<fII', 1.5, 0x7F800001, 0x7F800002)  # FV_EOV -> (1.5, None)
        )
        indiv = (
            b'\x11\x0a\x01'  # DP (int8, n_vals=0 across 2 samples) -> [nan, nan]
            b'\x11\x0b\x05'  # FF (float, n_vals=0 across 2 samples) -> [nan, nan]
            b'\x11\x0c\x07'  # FT (char, n_vals=0 across 2 samples) -> [('',), ('',)]
            b'\x11\x0e\x00'  # FZ (untyped, n_vals=0 across 2 samples) -> [nan, nan]
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'zero_and_eov.bcf'
            write_bcf(path, TYPED_HEADER,
                      [record(n_info=10, n_fmt=4, n_sample=2, info=info, indiv=indiv)])

            rec_pysam = next(VariantFile(path))
            rec_pybcf = next(BcfReader(path))
            self.assertEqual(list(rec_pysam.info), list(rec_pybcf.info))
            for field in rec_pysam.info:
                self.assertEqual(rec_pysam.info[field], rec_pybcf.info[field])
            self.assertEqual(rec_pybcf.samples['DP'].shape, (2,))
            self.assertTrue(np.isnan(rec_pybcf.samples['DP']).all())
            self.assertEqual(rec_pybcf.samples['FF'].shape, (2,))
            self.assertTrue(np.isnan(rec_pybcf.samples['FF']).all())
            self.assertEqual(rec_pybcf.samples['FT'], [('',), ('',)])
            self.assertEqual(rec_pybcf.samples['FZ'].shape, (2,))
            self.assertTrue(np.isnan(rec_pybcf.samples['FZ']).all())

    def test_nan_and_reserved_values(self):
        ''' check a real NaN is kept, while BCF's missing and reserved values are None or nan
        '''
        info = (
            b'\x11\x02\x11\x87'  # IS = 0x87 (reserved int8) -> None
            + b'\x11\x03\x41\x05\x80\x82\x81'  # IV = [5, missing, reserved 0x82, EOV] -> (5, None, None)
            + b'\x11\x04\x15' + struct.pack('<I', 0x7FC00000)  # FS = quiet NaN -> float('nan')
            + b'\x11\x05\x45' + struct.pack('<fIII', 1.5, 0x7FC00000, 0x7F800001, 0x7F800007)  # FV -> (1.5, nan, None, None)
        )
        indiv = (
            b'\x11\x0a\x11\x82\x14'  # DP = [reserved 0x82, 20] -> [nan, 20.0]
            + b'\x11\x0b\x15' + struct.pack('<If', 0x7F800005, 2.5)  # FF = [reserved sNaN, 2.5] -> [nan, 2.5]
        )
        # QUAL is a quiet NaN, and FILTER holds one missing int8 (0x80) -> []
        rec = record(qual=0x7FC00000, n_info=4, n_fmt=2, n_sample=2, body=b'\x07\x17A\x11\x80',
                     info=info, indiv=indiv)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'nan_and_reserved.bcf'
            write_bcf(path, TYPED_HEADER, [rec])

            var = next(BcfReader(path))
            self.assertTrue(np.isnan(var.qual))
            self.assertEqual(var.filter, [])
            self.assertIsNone(var.info['IS'])
            self.assertEqual(var.info['IV'], (5, None, None))
            self.assertTrue(np.isnan(var.info['FS']))
            fv = var.info['FV']
            self.assertEqual(len(fv), 4)
            self.assertEqual(fv[0], 1.5)
            self.assertTrue(np.isnan(fv[1]))
            self.assertIsNone(fv[2])
            self.assertIsNone(fv[3])
            self.assertTrue(np.isnan(var.samples['DP'][0]))
            self.assertEqual(var.samples['DP'][1], 20.0)
            self.assertTrue(np.isnan(var.samples['FF'][0]))
            self.assertEqual(var.samples['FF'][1], np.float32(2.5))

    def test_reserved_values_in_simd_paths(self):
        ''' check reserved int8/int16/int32/float values are nan with enough samples (36) for SIMD
        '''
        n_smp = 36
        smp_cols = '\t'.join(f'S{i}' for i in range(n_smp))
        header = (
            '##fileformat=VCFv4.2\n'
            '##contig=<ID=chr1,length=1000,IDX=0>\n'
            '##FORMAT=<ID=I8,Number=1,Type=Integer,Description="i8",IDX=1>\n'
            '##FORMAT=<ID=I16,Number=1,Type=Integer,Description="i16",IDX=2>\n'
            '##FORMAT=<ID=I32,Number=1,Type=Integer,Description="i32",IDX=3>\n'
            '##FORMAT=<ID=FL,Number=1,Type=Float,Description="fl",IDX=4>\n'
            f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{smp_cols}\n'
        )
        # even samples have a reserved value, odd samples have the lowest ordinary value
        indiv = (
            b'\x11\x01\x11' + b''.join(b'\x87' if i % 2 == 0 else b'\x88' for i in range(n_smp))
            + b'\x11\x02\x12' + b''.join(struct.pack('<h', -32761 if i % 2 == 0 else -32760) for i in range(n_smp))
            + b'\x11\x03\x13' + b''.join(struct.pack('<i', -2147483641 if i % 2 == 0 else -2147483640) for i in range(n_smp))
            + b'\x11\x04\x15' + b''.join(struct.pack('<I', 0x7F800007 if i % 2 == 0 else 0x7FC00000) for i in range(n_smp))
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'simd_reserved.bcf'
            write_bcf(path, header, [record(n_fmt=4, n_sample=n_smp, indiv=indiv)])

            var = next(BcfReader(path))
            is_even = [i % 2 == 0 for i in range(n_smp)]
            np.testing.assert_array_equal(np.isnan(var.samples['I8']), is_even)
            np.testing.assert_array_equal(var.samples['I8'][1::2], np.full(n_smp // 2, -120.0))
            np.testing.assert_array_equal(np.isnan(var.samples['I16']), is_even)
            np.testing.assert_array_equal(var.samples['I16'][1::2], np.full(n_smp // 2, -32760.0))
            np.testing.assert_array_equal(np.isnan(var.samples['I32']), is_even)
            np.testing.assert_array_equal(var.samples['I32'][1::2], np.full(n_smp // 2, -2147483640.0))
            self.assertTrue(np.isnan(var.samples['FL']).all())

    def test_corrupt_or_missing_bcf_raises_value_error(self):
        ''' check GzReader raises ValueError on missing or truncated/corrupt BCF streams
        '''
        with tempfile.TemporaryDirectory() as tmpdir:
            missing_path = Path(tmpdir) / 'does_not_exist.bcf'
            with self.assertRaises(ValueError):
                BcfReader(missing_path)

            src_bcf = Path(__file__).parent / 'data' / '1000G.shrunk.bcf'
            raw = src_bcf.read_bytes()
            truncated_path = Path(tmpdir) / 'truncated.bcf'
            # Truncate mid-file so header succeeds and variant/gzip reading detects corruption
            truncated_path.write_bytes(raw[:len(raw) // 2])
            with self.assertRaises(ValueError):
                list(BcfReader(truncated_path))

            # A header length near 4 GB with only a few bytes of text must fail at EOF,
            # not by allocating the declared length upfront
            huge_header = Path(tmpdir) / 'huge_header.bcf'
            with gzip.open(huge_header, 'wb') as gz:
                gz.write(b'BCF\x02\x02' + struct.pack('<I', 0xFFFFFFFF) + b'##fileformat=VCFv4.2\n')
            peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            with self.assertRaisesRegex(ValueError, 'truncated BCF header text'):
                BcfReader(huge_header)
            grown_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss - peak_kb
            self.assertLess(grown_kb, 100_000)

    def test_corrupt_bgzf_blocks_raise(self):
        ''' check a BGZF block with a bad checksum or length is rejected, as is data after the last block
        '''
        raw = (Path(__file__).parent / 'data' / 'bcfv2.2.bcf').read_bytes()
        # the first block's header holds its size minus 1, and the block ends with a
        # checksum, then the length of its uncompressed data
        size = struct.unpack_from('<H', raw, 16)[0] + 1
        length = struct.unpack_from('<I', raw, size - 4)[0]
        bad_checksum = bytearray(raw)
        bad_checksum[size - 8] ^= 0xff
        bad_length = bytearray(raw)
        struct.pack_into('<I', bad_length, size - 4, length + 1)

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'corrupt.bcf'
            for contents, error in ((bad_checksum, 'incorrect data check'),
                                    (bad_length, 'incorrect length check'),
                                    (raw + b'extra', 'invalid BGZF block header')):
                path.write_bytes(contents)
                with self.assertRaisesRegex(ValueError, error):
                    list(BcfReader(path))

    def test_reading_stops_after_corrupt_record(self):
        ''' check a reader stops after a failed read, rather than carrying on from mid-record
        '''
        header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'bad_contig.bcf'
            # The first record names a contig missing from the header. The bytes after
            # its fixed fields would pass for a record, as would the record after it.
            write_bcf(path, header, [record(chrom=7, body=record()), record(pos0=19)])
            bcf = BcfReader(path)
            with self.assertRaisesRegex(ValueError, 'unknown contig index'):
                next(bcf)
            self.assertEqual(list(bcf), [])

    def test_key_lookup(self):
        ''' check INFO/FORMAT keys that are absent give False, the default, or a KeyError
        '''
        sample_bcf = Path(__file__).parent / 'data' / '1000G.shrunk.bcf'
        var = next(BcfReader(sample_bcf))
        self.assertFalse('NO_SUCH_INFO' in var.info)
        self.assertIsNone(var.info.get('NO_SUCH_INFO'))
        self.assertEqual(var.info.get('NO_SUCH_INFO', 'fallback'), 'fallback')
        first_info_key = next(iter(var.info))
        self.assertEqual(var.info.get(first_info_key), var.info[first_info_key])
        with self.assertRaises(KeyError):
            _ = var.info['NO_SUCH_INFO']
        self.assertFalse('NO_SUCH_FMT' in var.samples)
        self.assertIsNone(var.samples.get('NO_SUCH_FMT'))
        self.assertEqual(var.samples.get('NO_SUCH_FMT', 'fallback'), 'fallback')
        with self.assertRaises(KeyError):
            _ = var.samples['NO_SUCH_FMT']

    def test_corrupt_info_raises(self):
        ''' check a record with a damaged INFO or shared block raises ValueError
        '''
        header = (
            '##fileformat=VCFv4.2\n'
            '##FILTER=<ID=PASS,Description="All filters passed">\n'
            '##contig=<ID=chr1,length=1000>\n'
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'bad_info.bcf'

            # INFO is only parsed when used, so a field numbered 99, which the header
            # lacks, raises on lookup
            write_bcf(path, header, [record(n_info=1, info=b'\x11\x63\x00')])
            var = next(BcfReader(path))
            with self.assertRaises(ValueError):
                _ = 'DP' in var.info
            with self.assertRaises(ValueError):
                _ = var.info['DP']

            # a REF that claims 5 bytes when only 1 is left must raise ValueError, rather
            # than end the iteration
            write_bcf(path, header, [record(body=b'\x07\x57A')])
            with self.assertRaises(ValueError):
                list(BcfReader(path))

    def test_corrupt_records_raise(self):
        ''' check records with impossible positions, IDs, genotypes or sample counts raise ValueError
        '''
        header = (
            '##fileformat=VCFv4.2\n'
            '##contig=<ID=chr1,length=1000>\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'bad_rec.bcf'
            for bad_record in (
                # pos0 == INT32_MAX (would overflow signed int32 on pos0 + 1)
                record(pos0=0x7FFFFFFF),
                # ID encoded as int8 instead of a string
                record(body=b'\x11\x05\x17A\x00'),
                # ID with an overflow length descriptor holding 0 integers instead of exactly 1
                record(body=b'\xf7\x01\x00\x17A\x00'),
                # n_sample (5) does not match the header's sample count (1)
                record(n_fmt=1, n_sample=5, indiv=b'\x11\x01\x01'),
            ):
                write_bcf(path, header, [bad_record])
                with self.assertRaises(ValueError):
                    next(BcfReader(path))

            # GT encoded as a string instead of integers. Sample data is only parsed
            # when used, so the record reads, then raises when GT is looked up
            write_bcf(path, header, [record(n_fmt=1, n_sample=1, indiv=b'\x11\x01\x17A')])
            var = next(BcfReader(path))
            with self.assertRaises(ValueError):
                _ = var.samples['GT']

    def test_repeated_keys_first_wins(self):
        ''' check that if a record repeats an INFO or FORMAT key, the first copy wins, like htslib
        '''
        header = (
            '##fileformat=VCFv4.2\n'
            '##contig=<ID=chr1,length=1000>\n'
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        )
        rec = record(n_info=2, n_fmt=2, n_sample=1,
                     info=b'\x11\x01\x11\x05\x11\x01\x11\x07',   # DP=5, then DP=7
                     indiv=b'\x11\x02\x11\x0a\x11\x02\x11\x0b')  # GQ=10, then GQ=11
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'dup_keys.bcf'
            write_bcf(path, header, [rec])

            var_pysam = next(VariantFile(path))
            var_pybcf = next(BcfReader(path))
            self.assertEqual(var_pybcf.info['DP'], var_pysam.info['DP'])
            self.assertEqual(var_pybcf.info['DP'], 5)
            self.assertEqual(var_pybcf.samples['GQ'][0], var_pysam.samples['S1']['GQ'])
            self.assertEqual(var_pybcf.samples['GQ'][0], 10)

class TestBcfv2_2(unittest.TestCase):
    ''' class to work through a BCFv2.2 file with some edge cases in data fields
    '''
    def test_missing_values(self):
        ''' check we can work with BCFs version 2.2
        '''
        path = Path(__file__).parent / 'data' / 'bcfv2.2.bcf'
        vcf_pysam = VariantFile(path)
        vcf_pybcf = BcfReader(path)
        
        for var_pysam, var_pybcf in zip(vcf_pysam, vcf_pybcf):
            # check all the info fields match
            self.assertEqual(list(var_pybcf.info), list(var_pysam.info))
            for field in var_pysam.info:
                self.assertEqual(var_pysam.info[field], var_pybcf.info[field])
            
            self.assertEqual(list(var_pysam.format), list(var_pybcf.samples))
            for key in var_pybcf.samples:
                data_pybcf = var_pybcf.samples[key]
                data_pysam = [x[key] for x in var_pysam.samples.itervalues()]
                
                if isinstance(data_pybcf, list):
                    for x, y in zip(data_pysam, data_pybcf):
                        if isinstance(x, str) and isinstance(y, tuple):
                            assert len(y) == 1
                            y = y[0]
                        self.assertEqual(x, y)
                
                # pysam and pybcf represent the data differently. pysam returns
                # a ragged list, whereas pybcf inserts nan values in the gaps
                elif len(data_pybcf.shape) > 1:
                    lens_pybcf = (~np.isnan(data_pybcf)).sum(axis=1)
                    lens_pysam = np.array([len(x) for x in data_pysam])
                    self.assertTrue((lens_pybcf == lens_pysam).all())
                    data_pybcf = [tuple(x[~np.isnan(x)].tolist()) for x in data_pybcf]
                    self.assertEqual(data_pybcf, data_pysam)
                else:
                    data_pysam = np.array(data_pysam)
                    if len(data_pysam.shape) > 1 and max(len(x) for x in data_pysam):
                        data_pysam = data_pysam.T
                    self.assertTrue((data_pybcf == data_pysam).all())
                
