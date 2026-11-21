
import struct
import tempfile
import unittest
from pathlib import Path

import numpy as np
from pysam import VariantFile
from pybcf import BcfReader

from .helpers import record, write_bcf

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

class TestSamples(unittest.TestCase):
    ''' class to make sure per-sample FORMAT fields are decoded correctly
    '''

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

    def test_phased_without_samples(self):
        ''' check phasing is an empty bool array for a BCF without samples
        '''
        sites_bcf = Path(__file__).parent / 'data' / 'hapmap_3.3.hg38.shrunk.bcf'
        var = next(BcfReader(sites_bcf))
        self.assertEqual(var.samples.phased.dtype, np.bool_)
        self.assertEqual(len(var.samples.phased), 0)

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
