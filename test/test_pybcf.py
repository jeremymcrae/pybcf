
import gc
import gzip
import os
import struct
import tempfile
import unittest
from pathlib import Path

import numpy as np
from pysam import VariantFile
from pybcf import BcfReader

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

    def test_header_quoted_commas_and_blank_lines(self):
        ''' check header parsing handles blank lines, quoted commas/equals/backslashes, and IDX errors
        '''
        header_text = (
            '##fileformat=VCFv4.2\n'
            '\n'
            '##FILTER=<ID=PASS,Description="All filters passed",IDX=0>\n'
            '##contig=<ID=chr1,length=1000>\n'
            '##INFO=<ID=REAL_ID,Number=1,Type=Integer,Description="has comma, ID=FAKE_ID, and slash \\\\",IDX=1>\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
            '\x00'
        ).encode('utf8')

        def write_header_bcf(path, hdr_bytes):
            with gzip.open(path, 'wb') as gz:
                gz.write(b'BCF\x02\x02')
                gz.write(struct.pack('<I', len(hdr_bytes)))
                gz.write(hdr_bytes)

        with tempfile.TemporaryDirectory() as tmpdir:
            bcf_path = Path(tmpdir) / 'hdr.bcf'
            write_header_bcf(bcf_path, header_text)

            bcf = BcfReader(bcf_path)
            self.assertEqual(bcf.header.contigs, ['chr1'])
            self.assertEqual(bcf.header.info, ['REAL_ID'])
            self.assertEqual(bcf.header.filters, ['PASS'])
            self.assertEqual(list(bcf), [])

            for bad_hdr in (
                '##fileformat=VCFv4.2\n#CHROM\tPOS\n#CHROM\tPOS\n\x00',
                '##fileformat=VCFv4.2\n##contig=<ID=chr1,IDX=-1>\n#CHROM\tPOS\n\x00',
                '##fileformat=VCFv4.2\n##contig=<ID=chr1,IDX=99999999999999999999>\n#CHROM\tPOS\n\x00',
            ):
                bad_path = Path(tmpdir) / 'bad_hdr.bcf'
                write_header_bcf(bad_path, bad_hdr.encode('utf8'))
                with self.assertRaises(ValueError):
                    BcfReader(bad_path)

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

    def test_sample_simd_int_widths_floats_and_phasing(self):
        ''' check SIMD + scalar FORMAT parsing for int8/int16/int32, floats, phasing, and 16-bit GT
        '''
        n_samples = 37  # exercises 32-byte / 16-sample AVX2 loops plus scalar tails
        sample_names = [f'S{i}' for i in range(n_samples)]

        with tempfile.TemporaryDirectory() as tmpdir:
            bcf_path = Path(tmpdir) / 'simd_samples.bcf'
            with VariantFile(bcf_path, 'wb') as out:
                out.header.contigs.add('chr1', length=10000)
                out.header.formats.add('GT', 1, 'String', 'Genotype')
                out.header.formats.add('I8', '.', 'Integer', '8-bit signed/missing ints')
                out.header.formats.add('I16', '.', 'Integer', '16-bit signed/missing ints')
                out.header.formats.add('I32', '.', 'Integer', '32-bit signed/missing ints')
                out.header.formats.add('FL', '.', 'Float', 'Float values with missing/EOV')
                for s in sample_names:
                    out.header.add_sample(s)

                # Variant 1: diploid/haploid/missing 8-bit GT + int8/int16/int32/float vectors
                r1 = out.new_record(contig='chr1', start=100, alleles=('A', 'C', 'G'))
                for i, s in enumerate(sample_names):
                    mod = i % 6
                    if mod == 0:
                        r1.samples[s]['GT'] = (0, 1)
                        r1.samples[s].phased = False
                    elif mod == 1:
                        r1.samples[s]['GT'] = (1, 2)
                        r1.samples[s].phased = True
                    elif mod == 2:
                        r1.samples[s]['GT'] = (None, None)
                        r1.samples[s].phased = False
                    elif mod == 3:
                        r1.samples[s]['GT'] = (None, None)
                        r1.samples[s].phased = True
                    elif mod == 4:
                        r1.samples[s]['GT'] = (1,)
                    else:
                        r1.samples[s]['GT'] = (None,)

                    if i % 5 == 0:
                        r1.samples[s]['I8'] = (None,)
                        r1.samples[s]['I16'] = (None,)
                        r1.samples[s]['I32'] = (None,)
                        r1.samples[s]['FL'] = (None,)
                    elif i % 5 == 1:
                        r1.samples[s]['I8'] = (-5,)
                        r1.samples[s]['I16'] = (-300,)
                        r1.samples[s]['I32'] = (1,)  # low 16 bits == 0x0001 regression test
                        r1.samples[s]['FL'] = (-1.25,)
                    else:
                        r1.samples[s]['I8'] = (-120 + i, 50)
                        r1.samples[s]['I16'] = (-2000 + i, 2000 + i)
                        r1.samples[s]['I32'] = (-100000 - i, 65537 + i)
                        r1.samples[s]['FL'] = (float(i) * 0.5, None)
                out.write(r1)

                # Variant 2: >63 alleles forces 16-bit GT encoding in BCF
                many_alleles = ('A',) + tuple(f'A{"C" * k}' for k in range(1, 70))
                r2 = out.new_record(contig='chr1', start=200, alleles=many_alleles)
                for i, s in enumerate(sample_names):
                    if i % 3 == 0:
                        r2.samples[s]['GT'] = (0, 68)
                        r2.samples[s].phased = True
                    elif i % 3 == 1:
                        r2.samples[s]['GT'] = (65,)
                    else:
                        r2.samples[s]['GT'] = (None, None)
                        r2.samples[s].phased = False
                out.write(r2)

                # Variant 3: samples present without GT field in FORMAT
                r3 = out.new_record(contig='chr1', start=300, alleles=('A', 'T'))
                for i, s in enumerate(sample_names):
                    r3.samples[s]['I8'] = (i % 10, 1)
                out.write(r3)

            vcf_pysam = VariantFile(bcf_path)
            vcf_pybcf = BcfReader(bcf_path)
            for var_pysam, var_pybcf in zip(vcf_pysam, vcf_pybcf):
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

    def test_lifetimes_context_manager_and_sites_only_samples(self):
        ''' check Cython wrapper lifetimes, context-manager close, FD cleanup, and 0-sample phased
        '''
        sample_bcf = Path(__file__).parent / 'data' / '1000G.shrunk.bcf'
        sites_bcf = Path(__file__).parent / 'data' / 'hapmap_3.3.hg38.shrunk.bcf'

        # Chained temporaries must keep parent C++ objects alive
        contigs = BcfReader(sample_bcf).header.contigs
        self.assertGreater(len(contigs), 0)
        info_keys = list(next(BcfReader(sample_bcf)).info)
        self.assertGreater(len(info_keys), 0)
        gt = next(BcfReader(sample_bcf)).samples['GT']
        self.assertEqual(gt.ndim, 2)

        # Variants and header remain usable after exiting `with`, while further reads fail cleanly
        with BcfReader(sample_bcf) as bcf:
            var = next(bcf)
        self.assertGreater(len(var.info), 0)
        self.assertEqual(len(var.info), len(list(var.info)))
        self.assertEqual(dict(var.info.items()), {k: var.info[k] for k in var.info})
        self.assertEqual(len(var.samples), len(list(var.samples)))
        sample_items = list(var.samples.items())
        self.assertEqual([k for k, _ in sample_items], list(var.samples))
        np.testing.assert_array_equal(var.samples.get('GT'), gt)
        self.assertEqual(bcf.header.contigs, contigs)
        with self.assertRaises(ValueError):
            next(bcf)

        # Explicit close() is idempotent and prevents further iteration or fetch
        explicit_bcf = BcfReader(sample_bcf)
        explicit_bcf.close()
        explicit_bcf.close()
        with self.assertRaises(ValueError):
            next(explicit_bcf)
        with self.assertRaises(ValueError):
            explicit_bcf.fetch(contigs[0])

        # Sites-only BCF has 0 samples; phased should return an empty bool array
        sites_var = next(BcfReader(sites_bcf))
        self.assertEqual(sites_var.samples.phased.dtype, np.bool_)
        self.assertEqual(len(sites_var.samples.phased), 0)

        # Creating and discarding BcfReader instances must not leak file descriptors
        if Path('/proc/self/fd').exists():
            gc.collect()
            before = len(os.listdir('/proc/self/fd'))
            for _ in range(25):
                r = BcfReader(sample_bcf)
                list(next(r).info)
                del r
            gc.collect()
            after = len(os.listdir('/proc/self/fd'))
            self.assertEqual(before, after)

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

            for var_pysam, var_pybcf in zip(VariantFile(bcf_path), BcfReader(bcf_path)):
                for field in var_pysam.info:
                    self.assertEqual(var_pysam.info[field], var_pybcf.info[field])

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

    def test_key_lookup_and_corrupt_record_errors(self):
        ''' check missing INFO/FORMAT keys return False/KeyError while corrupt records raise ValueError
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

        header_text = (
            '##fileformat=VCFv4.2\n'
            '##FILTER=<ID=PASS,Description="All filters passed">\n'
            '##contig=<ID=chr1,length=1000>\n'
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
            '\x00'
        ).encode('utf8')

        with tempfile.TemporaryDirectory() as tmpdir:
            # Record with corrupt INFO ID index (99) must raise ValueError on `in` and `[]`
            bad_info_bcf = Path(tmpdir) / 'bad_info.bcf'
            shared_info = (
                struct.pack('<iiiIII', 0, 9, 1, 0x7F800001, (1 << 16) | 1, 0)
                + b'\x07'      # empty varid
                + b'\x17A'     # ref = "A"
                + b'\x00'      # empty filters
                + b'\x11\x63\x00'  # INFO ID index 99 (unknown), flag value
            )
            with gzip.open(bad_info_bcf, 'wb') as gz:
                gz.write(b'BCF\x02\x02' + struct.pack('<I', len(header_text)) + header_text)
                gz.write(struct.pack('<II', len(shared_info), 0) + shared_info)

            bad_var = next(BcfReader(bad_info_bcf))
            with self.assertRaises(ValueError):
                _ = 'DP' in bad_var.info
            with self.assertRaises(ValueError):
                _ = bad_var.info['DP']

            # Record with overflowing string length in shared block must raise ValueError, not StopIteration
            bad_shared_bcf = Path(tmpdir) / 'bad_shared.bcf'
            shared_overflow = (
                struct.pack('<iiiIII', 0, 9, 1, 0x7F800001, (1 << 16) | 0, 0)
                + b'\x07'      # empty varid
                + b'\x57A'     # claims 5-byte ref string with only 1 byte left
            )
            with gzip.open(bad_shared_bcf, 'wb') as gz:
                gz.write(b'BCF\x02\x02' + struct.pack('<I', len(header_text)) + header_text)
                gz.write(struct.pack('<II', len(shared_overflow), 0) + shared_overflow)

            with self.assertRaises(ValueError):
                list(BcfReader(bad_shared_bcf))

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
            self.assertEqual(set(var_pybcf.info), set(var_pysam.info))
            for field in var_pysam.info:
                val = True if var_pysam.info[field] is None else var_pysam.info[field]
                self.assertEqual(val, var_pybcf.info[field])
            
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
                
