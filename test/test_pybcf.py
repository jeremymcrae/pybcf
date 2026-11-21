
import ast
import gzip
import os
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
                
