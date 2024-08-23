
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
            self.assertEqual(set(var_pybcf.info), set(var_pysam.info))
            
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
            self.assertEqual(set(var_pybcf.info), set(var_pysam.info))
            
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
            
            self.assertEqual(set(var_pysam.format), set(var_pybcf.samples))
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
            
            self.assertEqual(set(var_pysam.format), set(var_pybcf.samples))
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
                
