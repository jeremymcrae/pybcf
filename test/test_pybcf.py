
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
    is_nan = geno[:, 0] == -1
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
    
    def test_without_genotypes(self):
        ''' check this package matches pysam for info fields for BCF without genotypes
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
            
            # check all the info fields match
            for field in var_pysam.info:
                self.assertEqual(var_pysam.info[field], var_pybcf.info[field],)
    
    def test_with_genotypes(self):
        ''' check this package matches pysam for info fields for BCF
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
    
