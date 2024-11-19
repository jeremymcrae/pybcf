
import unittest
from pathlib import Path

from pybcf import BcfReader
from pysam import VariantFile

class TestIndexed(unittest.TestCase):
    ''' class to make sure BcfReader works correctly
    '''
    
    hapmap = Path(__file__).parent / 'data' / 'hapmap_3.3.hg38.shrunk.bcf'
    test_cases = Path(__file__).parent / 'data' / 'bcfv2.2.bcf'
    
    def test_fetch_exists(self):
        ''' check this package has a fetch method
        '''
        bcf = BcfReader(self.hapmap)
        self.assertTrue(hasattr(bcf, 'fetch') and callable(getattr(bcf, 'fetch')))
    
    def test_matches_pysam(self):
        ''' check this fetch gets the same variants as pysam fetch
        '''
        bcf = BcfReader(self.hapmap)
        vcf = VariantFile(self.hapmap)
        
        sites = [(x.chrom, x.pos) for x in bcf]
        chroms = set(x[0] for x in sites)
        self.assertTrue(len(chroms) == 1)
        chrom = next(iter(chroms))
        
        _, pos_50th = sites[len(sites) // 2]
        _, pos_75th = sites[(len(sites) // 4) * 3]
        
        # check that we get variants restricted to a region after fetching
        pybcf_fetched = [(x.chrom, x.pos) for x in bcf.fetch(chrom, pos_50th, pos_75th)]
        pysam_fetched = [(x.chrom, x.pos) for x in vcf.fetch(chrom, pos_50th-1, pos_75th)]
        self.assertEqual(pybcf_fetched, pysam_fetched)
    
    def test_fetch_after_bcf_iterated(self):
        ''' check that we can fetch variants after we reach the end of the bcf
        '''
        bcf = BcfReader(self.hapmap)
        chrom = next(bcf).chrom
        for x in bcf:
            pass
        
        # if we try to iterate now, there are no more variants
        self.assertEqual(list(bcf), [])
        
        # but we can still fetch variants even after hitting the end of the bcf
        self.assertTrue(len(list(bcf.fetch(chrom))) > 0)
    
    def test_fetch_one_chrom(self):
        ''' check that if we fetch one chrom, none of the other chroms are included
        '''
        bcf = BcfReader(self.test_cases)
        
        sites = [(x.chrom, x.pos) for x in bcf]
        chroms = set(x[0] for x in sites)
        self.assertTrue(len(chroms) > 1)
        chrom = sorted(chroms)[0]
        
        # we only get variants from the fetched chromosome
        fetched = [(x.chrom, x.pos) for x in bcf.fetch(chrom)]
        self.assertEqual([x for x in sites if x[0] <= chrom], fetched)
    
    def test_fetch_at_bcf_start(self):
        ''' check fetching a region at the start of a bcf
        '''
        bcf = BcfReader(self.test_cases)
        first = next(bcf)
        fetched = list(bcf.fetch(first.chrom, first.pos, first.pos))
        self.assertTrue(len(fetched) == 1)
    
    def test_iteration_stops_after_fetch(self):
        ''' check that we can't iterate after fetching variants
        '''
        bcf = BcfReader(self.test_cases)
        first = next(bcf)
        fetched = list(bcf.fetch(first.chrom, first.pos, first.pos))
        self.assertTrue(len(fetched) > 0)
        
        # check if we fetch variants, we can't iterate afterwards, even if the
        # fetching only finds variants early in the bcf
        iterated = [x for x in bcf]
        self.assertTrue(len(iterated) == 0)
    
    def test_fetch_later_bzgf_block(self):
        ''' check we can fetch that requires seeking to a later bzgf block
        
        Specifcally, fetch for a region that starts in a bzgf block different 
        to the block containing the first variant.
        '''
        
        # use the hapmap file, which is large enough to require different bzgf 
        # blocks, and variant positions which cover different index bins
        bcf = BcfReader(self.hapmap)
        
        sites = [(x.chrom, x.pos) for x in bcf]
        chroms = set(x[0] for x in sites)
        self.assertTrue(len(chroms) == 1)
        chrom = next(iter(chroms))
        
        _, pos_75th = sites[(len(sites) // 4) * 3]
        fetched = [(x.chrom, x.pos) for x in bcf.fetch(chrom, pos_75th)]
        self.assertEqual([x for x in sites if x[0] <= chrom and x[1] >= pos_75th], fetched)
    
    def test_fetch_middle_of_chrom(self):
        ''' check we can fetch variants in the middle of a chromosome
        '''
        bcf = BcfReader(self.hapmap)
        
        sites = [(x.chrom, x.pos) for x in bcf]
        chroms = set(x[0] for x in sites)
        self.assertTrue(len(chroms) == 1)
        chrom = next(iter(chroms))
        
        _, pos_25th = sites[len(sites) // 4]
        _, pos_50th = sites[len(sites) // 2]
        
        # check that we get variants restricted to a region after fetching
        fetched = [(x.chrom, x.pos) for x in bcf.fetch(chrom, pos_25th, pos_50th)]
        self.assertEqual([x for x in sites if pos_25th <= x[1] <= pos_50th], fetched)
    
    def test_raise_error_for_absent_chrom(self):
        ''' check we raise error if we try for a chromosome that doesn't exist
        '''
        bcf = BcfReader(self.test_cases)
        sites = [(x.chrom, x.pos) for x in bcf]
        chroms = set(x[0] for x in sites)
        
        missing_chrom = 'chrSDIHG'
        self.assertTrue(missing_chrom not in chroms)
        with self.assertRaises(ValueError):
            bcf.fetch(missing_chrom)
    
    def test_raise_error_for_swapped_coordinates(self):
        ''' check we raise error if the start position is after the end
        '''
        bcf = BcfReader(self.test_cases)
        chrom = next(bcf).chrom
        
        with self.assertRaises(ValueError):
            bcf.fetch(chrom, 10, 9)
        
        # but we can still iterate, since the fetch didn't succeed
        iterated = list(bcf)
        self.assertTrue(len(iterated) > 0)
    
    def test_fetch_chrom_without_variants(self):
        ''' check we don't raise an error if fetching for a chrom that exist
        '''
        # check we get an iterable without variants if we fetch a chrom that 
        # exists in the index, but doesn't have any variants in the VCF
        
        bcf = BcfReader(self.test_cases)
        sites = [(x.chrom, x.pos) for x in bcf]
        chroms_for_variants = set(x[0] for x in sites)
        chroms_for_bcf = bcf.header.contigs
        
        without_variants = set(chroms_for_bcf) - set(chroms_for_variants)
        chrom_without_variants = sorted(without_variants)[0]
        
        variants = bcf.fetch(chrom_without_variants)
        self.assertTrue(isinstance(variants, BcfReader))
        
        variants = list(variants)
        self.assertTrue(len(variants) == 0)
    
    def test_fetch_region_without_variants(self):
        ''' check we don't get any variants if we fetch a region without variants
        '''
        
        bcf = BcfReader(self.test_cases)
        sites = [(x.chrom, x.pos) for x in bcf]
        last_chrom, last_pos = sites[-1]
        
        variants = list(bcf.fetch(last_chrom, last_pos, last_pos))
        self.assertTrue(len(variants) == 1)
        
        variants = list(bcf.fetch(last_chrom, last_pos+1, last_pos+1))
        self.assertTrue(len(variants) == 0)
    
    def test_multiple_fetches(self):
        ''' check if we can fetch multiple times and get variants
        '''
        
        bcf = BcfReader(self.test_cases)
        sites = [(x.chrom, x.pos) for x in bcf]
        last_chrom, last_pos = sites[-1]
        
        variants = list(bcf.fetch(last_chrom, last_pos, last_pos))
        self.assertTrue(len(variants) == 1)
        
        variants = list(bcf.fetch(last_chrom, last_pos, last_pos))
        self.assertTrue(len(variants) == 1)
    
    def test_raise_error_without_index(self):
        ''' check we get an error if we fetch from a bcf without an index file
        '''
        index_path = self.test_cases.with_suffix('.bcf.csi')
        temp_path = self.test_cases.with_suffix('.bcf.BACKUP')
        
        try:
            index_path.rename(temp_path)
            bcf = BcfReader(self.test_cases)
            chrom = next(bcf).chrom
            with self.assertRaises(ValueError):
                bcf.fetch(chrom)
        finally:
            temp_path.rename(index_path) 
