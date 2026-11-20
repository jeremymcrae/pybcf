
import gzip
import shutil
import struct
import tempfile
import unittest
from pathlib import Path

from pybcf import BcfReader
from pysam import VariantFile, tabix_index

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
        with self.assertRaises(ValueError):
            bcf.fetch(chrom, -1, 10)
        
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

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(ValueError):
                BcfReader(self.hapmap, index_path=Path(tmpdir) / 'missing.csi')

            tmp_bcf = Path(tmpdir) / 'corrupt_idx.bcf'
            tmp_csi = Path(tmpdir) / 'corrupt_idx.bcf.csi'
            shutil.copy(self.test_cases, tmp_bcf)
            with gzip.open(tmp_csi, 'wb') as gz:
                gz.write(b'NOT_CSI')
            with self.assertRaises(ValueError):
                BcfReader(tmp_bcf)

    def test_fetch_empty_chrom_on_fresh_reader(self):
        ''' check fetching an empty contig works without prior iteration
        '''
        bcf = BcfReader(self.test_cases)
        self.assertEqual(list(bcf.fetch('chr3')), [])

    def test_csi_with_aux_data_and_fewer_contigs(self):
        ''' check CSI files with non-zero l_aux or fewer contigs than header work safely
        '''
        orig_csi = self.test_cases.with_suffix('.bcf.csi')
        raw = gzip.open(orig_csi, 'rb').read()
        magic = raw[:4]
        min_shift, depth, l_aux = struct.unpack('<iii', raw[4:16])
        self.assertEqual(l_aux, 0)
        rest = raw[16:]
        n_ref = struct.unpack('<i', rest[:4])[0]

        aux_bytes = b'AUX_METADATA_16B'
        # Rebuild CSI with l_aux = 16 and n_ref = 1 (fewer contigs than BCF header)
        # Keep only the first reference's bin table from rest[4:]
        pos = 4
        n_bins = struct.unpack('<I', rest[pos:pos + 4])[0]
        pos += 4
        for _ in range(n_bins):
            _, _, n_chunks = struct.unpack('<IQI', rest[pos:pos + 16])
            pos += 16 + n_chunks * 16
        first_contig_bytes = rest[4:pos]

        modified_raw = (
            magic
            + struct.pack('<iii', min_shift, depth, len(aux_bytes))
            + aux_bytes
            + struct.pack('<i', 1)
            + first_contig_bytes
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_bcf = Path(tmpdir) / 'test.bcf'
            tmp_csi = Path(tmpdir) / 'test.bcf.csi'
            shutil.copy(self.test_cases, tmp_bcf)
            with gzip.open(tmp_csi, 'wb') as gz:
                gz.write(modified_raw)

            bcf = BcfReader(tmp_bcf)
            self.assertEqual(len(list(bcf.fetch('chr1', 1, 3))), 3)
            # chr3 has contig_id >= indices.size() in this truncated CSI
            self.assertEqual(list(bcf.fetch('chr3')), [])

    def test_fetch_spanning_deletion_and_none_bounds(self):
        ''' check fetch includes upstream spanning deletions and supports start=None/stop=None
        '''
        with tempfile.TemporaryDirectory() as tmpdir:
            bcf_path = Path(tmpdir) / 'spanning.bcf'
            with VariantFile(bcf_path, 'wb') as out:
                out.header.contigs.add('chr1', length=1000000)
                # Deletion spanning 100..40099 (crosses 16 kb bin boundaries)
                out.write(out.new_record(contig='chr1', start=99, alleles=('A' * 40000, 'A')))
                # Short SNP at 200 (after deletion start, but upstream of 36000)
                out.write(out.new_record(contig='chr1', start=199, alleles=('C', 'G')))
                # SNP inside 36000..39000
                out.write(out.new_record(contig='chr1', start=36999, alleles=('G', 'T')))
                # SNP at exact 16 kb bin boundary (0-based 49152)
                out.write(out.new_record(contig='chr1', start=49152, alleles=('T', 'A')))
            tabix_index(str(bcf_path), preset='bcf', csi=True, force=True)

            vcf = VariantFile(bcf_path)
            bcf = BcfReader(bcf_path)

            for s, e in [(36000, 39000), (49153, 50000)]:
                expected = [v.pos for v in vcf.fetch('chr1', s - 1, e)]
                actual = [v.pos for v in bcf.fetch('chr1', s, e)]
                self.assertEqual(actual, expected)

            self.assertEqual([v.pos for v in bcf.fetch('chr1', None, None)], [100, 200, 37000, 49153])
            self.assertEqual([v.pos for v in bcf.fetch('chr1', 36000, None)], [100, 37000, 49153])
            self.assertEqual([v.pos for v in bcf.fetch('chr1', None, 300)], [100, 200])

