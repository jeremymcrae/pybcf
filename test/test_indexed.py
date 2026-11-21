
import gzip
import itertools
import multiprocessing
import os
import shutil
import struct
import sys
import tempfile
import unittest
from pathlib import Path

from pybcf import BcfReader
from pysam import VariantFile, tabix_index

from .helpers import record, write_bcf

class TestIndexed(unittest.TestCase):
    ''' class to make sure BcfReader works correctly
    '''
    
    hapmap = Path(__file__).parent / 'data' / 'hapmap_3.3.hg38.shrunk.bcf'
    test_cases = Path(__file__).parent / 'data' / 'bcfv2.2.bcf'
    many_samples = Path(__file__).parent / 'data' / '1000G.shrunk.bcf'
    
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
        with tempfile.TemporaryDirectory() as tmpdir:
            # a copy of the bcf, without the index file that sits next to the original
            path = Path(tmpdir) / 'no_index.bcf'
            shutil.copy(self.test_cases, path)
            bcf = BcfReader(path)
            chrom = next(bcf).chrom
            with self.assertRaises(ValueError):
                bcf.fetch(chrom)

            # an index path given by the user has to exist
            with self.assertRaises(ValueError):
                BcfReader(self.hapmap, index_path=Path(tmpdir) / 'missing.csi')

    def test_invalid_index_files_raise(self):
        ''' check index files that are not CSI, or have impossible values, are rejected
        '''
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_bcf = Path(tmpdir) / 'corrupt_idx.bcf'
            tmp_csi = Path(tmpdir) / 'corrupt_idx.bcf.csi'
            shutil.copy(self.test_cases, tmp_bcf)
            for contents in (
                b'NOT_CSI',
                # min_shift and depth that are too large, alone or together
                b'CSI\x01' + struct.pack('<iiii', 14, 0x40000000, 0, 0),
                b'CSI\x01' + struct.pack('<iiii', 0, 20, 0, 0),
                b'CSI\x01' + struct.pack('<iiii', 0, 11, 0, 0),
                # one contig, with a negative number of bins
                b'CSI\x01' + struct.pack('<iiiii', 14, 5, 0, 1, -1),
            ):
                with gzip.open(tmp_csi, 'wb') as gz:
                    gz.write(contents)
                with self.assertRaises(ValueError):
                    BcfReader(tmp_bcf)

    def test_fetch_after_changing_directory(self):
        ''' check a reader opened with a relative path can fetch after a change of directory
        '''
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                os.chdir(self.hapmap.parent)
                bcf = BcfReader(self.hapmap.name)
                chrom = next(bcf).chrom
                os.chdir(tmpdir)
                self.assertGreater(len(list(bcf.fetch(chrom))), 0)
            finally:
                os.chdir(cwd)

    @unittest.skipUnless(hasattr(os, 'fork'), 'needs forked processes')
    def test_fetch_in_forked_process(self):
        ''' check a forked copy of a reader can fetch, without disturbing the original

        The two copies share one open file. The file is much larger than the 64 kb
        read at a time, so the original reads more of it after the copy has fetched.
        '''
        sites = [(x.chrom, x.pos) for x in BcfReader(self.many_samples)]
        bcf = BcfReader(self.many_samples)
        first = next(bcf)

        def fetch_last_site():
            fetched = [(x.chrom, x.pos) for x in bcf.fetch(*sites[-1])]
            sys.exit(0 if fetched == sites[-1:] else 1)

        child = multiprocessing.get_context('fork').Process(target=fetch_last_site)
        child.start()
        child.join()
        self.assertEqual(child.exitcode, 0)

        rest = [(x.chrom, x.pos) for x in bcf]
        self.assertEqual([(first.chrom, first.pos)] + rest, sites)

    def test_fetch_without_bgzf_raises(self):
        ''' check fetching from a BCF that isn't bgzipped raises an error, even with an index
        '''
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'uncompressed.bcf'
            path.write_bytes(gzip.open(self.test_cases).read())
            bcf = BcfReader(path, index_path=self.test_cases.with_suffix('.bcf.csi'))
            with self.assertRaisesRegex(ValueError, 'bgzipped'):
                bcf.fetch('chr1')

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

    def test_csi_with_tiny_bins(self):
        ''' check a CSI with 1 bp leaf bins at maximum depth fetches the same variants quickly
        '''
        raw = gzip.open(self.test_cases.with_suffix('.bcf.csi'), 'rb').read()
        min_shift, depth, l_aux = struct.unpack_from('<iii', raw, 4)
        bin_limit = ((1 << (3 * (depth + 1))) - 1) // 7
        pos = 16 + l_aux
        n_ref = struct.unpack_from('<i', raw, pos)[0]
        pos += 4

        # Move every real chunk into the root bin of a min_shift=0, depth=10 index, so a
        # whole-contig query overlaps over a billion possible leaf bins
        tiny = b'CSI\x01' + struct.pack('<iiii', 0, 10, 0, n_ref)
        for _ in range(n_ref):
            n_bins = struct.unpack_from('<i', raw, pos)[0]
            pos += 4
            chunks = []
            for _ in range(n_bins):
                bin_idx, _, n_chunks = struct.unpack_from('<IQi', raw, pos)
                pos += 16
                if bin_idx < bin_limit:  # skip htslib's metadata pseudo-bin
                    chunks.append(raw[pos:pos + 16 * n_chunks])
                pos += 16 * n_chunks
            payload = b''.join(chunks)
            if payload:
                tiny += struct.pack('<iIQi', 1, 0, 0, len(payload) // 16) + payload
            else:
                tiny += struct.pack('<i', 0)

        with tempfile.TemporaryDirectory() as tmpdir:
            tiny_csi = Path(tmpdir) / 'tiny.csi'
            with gzip.open(tiny_csi, 'wb') as gz:
                gz.write(tiny)

            real = BcfReader(self.test_cases)
            tiny_bins = BcfReader(self.test_cases, index_path=tiny_csi)
            for chrom in real.header.contigs:
                expected = [(v.chrom, v.pos) for v in real.fetch(chrom)]
                self.assertEqual([(v.chrom, v.pos) for v in tiny_bins.fetch(chrom)], expected)
            self.assertEqual(len(list(tiny_bins.fetch('chr1', 1, 3))), 3)

    def test_fetch_starting_in_empty_window(self):
        ''' check a fetch that starts where there are no variants begins reading nearby

        A record that crosses a top-level bin boundary is stored in the root bin of
        the index. A query starting in a 16 kb window without variants used to begin
        reading from that record, however far upstream it was. The BGZF blocks in
        between are corrupted here, so reading through them raises an error.
        '''
        boundary = 1 << 26
        gap_start, gap_end = 100_000_000, 100_200_000
        positions = (list(range(1_000_000, boundary - 1000, 1_000_000))
                     + [boundary - 500]
                     + list(range(67_200_000, gap_start, 2000))
                     + list(range(gap_end, gap_end + 200_000, 2000)))
        header = ('##fileformat=VCFv4.2\n'
                  '##contig=<ID=chr1,length=248956422>\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        # the record at the boundary spans 1 kb, the others are a single base
        records = [record(pos0=pos - 1, rlen=1000 if pos == boundary - 500 else 1)
                   for pos in positions]
        regions = [(gap_start + 50_000, gap_start + 60_000),
                   (gap_start + 50_000, gap_end + 50_000)]

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'gap.bcf'
            write_bcf(path, header, records)
            tabix_index(str(path), preset='bcf', csi=True, force=True)

            with VariantFile(str(path)) as vcf:
                expected = [[v.pos for v in vcf.fetch('chr1', start - 1, end)]
                            for start, end in regions]
            self.assertEqual(expected[0], [])
            self.assertEqual(len(expected[1]), 26)

            # find the uncompressed bytes from the boundary record to a record 100 kb
            # before the gap, then corrupt every BGZF block inside that. Opening the
            # file reads 64 kb ahead, so blocks that close to the start are left alone.
            starts = list(itertools.accumulate([9 + len(header) + 1] + [len(x) for x in records]))
            lo = starts[positions.index(boundary - 500) + 1] + 65536
            hi = starts[positions.index(gap_start - 100_000)]
            raw = bytearray(path.read_bytes())
            offset = uncompressed = n_corrupted = 0
            while offset < len(raw):
                size = struct.unpack_from('<H', raw, offset + 16)[0] + 1
                length = struct.unpack_from('<I', raw, offset + size - 4)[0]
                if lo <= uncompressed and uncompressed + length <= hi:
                    raw[offset + size // 2] ^= 0xff
                    n_corrupted += 1
                offset += size
                uncompressed += length
            self.assertGreater(n_corrupted, 0)
            path.write_bytes(raw)

            bcf = BcfReader(path)
            for (start, end), sites in zip(regions, expected):
                self.assertEqual([v.pos for v in bcf.fetch('chr1', start, end)], sites)

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

    def test_fetch_telomeric_position(self):
        ''' check records at POS 0 (the telomere) are fetched when start is None or 0
        '''
        header = ('##fileformat=VCFv4.2\n'
                  '##contig=<ID=chr1,length=1000,IDX=0>\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'telomeric.bcf'
            # POS 0 is stored as a 0-based position of -1
            write_bcf(path, header, [record(pos0=-1, body=b'\x07\x17N\x00'), record()])
            tabix_index(str(path), preset='bcf', csi=True, force=True)

            bcf = BcfReader(path)
            self.assertEqual([v.pos for v in bcf.fetch('chr1')], [0, 10])
            self.assertEqual([v.pos for v in bcf.fetch('chr1', 0, 5)], [0])
