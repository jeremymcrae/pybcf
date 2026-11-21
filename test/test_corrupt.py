
import gzip
import resource
import struct
import tempfile
import unittest
from pathlib import Path

from pybcf import BcfReader

from .helpers import record, write_bcf

class TestCorrupt(unittest.TestCase):
    ''' class to make sure missing, truncated or corrupt BCFs raise errors
    '''

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
