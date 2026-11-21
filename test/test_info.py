
import struct
import tempfile
import unittest
from pathlib import Path

import numpy as np
from pysam import VariantFile
from pybcf import BcfReader

from .helpers import record, write_bcf

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

class TestInfo(unittest.TestCase):
    ''' class to make sure INFO fields are decoded correctly

    Tests of values that INFO and FORMAT fields store the same way, such as missing
    and reserved values, are also here, and check both.
    '''

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
