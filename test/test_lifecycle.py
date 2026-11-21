
import gc
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
from pybcf import BcfReader

from .helpers import write_bcf

class TestLifecycle(unittest.TestCase):
    ''' class to make sure readers, and the objects made from them, are created and closed correctly
    '''

    sample_bcf = Path(__file__).parent / 'data' / '1000G.shrunk.bcf'

    def test_wrappers_keep_parents_alive(self):
        ''' check header, info and sample objects work when nothing else holds their reader or variant
        '''
        # the readers and variants here are temporaries, which each result has to keep alive
        contigs = BcfReader(self.sample_bcf).header.contigs
        self.assertGreater(len(contigs), 0)
        info_keys = list(next(BcfReader(self.sample_bcf)).info)
        self.assertGreater(len(info_keys), 0)
        gt = next(BcfReader(self.sample_bcf)).samples['GT']
        self.assertEqual(gt.ndim, 2)

        # variants and the header stay usable after the `with` block closes the reader
        with BcfReader(self.sample_bcf) as bcf:
            var = next(bcf)
        self.assertGreater(len(var.info), 0)
        self.assertEqual(len(var.info), len(list(var.info)))
        self.assertEqual(dict(var.info.items()), {k: var.info[k] for k in var.info})
        self.assertEqual(len(var.samples), len(list(var.samples)))
        sample_items = list(var.samples.items())
        self.assertEqual([k for k, _ in sample_items], list(var.samples))
        np.testing.assert_array_equal(var.samples.get('GT'), gt)
        self.assertEqual(bcf.header.contigs, contigs)

    def test_closed_reader_raises(self):
        ''' check a closed reader raises ValueError when asked for more variants
        '''
        with BcfReader(self.sample_bcf) as bcf:
            chrom = next(bcf).chrom
        with self.assertRaises(ValueError):
            next(bcf)

        # close() can be called more than once
        bcf = BcfReader(self.sample_bcf)
        bcf.close()
        bcf.close()
        with self.assertRaises(ValueError):
            next(bcf)
        with self.assertRaises(ValueError):
            bcf.fetch(chrom)

    def test_wrappers_cannot_be_created_directly(self):
        ''' check the header, variant, info and sample classes raise TypeError if created by hand
        '''
        bcf = BcfReader(self.sample_bcf)
        var = next(bcf)
        for cls in (type(bcf.header), type(var.info), type(var.samples)):
            with self.assertRaises(TypeError):
                cls()
        with self.assertRaises(TypeError):
            type(var)(bcf)
        # `__new__` skips `__init__`, so the wrappers must also refuse to exist without
        # the object they read from
        for cls in (type(var.info), type(var.samples), type(var)):
            for args in ((), (None,)):
                with self.assertRaises(TypeError):
                    cls.__new__(cls, *args)

        # A variant made with `__new__` is empty, and has no contig to name if the
        # header lacks contigs
        with tempfile.TemporaryDirectory() as tmpdir:
            no_contigs = Path(tmpdir) / 'no_contigs.bcf'
            write_bcf(no_contigs, '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            empty = type(var).__new__(type(var), BcfReader(no_contigs))
            self.assertEqual((empty.pos, empty.ref, empty.alts, empty.filter), (0, '', None, []))
            self.assertEqual((len(empty.info), len(empty.samples)), (0, 0))
            with self.assertRaises(IndexError):
                empty.chrom

    @unittest.skipUnless(Path('/proc/self/fd').exists(), 'needs /proc to count open files')
    def test_discarded_reader_closes_file(self):
        ''' check a reader with no references left closes its file at once

        This fails if the reader and its header or variants refer to each other in
        a cycle, as the file then stays open until the garbage collector runs.
        '''
        gc.collect()
        gc.disable()
        try:
            before = len(os.listdir('/proc/self/fd'))
            for _ in range(25):
                r = BcfReader(self.sample_bcf)
                r.header.contigs
                r.samples
                var = next(r)
                list(var.info)
                var.samples['GT']
                del r, var
            after = len(os.listdir('/proc/self/fd'))
        finally:
            gc.enable()
        self.assertEqual(before, after)
