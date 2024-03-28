
import logging
from pathlib import Path

from libcpp.string cimport string

cdef extern from 'bcf.h' namespace 'bcf':
    cdef cppclass BCF:
        # declare class constructor and methods
        BCF(string path) except +
        # Variant & next()

cdef class BcfReader:
    ''' class to open bcf files from disk, and access variant data within
    '''
    cdef BCF * thisptr
    cdef string path
    def __cinit__(self, path):
        self.path = str(path) if isinstance(path, Path) else path
        self.path = self.path.encode('utf8')
        
        logging.debug(f'opening BcfReader from {self.path.decode("utf")}')
        self.thisptr = new BCF(self.path)
    
    def __dealloc__(self):
        if self.is_open:
          del self.thisptr
        self.is_open = False
    
    def __repr__(self):
        return f'BcfReader("{self.path.decode("utf8")}")'
    
    def __iter__(self):
        return self
    
    def __next__(self):
        ''' iterate through all variants in the bcf file
        '''
        raise NotImplementedError
    
    @property
    def header(self):
      ''' get header info from bcf file
      '''
      raise NotImplementedError
    
    @property
    def samples(self):
      ''' get list of samples in the bcf file
      '''
      raise NotImplementedError
    
    def fetch(self, chrom, start=None, stop=None):
        ''' fetches all variants within a genomic region
        
        Args:
            chrom: chromosome that variants must be on
            start: start nucleotide of region. If None, gets all variants on chromosome
            stop: end nucleotide of region. If None, gets variants with positions after start
        
        Yields:
            BcfVars for variants within the genome region
        '''
        raise NotImplementedError
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.__dealloc__()
        return False
