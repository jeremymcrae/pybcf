
### pybcf

This is a package for reading BCF files (binary call format). Similar to pysam/cyvcf2, but limited to reading from BCF files only, and oriented around getting genotype data into numpy arrays quickly.

'''py
from pybcf import BcfReader

bcf = BcfReader(bcf_path)

for var in bcf:
    # the usual attributes are available e.g.
    # var.chrom, var.pos, var.ref, var.alts, var.info['AF']
    
    # get genotypes as numpy array
    genotypes = var.samples['GT']  # as numpy array, missing=-1
'''

### Limitations
 - can't fetch variants from random regions yet
 - not extensively tested yet
 - extracting info fields is slow at the moment
