
### pybcf

This is a package for reading BCF files (binary call format). Similar to pysam/cyvcf2, 
but limited to reading from BCF files only, and oriented around getting genotype 
data into numpy arrays quickly. 

This performs relatively better at larger sample sizes e.g. it's ~15X faster than
pysam with 100 samples, but improves to ~200X faster than pysam with 100k samples. 

```py
from pybcf import BcfReader

bcf = BcfReader(bcf_path)

sample_ids = bcf.samples
contigs = bcf.header.contigs

for var in bcf:
    # the usual attributes are available e.g.
    # var.chrom, var.pos, var.ref, var.alts, var.info['AF']
    
    # sample data is accessed as numpy arrays via the format keys
    keys = list(var.samples)
    genotypes = var.samples['GT']  # as n x 2 numpy array, missing=nan

# or fetch from random regions if bcf is indexed. Positions are 1-based and
# include both ends (pysam's fetch is 0-based and excludes the end)
for var in bcf.fetch('chr1', 10000, 200000):
    print(var.chrom, var.pos)

# you can read from stdin, but that doesn't allow fetching random regions
import sys
bcf = BcfReader(sys.stdin) 
```

### Limitations
 - `fetch()` needs a bgzip-compressed BCF with a `.csi` index
 - a reader has one position in the file, so `fetch()` inside a loop over the
   same reader ends that loop. Open a second `BcfReader` to do both at once
 - extracting info fields is a little slow
