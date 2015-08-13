# na2pdb
Convert DNA sequencesto PDB files and spatially wrangle them.  
Useful for simulation input data.
Conceived to be the base of a [cadnano 2.5](https://github.com/cadnano/cadnano2.5) 
*.nno converter, which should show up in that repo soon.

dependencies: [ProDy](http://prody.csb.pitt.edu/index.html), numpy

Should be Python 2/3 compatible

RNA is not supported yet, but totally could be.

Base data pdb's were sourced from [Avogadro 1.1.1](http://avogadro.cc/wiki/Main_Page) / 
OpenBabel output and an effort is made to come as close to the spatial output 
for double stranded DNA as would be generated from Avogadro 1.1.1.

```python
""" create a strand that double backs on itself
"""
from na2pdb import AtomicSequence
aseq = AtomicSequence("ACGTACGT", name="not_useful")
aseq.transformBases(4, 8, 0, 0, 0, False)
# 1. Get base separation
aseq.linearize()
# 2. do all rotations
aseq.applyReverseQueue()
aseq.applyTwist()
# 3. move to position
aseq.applyTransformQueue()

out_file = 'test_file.pdb'
aseq.toPDB(out_file)
```

