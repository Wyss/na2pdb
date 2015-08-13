# na2pdb
Convert DNA and RNA to PDB files.  Useful for simulation input data.
Conceived to be the base of a [cadnano 2.5](https://github.com/cadnano/cadnano2.5) 
*.nno converter, which should show up in that repo soon.

dependencies: ProDy, numpy

Should be Python 2/3 compatible

RNA is not supported yet, but totally could be.

Base data pdb's were sourced from Avogadro 1.1.1/OpenBabel Output and an 
effort is made to come as close to the spatial output for double stranded DNA as
would be generated from Avogadro 1.1.1.

```python
""" create a strand that double backs on itself
"""
from na2pdb import AtomicSequence
as = AtomicSequence("ACGTACGT", name="not_useful")
as.transformBases(4, 8, 0, 0, 0, False)
# 1. Get base separation
as.linearize()
# 2. do all rotations
as.applyReverseQueue()
as.applyTwist()
# 3. move to position
as.applyTransformQueue()

out_file = 'test_file.pdb'
as.toPDB(out_file)
```

