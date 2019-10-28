
### ToDos
- [chemfp](https://chemfp.readthedocs.io/en/latest/using-api.html)
- add caching of fingerprints
- test on server
- check [fragmentation algorithm](https://www.rdkit.org/docs/source/rdkit.Chem.Fraggle.FraggleSim.html#fragmentation-algorithm) on results?

## Install rdkit using anaconda
- install [rdkit](https://www.rdkit.org/docs/Install.html#cross-platform-under-anaconda-python-fastest-install) using anaconda
- add [postgres](https://www.rdkit.org/docs/Install.html#installing-and-using-postgresql-and-the-rdkit-postgresql-cartridge-from-a-conda-environment) installation using anaconda

## Solution using rdkit functionality
1. Create Postgres Database [Cartridge](https://www.rdkit.org/docs/Cartridge.html)
2. Use [SimilarityScreener](https://www.rdkit.org/docs/source/
rdkit.Chem.Fingerprints.SimilarityScreener.html#rdkit.Chem.Fingerprints.SimilarityScreener.SimilarityScreener)

### Data

Download the Enamine REAL files (log in needed):
```bash
>>> wc -l 2019q1-2_Enamine_REAL_723M_SMILES_Part_01.smiles
60262532 2019q1-2_Enamine_REAL_723M_SMILES_Part_01.smiles
```

List of available Fingerprints in rdkit: [here](https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-fingerprints)

### Similarity metrics between Fingerprints
- Available similarity metrics include Tanimoto, Dice, Cosine, Sokal, Russel, Kulczynski, McConnaughey, and Tversky.


## Excecution Time
- Read 100000 lines
- Compare to 1 reference molecule

```
C:\Users\kzl465\Documents\rotation_3>ipython search_database.py -i
Python 3.6.9 |Anaconda, Inc.| (default, Jul 30 2019, 14:00:49) [MSC v.1915 64 bit (AMD64)]
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.
Fingerprint lenght: 2048
Script is executed from: C:\Users\kzl465\Documents\rotation_3
Found 1 input files:
- data\2019q1-2_Enamine_REAL_723M_SMILES_Part_01.smiles
[11:43:50] SMILES Parse Error: syntax error for input: 'smiles'
Failed to read line 0
[]

         703686 function calls in 223.656 seconds

   Ordered by: standard name

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000    0.000    0.000 :0(_getdefaultlocale)
     1809    0.047    0.000    0.047    0.000 :0(charmap_decode)
        1    0.000    0.000  223.656  223.656 :0(exec)
       10    0.000    0.000    0.000    0.000 :0(finditer)
        5    0.000    0.000    0.000    0.000 :0(flush)
        5    0.000    0.000    0.000    0.000 :0(len)
        1    0.000    0.000    0.000    0.000 :0(open)
        3    0.000    0.000    0.000    0.000 :0(print)
        1    0.000    0.000    0.000    0.000 :0(setprofile)
   100001    0.516    0.000    0.516    0.000 :0(split)
        5    0.000    0.000    0.000    0.000 :0(write)
        1    0.000    0.000  223.656  223.656 <string>:1(<module>)
   100001    0.484    0.000    1.328    0.000 FingerprintMols.py:281(__init__)
   100001    0.406    0.000    0.406    0.000 FingerprintMols.py:286(_fingerprinterInit)
   100001    0.234    0.000    0.234    0.000 FingerprintMols.py:312(_screenerInit)
   100001    0.203    0.000    0.203    0.000 FingerprintMols.py:323(_clusterInit)
   100001  180.016    0.002  181.344    0.002 FingerprintMols.py:63(FingerprintMol)
   100000    1.938    0.000    1.938    0.000 __init__.py:33(FingerprintSimilarity)
        1    0.000    0.000    0.000    0.000 _bootlocale.py:11(getpreferredencoding)
        5    0.000    0.000    0.000    0.000 ansitowin32.py:160(write)
        5    0.000    0.000    0.000    0.000 ansitowin32.py:177(write_and_convert)
        5    0.000    0.000    0.000    0.000 ansitowin32.py:193(write_plain_text)
        5    0.000    0.000    0.000    0.000 ansitowin32.py:245(convert_osc)
        5    0.000    0.000    0.000    0.000 ansitowin32.py:40(write)
        1    0.000    0.000    0.000    0.000 codecs.py:259(__init__)
     1809    0.031    0.000    0.078    0.000 cp1252.py:22(decode)
        1    0.000    0.000  223.656  223.656 profile:0(print(check_file(filename)); print())
        0    0.000             0.000          profile:0(profiler)
        1   39.781   39.781  223.656  223.656 search_database.py:60(check_file)


Time elapsed 00:03:50
```


## Generators
- Since Python 3.7 Generators are stop using `return` instead of `stopIteration`
