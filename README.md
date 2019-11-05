# Similarity Search of Enamine REAL
In a first approach multitasking was tried to effectively find similar fingerprints. However, this did not scale good enough. Using the approach of caching the fingerprints before quering the ENAMINE REAL database brought down the search time on 12 CPUs to less than 30mins.

The python scripts can be executed using two shell scripts.

```
user@server: ls bin/
create_blobs.sh  
run_search.sh  
```
After the setup (see below), execute both scripts using up to 12 CPUs from the folder you cloned this repository to, e.g. `~/search_enamine_real_db`:
```
bash bin/create_blobs.sh --cpus 12
bash bin/run_search --cpus 12 --query 'C1C=CC....' --threshold 0.8 --force false
``` 
> NOTE: This will create blobs of roughly 590GBs!

Currently the query uses the rdkit fingerprint (daylight), but you can in principle pick any in [rdkit available fingerprint](https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-fingerprints). As metric the Tanimoto coefficient is used, but please feel to select any in [rdkit available](https://www.rdkit.org/docs/GettingStartedInPython.html#fingerprinting-and-molecular-similarity), e.g.: Tanimoto, Dice, Cosine, Sokal, Russel, Kulczynski, McConnaughey, and Tversky.

For now only single molecules can be queried, but it would be straight forward to implent multiple ligand queries.

The system was used (base) user@server:~/search_enamine_real_db$ uname -mrsv                                                                 Linux 4.15.0-45-generic #48-Ubuntu SMP Tue Jan 29 16:28:13 UTC 2019 x86_64  

## Setup
The functionality was tested on a custom linux server:
```
(base) user@server:~/search_enamine_real_db$ uname -mrsv                                                                 Linux 4.15.0-45-generic #48-Ubuntu SMP Tue Jan 29 16:28:13 UTC 2019 x86_64  
```
using `Intel(R) Xeon(R) Platinum 8180 CPU @ 2.50GHz` CPUs

### Install rdkit using anaconda
- install [rdkit](https://www.rdkit.org/docs/Install.html#cross-platform-under-anaconda-python-fastest-install) using anaconda

Optional:
- add [postgres](https://www.rdkit.org/docs/Install.html#installing-and-using-postgresql-and-the-rdkit-postgresql-cartridge-from-a-conda-environment) installation using anaconda

### Data
Download the [Enamine REAL](https://enamine.net/library-synthesis/real-compounds/real-database#) files (log in needed):
```bash
>>> wc -l 2019q1-2_Enamine_REAL_723M_SMILES_Part_01.smiles
60262532 2019q1-2_Enamine_REAL_723M_SMILES_Part_01.smiles
```

List of available Fingerprints in rdkit: [here]

## Usage Details
There are several keyword arguments you can pass to the python scripts, of which some can be further passed to 
the available shell scripts.

Arguments of `src/create_blobs.py`:
```cmd
user@server:~/search_enamine_real_db$ python create_blobs.py --help
usage: create_blobs.py [-h] [--reference_mol REFERENCE_MOL]
                       [--input_folder INPUT_FOLDER] [--pattern PATTERN]
                       [--tanimoto_threshold TANIMOTO_THRESHOLD] [-f FORCE]
                       [--outpath OUTPATH] [-n NUMBER]

optional arguments:
  -h, --help            show this help message and exit
  --reference_mol REFERENCE_MOL
                        Molecule(s) of interest for which analoges are
                        required.
  --input_folder INPUT_FOLDER
                        Search path for inputs. Subdirectories will be
                        included!
  --pattern PATTERN     pattern of files which should be
  --tanimoto_threshold TANIMOTO_THRESHOLD
                        Include virtual molecules with a similarity to
                        targets.
  -f FORCE, --force FORCE
                        Set flag to overwrite previous results
  --outpath OUTPATH     Path to save blobs of fingerprints.
  -n NUMBER, --number NUMBER
                        Select part of enamine real database. 1-12
```

None of the arguments is currently passed to the creation of blobs in `bin/create_blobs.sh`:
```cmd
user@server:~/search_enamine_real_db$ bash bin/create_blobs.sh --help
Creates binary files for Enamine project.
./create_blobs.sh --cpus 4 | [-h]]
--cpus: 2, 3, 4, 6 or 12
```


Arugments of `src/search_database_singleprocess.py`:
```
user@server:~/search_enamine_real_db$ python search_database_singleprocess.py --help
usage: search_database_singleprocess.py [-h] [--reference_mol REFERENCE_MOL]
                                        [--input_folder INPUT_FOLDER]
                                        [--pattern PATTERN]
                                        [--tanimoto_threshold TANIMOTO_THRESHOLD]
                                        [-f FORCE] [--outpath OUTPATH]

optional arguments:
  -h, --help            show this help message and exit
  --reference_mol REFERENCE_MOL
                        Molecule(s) of interest for which analoges are
                        required.
  --input_folder INPUT_FOLDER
                        Search path for inputs. Subdirectories will be
                        included!
  --pattern PATTERN     pattern of files which should be
  --tanimoto_threshold TANIMOTO_THRESHOLD
                        Include virtual molecules with a similarity to
                        targets.
  -f FORCE, --force FORCE
                        Set flag to overwrite previous results
  --outpath OUTPATH     Path to save results.
```

which are partially available in `bin/run_search.sh`:
```
user@server:~/search_enamine_real_db$ bash bin/run_search.sh --help
Search Enamine database using stored fingerprints in ./blobs/
./run_search.sh --cpus 4 --query [SMILES] --threshold 0.8 | [-h]]

--cpus                  # of CPUs: 1, 2, 3, 4, 6 or 12
--query                 Query molecule as SMILES
-t, --threshold         Similarity threshold to keep results, default 0.7
--force                 Overwrite previous results, default false

-h, --help              Display this help message.

```

## Excecution Time of on the fly calculation of fingerprints
It is also possible to execute the script calculating fingerprints on the fly.

```
user@server: ~/search_enamine_real_db$ python search_database_singleprocess.py 
```
This will take less than 4 minutes and does
- read 100000 lines
- compares to 1 default molecule 

Profiling of the script yields that the calculation of the fingerprint is the most time consuming operation.
```cmd
(user@server):~/search_enamine_real_db$ python search_database.py 

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
> Main bottleneck is the calculation of fingerprints which seems to be CPU intensive 

## Notes
### Generators
- Since Python 3.7 Generators are stop using `return` instead of `stopIteration`

### Moding and integers in bash
On the server an integer with a leading error throws an error if it is larger than seven:
```cmd
user@cosmos: echo $((07 % 12))                                                                                                                                                                      7           7
user@cosmos: echo $((08 % 12))                                                                                                                                                                      7          
-bash: 08: value too great for base (error token is "08")
```
- Wanted: looping with integers of two digits, [see] (https://stackoverflow.com/questions/5099119/adding-a-zero-to-single-digit-variable)


## ToDos
### Self-build solution
- check [fragmentation algorithm](https://www.rdkit.org/docs/source/rdkit.Chem.Fraggle.FraggleSim.html#fragmentation-algorithm) on results?
- Speed up building of fingerprints: Ideas
     -  [load chuncs](https://stackoverflow.com/questions/49752452/using-a-python-generator-to-process-large-text-files)

### Solution using chemfp
- [chemfp](https://chemfp.readthedocs.io/en/latest/using-api.html)

### Solution using rdkit functionality
1. Create Postgres Database [Cartridge](https://www.rdkit.org/docs/Cartridge.html)
2. Use [SimilarityScreener](https://www.rdkit.org/docs/source/rdkit.Chem.Fingerprints.SimilarityScreener.html#rdkit.Chem.Fingerprints.SimilarityScreener.SimilarityScreener)

