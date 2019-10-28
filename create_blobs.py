import os
import pickle
import argparse
import profile
import fnmatch

from datetime import datetime

import rdkit
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

from utils import add_arguments
from utils import process_args
from utils import format_duration
from utils import write_blob


parser = argparse.ArgumentParser()
add_arguments(parser)
parser.add_argument('-n', '--number', type=int, default=2,
                    help='Select part of enamine real database. 1-12')
args = parser.parse_args()

starttime= datetime.now()

# AllChem.GetMorganFingerprint(reference,2)
# FingerprintMols.FingerprintMol(reference)
# FingerprintMols.GetRDKFingerprint

Fingerprint = FingerprintMols.FingerprintMol

_cwd, inpath, outpath, inputs = process_args(args)

assert args.number in range(1,13), "Select a number from 1 to 12 for flag --number."
#This will pick up the last filename, if it exits, but even it is not matched
for filename in inputs:
    if fnmatch.fnmatch(filename, '{}.smiles'):
        break
print("Selected file:", filename)

i=1
k=1
cache = []
blobs = []


def process_row(row:str):
    """
    Inputs:
        row: str
        First two tab-seperated entries of string are the 
        SMILES and ID of the row.
    Return:
        str: 
        If the tanimoto threshold is passed, return the 
        smile, the id and the tanimoto coef. in smiles-format:
        "SMILES\tID\tTANIMOTOCOEF"
    """
    try:
        smiles, idx, *_ = row.split("\t")
        mol = Chem.MolFromSmiles(smiles)
        fp = Fingerprint(mol) #this takes the bulk of time
        return (idx, fp, smiles)
        # print(f"Too low Tanimoto Similarity for mol with id {idx}:\t {tanimoto_sim}")
    except:
        print("Failed to read row:", row)

with open(filename) as f:
    
    for line_no, row in enumerate(f):
        t_result = process_row(row)
        if t_result is not None:
            cache.append(t_result)
        if (i % 1000001) == 0:
            fname = write_blob(k, cache)
            blobs.append(fname)
            cache = []
            k += 1
        i += 1
    fname = write_blob(k, cache)
    blobs.append(fname)
    cache = []

endtime = datetime.now()
print(format_duration(starttime, endtime))


# fname = 
# gen = read_blob(fname)
# for r in gen: print(r)



# ref = 'COC1=C(OCCCN2CCOCC2)C=C2C(NC3=CC(Cl)=C(F)C=C3)=NC=NC2=C1'
# ref = Chem.MolFromSmiles(ref)
# fp = Fingerprint(ref)

# fp_explicit = rdkit.DataStructs.cDataStructs.ExplicitBitVect
