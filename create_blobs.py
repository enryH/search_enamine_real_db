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

from tqdm import tqdm

import utils
from utils import add_arguments
from utils import process_args
from utils import format_duration
from utils import BlobsIO


N_PER_FILE = 60262532


parser = argparse.ArgumentParser()
add_arguments(parser)
parser.add_argument("--outpath", type=str, default='./blobs',
                    help="Path to save blobs of fingerprints.")
parser.add_argument('-n', '--number', type=int, default=2,
                    help='Select part of enamine real database. 1-12')



args = parser.parse_args()

starttime= datetime.now()

# AllChem.GetMorganFingerprint(reference,2)
# FingerprintMols.FingerprintMol(reference)
# FingerprintMols.GetRDKFingerprint

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


Fingerprint = FingerprintMols.FingerprintMol

_cwd, inpath, outpath, inputs = process_args(args)

assert args.number in range(1,13), "Select a number from 1 to 12 for flag --number."
#This will pick up the last filename, if it exits, but even it is not matched
for filename in inputs:
    # if fnmatch.fnmatch(filename, f'{args.number}.smiles'):
    if f'{args.number}.smiles' in filename:
        print("Selected:", filename)
        break
print("Selected file:", filename)

part_id = utils.find_int(filename)


i=1
cache = []
blobs = []



with open(filename) as f:
    blob_io = BlobsIO(part=part_id, path=outpath,
                      overwrite=True if args.force else False)
    for row in tqdm(f, desc=f'Part {part_id:02}', ascii=True, total=N_PER_FILE):
        t_result = process_row(row)
        if t_result is not None:
            cache.append(t_result)
        if (i % 1000001) == 0:
            fname = blob_io.write_blob(cache=cache)
            blobs.append(fname)
            cache = []
        i += 1
    fname = blob_io.write_blob(cache=cache)
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
