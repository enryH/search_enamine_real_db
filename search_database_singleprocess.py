import os
import argparse
import profile

from datetime import datetime

import rdkit
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

from utils import add_arguments
from utils import process_args
from utils import format_duration

parser = argparse.ArgumentParser()
add_arguments(parser)
args = parser.parse_args()

starttime= datetime.now()

# AllChem.GetMorganFingerprint(reference,2)
# FingerprintMols.FingerprintMol(reference)
# FingerprintMols.GetRDKFingerprint

Fingerprint = FingerprintMols.FingerprintMol

reference = Chem.MolFromSmiles(args.reference_mol)
fp_reference= Fingerprint(reference)
ref = fp_reference

print(f"Fingerprint lenght: {ref.GetNumBits()}")

_cwd, inpath, outpath, inputs = process_args(args)

filename=inputs[0]

def process_row(row:str, ref=ref):
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
    smiles, id, *_ = row.split("\t")
    mol = Chem.MolFromSmiles(smiles)
    mol = Fingerprint(mol) #this takes the bulk of time
    # tanimoto_sim = DataStructs.TanimotoSimilarity(reference, mol)
    tanimoto_sim = TanimotoSimilarity(ref, mol)
    if tanimoto_sim >= args.tanimoto_threshold:
        print(f"Added Molecule to results from line: {line}")
        return '\t'.join(smiles, id, tanimoto_sim)

if __name__=='__main__':
    LIMIT = 1000
    def check_file(filename, ref=ref, limit=LIMIT):
        # global results
        results = []
        with open(filename) as f:
            for line_no, row in enumerate(f):
                try:
                    result = process_row(row)
                    if result is not None:
                        results.append(result)
                except:
                    print(f"Failed to read line {line_no}")
                if line_no >= limit:
                    break
        return results
    # check_file(filename)
    profile.run('print(check_file(filename)); print()')


    endtime = datetime.now()
    print(format_duration(starttime, endtime))

