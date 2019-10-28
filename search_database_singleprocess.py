"""
    ipython search_database_singleprocess.py -i -- --input_folder blobs --pattern pkl
"""


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
from utils import read_blob
from utils import make_string
from utils import DumpsResults

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


# Not cached:
def process_row_str(row:str, ref=ref):
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
    fp = Fingerprint(mol) #this takes the bulk of time
    # tanimoto_sim = DataStructs.TanimotoSimilarity(reference, mol)
    tanimoto_sim = TanimotoSimilarity(ref, fp)
    if tanimoto_sim >= args.tanimoto_threshold:
        print(f"Added Molecule to results from line: {line}")
        return make_string(smiles, id, tanimoto_sim)


if __name__=='__main__':
    def check_file(filename, ref=ref, limit=100000):
        # global results
        results = []

        file_type= filename.split(".")[-1]
        if file_type == "smiles":
            with open(filename) as f:
                for line_no, row in enumerate(f):
                    try:
                        result = process_row_str(row)
                        if result is not None:
                            results.append(result)
                            # if len(results) > 1000:
                            #     dump(results)
                            #     results = []
                    except:
                        print(f"Failed to read line {line_no}")
                    if line_no >= limit:
                        break
            return results
        elif file_type == "pkl":
            fp_in_mem = read_blob(filename)
            blob_id = int(filename.split('.')[0].split('_')[-1])
            fname_result = 'results_blob_{:03}_no'.format(blob_id)
            dump = DumpsResults(folder=args.reference_mol, path=outpath, 
                                fname=fname_result, overwrite=args.force)
            for _tuple in fp_in_mem:
                try: 
                    idx, fp, smiles = _tuple
                    tanimoto_sim = TanimotoSimilarity(ref, fp)
                    if tanimoto_sim >= args.tanimoto_threshold:
                        print(f"Added Molecule to results: {idx}")
                        result = make_string(smiles, idx, tanimoto_sim)
                        results.append(result)
                        if len(results) > 1000:
                            dump(results)
                            results = []
                except TypeError as e:
                    print(e)
            dump(results)
            print(f"Checked for {len(fp_in_mem)} molecules.")
        return 

    dump = DumpsResults(folder=args.reference_mol, path=outpath, overwrite=args.force)
    for filename in inputs:
        check_file(filename)

    endtime = datetime.now()
    print(format_duration(starttime, endtime))

