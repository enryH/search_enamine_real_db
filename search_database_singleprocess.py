"""
Screening molecules against REAL Enamine database. 
Giving `.smiles`-files as input will read the file(s) in the input-folder line by line,
compute the SMILES fingerprint representation and then compare to the reference.
Specifying as pattern pickeld files, ending on `.pkl`, will read in preprocessed files.
The current fingerprint for comparing molecules is the RDKit-Fingerprint (Daylight). 

python search_database_singleprocess.py --input_folder blobs --pattern pkl --reference_mol CCC1CC
"""
import os
import argparse

from datetime import datetime

import pathlib

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
from utils import find_int

parser = argparse.ArgumentParser()
add_arguments(parser)
parser.add_argument("--outpath", type=str, default='./results',
                    help="Path to save results.")
args = parser.parse_args()

starttime= datetime.now()

# AllChem.GetMorganFingerprint(reference,2)
# FingerprintMols.FingerprintMol(reference)
# FingerprintMols.GetRDKFingerprint

Fingerprint = FingerprintMols.FingerprintMol


reference = Chem.MolFromSmiles(args.reference_mol)
fp_reference= Fingerprint(reference)
ref = fp_reference

print(f"Fingerprint length: {ref.GetNumBits()}")

_cwd, inpath, outpath, inputs = process_args(args)

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

def check_for_part_in_path(path:pathlib.Path):
    folder = path.parent.stem
    if "part" in folder:
        return folder
    else:
        return ''

def check_file(filename, ref=ref, limit=100000):
    results = []
    path = pathlib.Path(filename)
    file_type=path.suffix
    if file_type == ".smiles":
        with open(filename) as f:
            for line_no, row in enumerate(f):
                try:
                    result = process_row_str(row)
                    if result is not None:
                        results.append(result)
                        # # Save results every N found molecules
                        if len(results) > 1000:
                            dump(results)
                            results = []
                except:
                    print(f"Failed to read line {line_no}")
                if line_no >= limit:
                    break
        return results
    elif file_type == ".pkl":
        id = find_int(filename)
        fname_result = 'results_blob_{:02}'.format(id)
        subfolder = check_for_part_in_path(path)
        try:
            #abort in case file already exists.
            dump = DumpsResults(folder=os.path.join(args.reference_mol, subfolder), 
                            path=outpath, 
                            fname=fname_result, overwrite=args.force)           
            fp_in_mem = read_blob(filename)
            for _tuple in fp_in_mem:
                try: 
                    idx, fp, smiles = _tuple
                    tanimoto_sim = TanimotoSimilarity(ref, fp)
                    if tanimoto_sim >= args.tanimoto_threshold:
                        result = make_string(smiles, idx, tanimoto_sim)
                        results.append(result)
                        # # Save results every N found molecules
                        if len(results) >= 1000:
                            dump(results)
                            results = []
                except TypeError as e:
                    print(e)
            dump(results)
            print(f"Checked for {len(fp_in_mem)} molecules.")
        except FileExistsError:
            print(f"Files for blob already exist: f{filename} ")
    else:
        raise ValueError("Filetype unkown: {}".format(file_type))
    return 


for filename in inputs:
    print("Check fps in blob:", filename)
    check_file(filename)

endtime = datetime.now()
print(format_duration(starttime, endtime))
