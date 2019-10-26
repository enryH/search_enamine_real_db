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

LIMIT = os.environ.get("LIMIT", 200000)

parser = argparse.ArgumentParser()
add_arguments(parser)
parser.add_argument("--cpus", type=int, default=os.cpu_count()-3,
                    help='How many CPU to run on in parallel. '
                         'Default: All available minus 3')
parser.add_argument("--limit", type=int, default=LIMIT,
                    help='Maximum number of molecules to compare to. '
                         'Default is set by LIMIT environment variable or 200.000')
parser.add_argument("-f", "--force", type=bool, default=False,  
                    help="Set flag to overwrite previous results")
parser.add_argument("--fp_type", default='rdkit',
                    help="Choose fingerprint type: rdkit, morgan2")
args = parser.parse_args()

print("Start process using: ", args)   

starttime= datetime.now()

# Select Fingerprint type
args.fp_type = args.fp_type.lower()
if args.fp_type == 'morgan2':
    from rdkit.Chem import AllChem
    def morgen_at(radius):
        """
        Return unmodded fingerprint with int-hash
        for atom environments up to radius 2.
        """
        def morgan_at_radius(mol):
            return AllChem.GetMorganFingerprint(mol, radius)
        return morgan_at_radius
    Fingerprint = morgen_at(2)
elif args.fp_type == 'rdkit':
    Fingerprint = FingerprintMols.FingerprintMol
else:
    raise IOError(f"Fingerprint argument did not match any known kind: {args.fp_type}")
# Calculate fp for reference for later comparison
reference = Chem.MolFromSmiles(args.reference_mol)
fp_reference= Fingerprint(reference)
ref = fp_reference

# parse command-line arguments
_cwd, inpath, outpath, inputs = process_args(args)

filename=inputs[0]
print(args)

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
    try:
        smiles, idx, *_ = row.split("\t")
        mol = Chem.MolFromSmiles(smiles)
        mol = Fingerprint(mol) #this takes the bulk of time
        # tanimoto_sim = DataStructs.TanimotoSimilarity(reference, mol)
        tanimoto_sim = TanimotoSimilarity(ref, mol)
        if tanimoto_sim >= args.tanimoto_threshold:
            print(f"Added Molecule with id: {idx}")
            l_result = [str(x) for x in [smiles, idx, tanimoto_sim]]
            return '\t'.join(l_result)
        # print(f"Too low Tanimoto Similarity for mol with id {idx}:\t {tanimoto_sim}")
    except:
        print("Failed to read row:", row)

def generator_from_file_in_blocks(f, limit=10000):
    if (limit == -1) or (limit is None):
        for row in f:
            yield row
    else:
        assert (limit > 0) and isinstance(limit, int), "Please specify an positive integer."
        for i, row in enumerate(f):
            if i > limit:
                # raise StopIteration(f"Reached Limit of {limit}") # Python <= 3.6
                print("Reached specified limit")
                return
            else:
                yield row

if __name__=='__main__':
    from concurrent import futures
    from itertools import takewhile
    import sys
    
    from utils import DumpsResults
    dump = DumpsResults(folder=args.reference_mol, path=outpath, overwrite=args.force)
    
    def check_file(filename, ref=ref, limit=args.limit):
        executor = futures.ProcessPoolExecutor(max_workers=args.cpus)
        with open(filename) as f:
            gen = generator_from_file_in_blocks(f, limit=limit)
            results = executor.map(process_row, gen, chunksize=1, timeout=None)
        results_filtered = []
        for result in results:
            if result is not None:
                results_filtered.append(result)
            if len(results_filtered) > 1000:
                dump(results_filtered)
                results_filtered = []
        return results_filtered
    
    results = check_file(filename)




    # data = [x.split("\t") for x in results]
    # import pandas as pd
    # df = pd.DataFrame(data, columns=["SMILES", "ID", "TanimotoSim"])
    # df = df.set_index("ID")
    # df.sort_values("TanimotoSim", ascending=False, inplace=True)
    # df.head()
    
    # # profile.run('print(check_file(filename)); print()')


    endtime = datetime.now()
    print(format_duration(starttime, endtime))

