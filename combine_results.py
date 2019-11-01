import os
import glob
import pickle
import argparse

from utils import get_files

def combine_results(folder, threshold=None):
    print(os.getcwd())
    files = get_files(os.path.normpath(folder), 'results_blob*')
    assert len(files)>0, f"No files found in {folder}"
    print("Found Files: ")
    print('\n'.join(files))
    counter = 0
    if threshold is not None:
        assert isinstance(threshold,float), f"Threshold has to be float not, {type(threshold)} containing {threshold}"
    fname = f'results_threshold_{threshold}.txt'
    with open(fname, 'w') as out:
        for _file in files:
            print("Read:",_file)
            with open(_file, 'r') as f:
                l_result_str = f.readlines()
            counter += len(l_result_str)
            out.writelines(l_result_str)
    return {'N': counter, 'fname': fname}
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--resultpath", type=str, 
                        default='./results',
                        help="Path to save results.")
    parser.add_argument("--threshold", type=float,
                        default=None)
    args = parser.parse_args()
    info = combine_results(args.resultpath, threshold=args.threshold)
    from pprint import pprint
    pprint(info)
