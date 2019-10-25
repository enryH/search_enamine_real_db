import os
import glob
import pickle

def add_arguments(parser):
    parser.add_argument("--reference_mol", type=str, default='COC1=C(OCCCN2CCOCC2)C=C2C(NC3=CC(Cl)=C(F)C=C3)=NC=NC2=C1',  # Gefitinib, Iressa
                        help='Molecule(s) of interest for which analoges are required.')
    parser.add_argument("--outpath", type=str, default='./results',
                        help="Path to save results.")
    parser.add_argument('--input_folder',  type=str, default='.', 
                        help='Search path for inputs. Subdirectories will be included!')
    parser.add_argument('--pattern', type=str, default='.smiles',
                        help='pattern of files which should be ')
    parser.add_argument('--tanimoto_threshold', type=float, 
                        default=0.5,
                        help=('Include virtual molecules with a similarity to targets.')
                        )

def process_args(args):
    """
    Process some command-line arguments:

    inpath, outpath: normalize to current os
    pattern: search on inpath for files with specified pattern

    """
    
    _cwd = glob.os.getcwd()
    print(f"Script is executed from: {_cwd}")

    inpath=glob.os.path.abspath(args.input_folder)
    outpath=glob.os.path.abspath(args.outpath)

    try:
        glob.os.chdir(inpath)
    except:
        print(f"Could change to {inpath}")
        print(f"Search for files in: {_cwd}")
        inpath=_cwd

    inputs = glob.glob('*/*'+args.pattern, recursive=True)

    if len(inputs) > 0:
        print(f"Found {len(inputs)} input files:")
        print("-",'\n- '.join(inputs))
    else:
        raise FileNotFoundError(f"No input files found using pattern: {args.pattern} "
                                f"in {_cwd} folder and all its subfolders.")
    return _cwd, inpath, outpath, inputs


def format_duration(starttime, endtime):
    seconds = int((endtime-starttime).total_seconds())
    hours, seconds   = divmod(seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    return "Time elapsed {:02}:{:02}:{:02}".format(hours, minutes, seconds)



class DumpsResults():
    """
    A call pickles results to a specified location by its
    folder and name. Repeated calls dump results to consecutely
    numbered files.
    """
    def __init__(self, folder:str, path='.', fname="results", overwrite=False):
        self.parentfolder = path
        self.folder = folder
        self.fname = fname + "_{}"
        self.path = os.path.join(self.parentfolder, self.folder)
        self.fpath = os.path.join(self.path, fname)
        try:
            os.makedirs(self.path)
        except FileExistsError:
            if overwrite:
                pass
            else:
                raise FileExistsError("Run previously executed. Abort")
        self.i = 1

    def __call__(self, obj):
        fname = self.fpath.format(self.i)
        with open(fname, "wb") as f:
            pickle.dump(obj, file=f)
            print("Saved results to:", fname)
        self.i += 1
