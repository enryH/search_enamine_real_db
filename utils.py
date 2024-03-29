import os
import glob
import pickle

def add_arguments(parser):
    parser.add_argument("--reference_mol", type=str, default='COC1=C(OCCCN2CCOCC2)C=C2C(NC3=CC(Cl)=C(F)C=C3)=NC=NC2=C1',  # Gefitinib, Iressa
                        help='Molecule(s) of interest for which analoges are required.')
    parser.add_argument('--input_folder',  type=str, default='.', 
                        help='Search path for inputs. Subdirectories will be included!')
    parser.add_argument('--pattern', type=str, default='.smiles',
                        help='pattern of files which should be ')
    parser.add_argument('--tanimoto_threshold', type=float, 
                        default=0.7,
                        help=('Include virtual molecules with a similarity to targets.')
                        )
    parser.add_argument("-f", "--force", type=bool, default=False,  
                    help="Set flag to overwrite previous results")

def get_files(inpath, pattern):
    try:
        os.chdir(inpath)
    except FileNotFoundError as err:
        print(f"Could not find folder: {inpath}")
        raise FileNotFoundError("Could not find processed files:\n"
                    f"in relativ path: {inpath}\n"
                    f"in working Dir: {os.getcwd()}\n"
                    f"giving absolut path: {os.path.normpath(inpath)}\n")

    inputs = glob.glob('**/*'+pattern, recursive=True)
    return inputs

def process_args(args):
    """
    Process some command-line arguments:

    inpath, outpath: normalize to current os
    pattern: search on inpath for files with specified pattern

    """
    
    _cwd = glob.os.getcwd()
    print(f"Script is executed from: {_cwd}")

    inpath=glob.os.path.abspath(args.input_folder)
    # inpath=os.path.normpath(inpath)

    outpath=glob.os.path.abspath(args.outpath)

    inputs = get_files(inpath, args.pattern)    

    if len(inputs) > 0:
        inputs = [os.path.abspath(f) for f in inputs]
        print(f"Found {len(inputs)} input files:")
        print("-",'\n- '.join(inputs))
    else:
        raise FileNotFoundError(f"No input files found using pattern: {args.pattern} "
                                f"in {inpath} folder and all its subfolders.")
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
        self.fname = fname + "_{:02}"
        self.path = os.path.join(self.parentfolder, self.folder)
        self.fpath = os.path.join(self.path, self.fname)
        
        os.makedirs(self.path, exist_ok=True) # folder is allowed to exist
        self.i = 1
        if self.fname.format(self.i) in os.listdir(self.path):  
            if overwrite:
                pass
            else:
                raise FileExistsError("Run previously executed. Abort")
        

    def __call__(self, obj):
        fname = self.fpath.format(self.i)
        with open(fname, "w") as f:
            f.writelines(obj)
            print("Saved results to:", fname)
        self.i += 1


class BlobsIO():
    def __init__(self, part, path='./blobs', overwrite=False):
        self.path = path
        folder = 'part_{:02}'.format(part)
        self.folder = os.path.join(self.path, folder)
        try:
            os.makedirs(self.folder, exist_ok=False)
        except IOError:
            if overwrite:
                pass
            else:
                raise IOError('Path already exist.'
                              f'Please delete before continuing: {self.folder}')
        self.k = 0

    def write_blob(self, cache):
        fname = os.path.join(self.folder, f"real_blob_{self.k:02}.pkl")
        with open(fname, 'wb') as f:
            pickle.dump(cache, f)
        print("Wrote file:", fname)
        self.k += 1
        return fname

    def get_list_of_blobs(self, part=None):
        raise NotImplementedError

def read_blob(file):
    with open(file, 'rb') as f:
        l_in_mem = pickle.load(f)
    return l_in_mem

def make_string(*args):
    l_result = [str(x) for x in args]
    return '\t'.join(l_result) + '\n'

def find_int(filename):
    "Return last integer in fname."
    return int(filename.split('.')[0].split('_')[-1])

