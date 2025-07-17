import sys
import os
from os.path import isdir, isfile
from .virnatrap import run_virna_pred
import argparse
import virnatrap

# Constants ------------------------------------------------------------------------------------------------------------
DESCRIPTION = (
    "Extract viral contigs from a directory with unmapped RNAseq reads fastq files "
    "and saves a file with contigs for each fastq in an output directory"
)
mpath = virnatrap.__file__
PWD = mpath[:-12]

# Terminal functions ---------------------------------------------------------------------------------------------------
def virnatrap_predict():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--input", type=str, help="input directory", required=True
    )
    parser.add_argument(
        "--output", type=str, help="output directory", required=True
    )
    parser.add_argument(
        "--fastmode", type=int,
        help="1 to run in fastmode and call C function to assemble reads",
        required=False,
        default=0,
    )
    parser.add_argument(
        "--multi_proc", type=int,
        help="run with pool multi processing, if many input files",
        required=False,
        default=1,
    )
    parser.add_argument(
        "--num_threads", type=int,
        help="number of threads to run with pool multi processing",
        required=False,
        default=48,
    )
    parser.add_argument(
        "--model_path", type=str,
        help="path to Tensorflow model to predict whether reads come from viruses",
        required=False,
    )
    parser.add_argument(
        "--dump_seeds", action="store_true",
        help="write raw 48bp seeds to a separate *_seeds.txt file",
    )

    args = parser.parse_args()

    inpath = args.input
    outpath = args.output
    fastmode = bool(args.fastmode)
    multi_proc = bool(args.multi_proc)
    num_threads = args.num_threads
    dump_seeds = args.dump_seeds

    if not isdir(inpath):
        print(f"input directory {inpath} not found")
        sys.exit(1)
    if not isdir(outpath):
        print(f"output directory {outpath} not found")
        sys.exit(1)

    if args.model_path and isfile(args.model_path):
        model_path = args.model_path
    else:
        if args.model_path:
            print(f"model {args.model_path} not found, using default viRNAtrap model")
        model_path = os.path.join(PWD, 'model', 'model_lr_0.005_pool_5_emb_25_l2_0.02_64.hdf5')

    print(f"Reading fastq at {inpath}...")
    run_virna_pred(
        inpath,
        outpath,
        fastmode,
        multi_proc,
        model_path,
        num_threads,
        dump_seeds
    )
    print("done")

if __name__ == '__main__':
    virnatrap_predict()
