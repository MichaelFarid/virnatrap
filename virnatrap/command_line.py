import sys
import os
from os.path import isdir, isfile
from virnatrap.virnatrap import run_virna_pred
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
        "--input", "-i",
        type=str,
        help="input directory containing *.fastq files",
        required=True
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="output directory for *_contigs.txt and *_seeds.txt",
        required=True
    )
    parser.add_argument(
        "--fastmode", "-f",
        type=int,
        choices=[0,1],
        default=0,
        help="1 to use C‐accelerated assembler, 0 to use pure‐Python fallback"
    )
    parser.add_argument(
        "--multi_proc", "-m",
        type=int,
        choices=[0,1],
        default=1,
        help="1 to parallelize across FASTQ files with multiprocessing, 0 to run serially"
    )
    parser.add_argument(
        "--num_threads", "-t",
        type=int,
        default=48,
        help="number of worker processes if --multi_proc=1"
    )
    parser.add_argument(
        "--model_path",
        type=str,
        help="path to your Keras/TensorFlow model (.h5)",
        required=False
    )
    parser.add_argument(
        "--dump_seeds", 
        action="store_true",
        help="just write the 48 bp seeds above threshold to *_seeds.txt and skip contig assembly"
    )

    args = parser.parse_args()

    inpath     = args.input
    outpath    = args.output
    fastmode   = bool(args.fastmode)
    multi_proc = bool(args.multi_proc)
    num_threads= args.num_threads
    dump_seeds = args.dump_seeds

    # Validate directories
    if not isdir(inpath):
        print(f"ERROR: input directory {inpath} not found", file=sys.stderr)
        sys.exit(1)
    if not isdir(outpath):
        print(f"ERROR: output directory {outpath} not found", file=sys.stderr)
        sys.exit(1)

    # Model location
    if args.model_path and isfile(args.model_path):
        model_path = args.model_path
    else:
        if args.model_path:
            print(f"WARNING: model {args.model_path} not found—using default", file=sys.stderr)
        model_path = os.path.join(PWD, 'model',
            'model_lr_0.005_pool_5_emb_25_l2_0.02_64.hdf5'
        )

    print(f"Reading FASTQ files from {inpath}…", file=sys.stderr)

    # Call the main driver
    run_virna_pred(
        inpath,
        outpath,
        fastmode,
        multi_proc,
        model_path,
        num_threads,
        dump_seeds
    )

    print("All done.", file=sys.stderr)


if __name__ == '__main__':
    virnatrap_predict()
