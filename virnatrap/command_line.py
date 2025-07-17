import sys
import os
from os.path import isdir, isfile
import argparse
import virnatrap
from .virnatrap import run_virna_pred

# CLI script for dumping seed reads before contig assembly
DESCRIPTION = (
    "Dump seed reads for external contig assembly from unmapped RNAseq FASTQ files"
)
# Determine package directory
mpath = virnatrap.__file__
PWD = os.path.dirname(mpath)


def virnatrap_predict():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--input", type=str, help="Input directory containing FASTQ files", required=True
    )
    parser.add_argument(
        "--output", type=str, help="Output directory for seed files", required=True
    )
    parser.add_argument(
        "--fastmode", type=int, choices=[0, 1],
        help="(Ignored) seed dumping does not perform contig assembly", default=0
    )
    parser.add_argument(
        "--multi_proc", type=int, choices=[0, 1],
        help="Run in parallel using multiprocessing", default=1
    )
    parser.add_argument(
        "--num_threads", type=int,
        help="Number of worker threads for multiprocessing", default=48
    )
    parser.add_argument(
        "--model_path", type=str,
        help="Path to TensorFlow model for scoring reads", required=False
    )
    parser.add_argument(
        "--threshold", type=float,
        help="Prediction score threshold for selecting seed reads", default=0.7
    )

    args = parser.parse_args()
    inpath = args.input
    outpath = args.output

    # Validate paths
    if not isdir(inpath):
        print(f"Input directory '{inpath}' not found.")
        sys.exit(1)
    if not isdir(outpath):
        print(f"Output directory '{outpath}' not found.")
        sys.exit(1)

    # Parse flags
    fastmode = bool(args.fastmode)
    multi_proc = bool(args.multi_proc)
    num_threads = args.num_threads
    threshold = args.threshold

    # Determine model path
    if args.model_path and isfile(args.model_path):
        model_path = args.model_path
    else:
        default_model = os.path.join(PWD, 'model', 'model_lr_0.005_pool_5_emb_25_l2_0.02_64.hdf5')
        model_path = default_model
        if args.model_path:
            print(f"Model '{args.model_path}' not found; using default model at {default_model}.")

    # Run the pipeline with threshold
    print(f"Reading FASTQ files from '{inpath}' with threshold {threshold}...")
    run_virna_pred(inpath, outpath, fastmode, multi_proc, model_path, num_threads, threshold)
    print("Done.")


if __name__ == '__main__':
    virnatrap_predict()
