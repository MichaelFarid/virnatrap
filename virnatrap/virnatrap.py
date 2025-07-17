"""
viRNAtrap Seed Dumper

This script processes unmapped RNAseq FASTQ files, identifies viral reads using a trained model,
extracts the first 48 bp of reads above a score threshold, and writes them to a *_seeds.txt file.
"""

import os
import re
import sys
import glob
import random
import argparse
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import load_model
from tensorflow.keras.preprocessing.sequence import pad_sequences

# Suppress TF warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
from tensorflow import get_logger, autograph
get_logger().setLevel('ERROR')
autograph.set_verbosity(0)

# Constants
SEGMENT_LENGTH = 48
SEARCHSUBLEN   = 24
DEFAULT_NUC_ORDER = {y: x for x, y in enumerate(["A", "T", "C", "G"])}
NUCLEOTIDES = list(DEFAULT_NUC_ORDER.keys())

# Utility

def random_base(_=None):
    return random.choice(NUCLEOTIDES)


def handle_non_ATGC(seq):
    return re.sub('[^ATCG]', random_base, seq)


def encode_sequences(seqs):
    arr = []
    for s in seqs:
        s = s.upper()[:SEGMENT_LENGTH]
        encoded = [DEFAULT_NUC_ORDER.get(b,0) + 1 for b in s]
        arr.append(encoded)
    return np.array(pad_sequences(arr, maxlen=SEGMENT_LENGTH, padding='post'))

# FASTQ processing

def filter_sequences(seqs):
    bad = ['A'*SEARCHSUBLEN, 'C'*SEARCHSUBLEN, 'G'*SEARCHSUBLEN, 'T'*SEARCHSUBLEN]
    return [s for s in seqs if not any(b in s for b in bad)]


def proc_fastq(fastq_file):
    with open(fastq_file) as f:
        lines = f.readlines()
    reads = []
    for i in range(1, len(lines), 4):
        seq = handle_non_ATGC(lines[i].strip())
        reads.append(seq)
    unique = list(dict.fromkeys(reads))
    filtered = filter_sequences(unique)
    medlen = int(np.median([len(s) for s in filtered]))
    trimmed = [s[:medlen] for s in filtered]
    return encode_sequences(trimmed), trimmed

# Main

def dump_seeds(fastq_dir, out_dir, model_path, threshold):
    # Load model
    model = load_model(model_path)

    # Ensure output dir exists
    os.makedirs(out_dir, exist_ok=True)

    # Process all .fastq files
    for fq in glob.glob(os.path.join(fastq_dir, '*.fastq')):
        base = os.path.basename(fq).replace('.fastq','')
        seed_file = os.path.join(out_dir, f"{base}_seeds.txt")
        # Skip if already exists
        if os.path.exists(seed_file):
            continue

        # Encode and score
        X, seqs = proc_fastq(fq)
        scores = model.predict(X).flatten()

        # Select seeds
        seeds = []
        for s, sc in zip(seqs, scores):
            if sc > threshold and len(s) >= SEGMENT_LENGTH:
                seeds.append((s[:SEGMENT_LENGTH], sc))

        # Write seeds
        with open(seed_file, 'w') as outf:
            for idx, (seq, sc) in enumerate(seeds):
                outf.write(f">seed{idx}[{sc:.3f}]\n{seq}\n")

        print(f"Wrote {len(seeds)} seeds to {seed_file}")

# CLI

if __name__ == '__main__':
    p = argparse.ArgumentParser(description="Dump viral 48 bp seeds from unmapped FASTQ")
    p.add_argument('--input',  '-i', required=True, help='Directory of unmapped FASTQ files')
    p.add_argument('--output', '-o', required=True, help='Directory to write *_seeds.txt')
    p.add_argument('--model',  '-m', required=True, help='Path to Keras model (.h5)')
    p.add_argument('--threshold', '-t', type=float, default=0.7, help='Score threshold for seeds')
    args = p.parse_args()

    dump_seeds(args.input, args.output, args.model, args.threshold)
