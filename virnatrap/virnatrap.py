"""
Model classes and functions to identify viral reads and assemble viral contigs from input fastq
"""
# Imports --------------------------------------------------------------------------------------------------------------
import os
import re
import sys
from multiprocessing import Pool, freeze_support
import multiprocessing as mp
import numpy as np
from pkg_resources import resource_filename
from collections import OrderedDict, defaultdict
import random
import tensorflow as tf
import ctypes
from ctypes import *
from tensorflow import get_logger, autograph
from tensorflow import keras
from tensorflow.keras import Sequential
from tensorflow.keras.layers import *
from tensorflow.keras.models import load_model
from tensorflow.keras.preprocessing.sequence import pad_sequences
import glob
import argparse
from os.path import exists

DESCRIPTION = (
    "Extract viral contigs from a directory with unmapped RNAseq reads fastq files "
    "and saves a file with contigs for each fastq in an output directory"
)

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
get_logger().setLevel('ERROR')
autograph.set_verbosity(0)

# Constants ------------------------------------------------------------------------------------------------------------
DEFAULT_NUC_ORDER = {y: x for x, y in enumerate(["A", "T", "C", "G"])}
NUCLEOTIDES = sorted(DEFAULT_NUC_ORDER.keys())
SEGMENT_LENGTH = 48
SEARCHSUBLEN = 24

PWD = os.getcwd()

# Utility functions ---------------------------------------------------------------------------------------------------
def random_base(seq=None):
    return random.choice(NUCLEOTIDES)

def handle_non_ATGC(sequence):
    ret = re.sub('[^ATCG]', random_base, sequence)
    assert len(ret) == len(sequence)
    return ret

def encode_sequence(sequence, nuc_order=None):
    if nuc_order is None:
        nuc_order = DEFAULT_NUC_ORDER
    seq = sequence.upper()[:SEGMENT_LENGTH]
    assert re.match(f'^[{"".join(nuc_order.keys())}]+$', seq)
    encoded = [nuc_order[b] + 1 for b in seq]
    return np.array([encoded])

def encode_sequences(seqs, nuc_order=None, segment_length=SEGMENT_LENGTH):
    arr = [encode_sequence(s, nuc_order)[0] for s in seqs]
    return np.array(pad_sequences(arr, maxlen=segment_length, padding='post'))

# Model utils ---------------------------------------------------------------------------------------------------------
def save_model(model, model_path):
    model.save(model_path + '.h5')

def load_model_keras(model_path):
    return load_model(model_path)

# C-accelerated assembly wrapper --------------------------------------------------------------------------------------
def make_clist(lst):
    return (c_char_p * len(lst))(*[x.encode() for x in lst])

def assemble_read_call_c(reads_seed, reads_all, scores_all, seed_scores, out_file):
    so = CDLL(os.path.join(PWD, 'src', 'assemble_read_c.so'))
    arr_all = (ctypes.c_float * len(scores_all))(*scores_all)
    arr_seed = (ctypes.c_float * len(seed_scores))(*seed_scores)
    out_c = c_char_p(out_file.encode())
    m = c_int(len(reads_all))
    n = c_int(len(reads_seed[0]))
    clist_all = make_clist(reads_all)
    clist_seed = make_clist(reads_seed)
    so.assemble_read_loop.restype = ctypes.c_int
    return so.assemble_read_loop(arr_all, arr_seed, clist_all, clist_seed,
                                  n, m, len(seed_scores), out_c)

# Pure-Python assembly (fallback) --------------------------------------------------------------------------------------
def assemble_right(read, read_list, score_list, score_read):
    # ... existing code ...
    pass

def assemble_left(read, read_list, score_list, score_read):
    # ... existing code ...
    pass

def assemble_read(read, read_list, score_list, score_read):
    # ... existing code ...
    pass

def assemble_read_loop(reads_seed, reads_all, scores_all, seed_scores, out_file):
    # ... existing code ...
    pass

# Main extraction -----------------------------------------------------------------------------------------------
def extract_contigs(invars, large_file_thr=1000000, dump_seeds=False):
    inpath, outpath, fastmode, model_path = invars
    base = os.path.basename(inpath).replace('_unmapped.fastq', '')
    fn_contigs = os.path.join(outpath, f"{base}_contigs.txt")
    fn_seeds = os.path.join(outpath, f"{base}_seeds.txt")
    if exists(fn_contigs):
        return 0

    model = load_model_keras(model_path)
    encoded, seqs = proc_fastq(inpath)
    scores = model.predict(encoded).flatten()

    # select seeds
    thr = 0.7
    seeds = [s[:SEGMENT_LENGTH] for s, sc in zip(seqs, scores) if sc > thr and len(s)>=SEGMENT_LENGTH]
    seed_scores = [sc for sc in scores if sc > thr]

    if dump_seeds:
        with open(fn_seeds, 'w') as f:
            for i, (s, sc) in enumerate(zip(seeds, seed_scores)):
                f.write(f">seed{i}[{sc:.3f}]\n{s}\n")

    if fastmode:
        assemble_read_call_c(seeds, seqs, scores.tolist(), seed_scores, fn_contigs)
    else:
        assemble_read_loop(seeds, seqs, scores.tolist(), seed_scores, fn_contigs)

    return 1

# Fastq processing -----------------------------------------------------------------------------------------------
def filter_sequences(seqs):
    bad = ['A'*SEARCHSUBLEN, 'C'*SEARCHSUBLEN, 'G'*SEARCHSUBLEN, 'T'*SEARCHSUBLEN]
    return [s for s in seqs if all(b not in s for b in bad)]

def proc_fastq(infile):
    with open(infile) as f:
        lines = f.readlines()
    seqs = [handle_non_ATGC(lines[i].strip()) for i in range(1, len(lines), 4)]
    seqs = list({handle_non_ATGC(s) for s in seqs})
    seqs = filter_sequences(seqs)
    med = int(np.median([len(s) for s in seqs]))
    seqs = [s[:med] for s in seqs]
    return encode_sequences(seqs), seqs

from functools import partial

# Updated run_virna_pred to avoid pickling lambda
def run_virna_pred(inpath, outpath, fastmode, multi_proc, model_path, num_threads, dump_seeds=False):
    print("starting now")
    # Gather input and output files
    infs = list(set([f for f in glob.glob(inpath + '/*.fastq')]))
    outs = glob.glob(outpath + '/*.txt')
    bases_in = {os.path.basename(f).replace('_unmapped.fastq','') for f in infs}
    bases_out = {os.path.basename(f).replace('_contigs.txt','') for f in outs}
    to_process = [f for f in infs if os.path.basename(f).replace('_unmapped.fastq','') not in bases_out]

    args = [[f, outpath, fastmode, model_path] for f in to_process]

    print('starting_prediction...')
    if multi_proc:
        freeze_support()
        pool = mp.Pool(processes=num_threads)
        # Use partial instead of lambda for pickleability
        func = partial(extract_contigs, dump_seeds=dump_seeds)
        pool.map(func, args)
        pool.close()
        pool.join()
    else:
        for iv in args:
            extract_contigs(iv, dump_seeds=dump_seeds)

    print("Done processing")
