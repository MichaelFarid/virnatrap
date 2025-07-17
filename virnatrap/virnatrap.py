"""
Model classes and functions to identify viral reads and dump seed reads before contig assembly
"""
# Imports --------------------------------------------------------------------------------------------------------------
import os
import re
import sys
from multiprocessing import Pool, freeze_support
import multiprocessing as mp
import numpy as np
from pkg_resources import resource_filename
from collections import OrderedDict
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

# Suppress TensorFlow logs
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
get_logger().setLevel('ERROR')
autograph.set_verbosity(0)

# Globals -------------------------------------------------------------------------------------------------------------
DEFAULT_NUC_ORDER = {y: x for x, y in enumerate(["A", "T", "C", "G"])}
NUCLEOTIDES = sorted([x for x in DEFAULT_NUC_ORDER.keys()])
SEGMENT_LENGTH = 48
SEARCHSUBLEN = 24
PWD = os.getcwd()

# Utility functions ---------------------------------------------------------------------------------------------------
flatten = lambda l: [item for sublist in l for item in sublist]

def random_base(seq=None):
    return random.choice(NUCLEOTIDES)

def handle_non_ATGC(sequence):
    ret = re.sub('[^ATCG]', random_base, sequence)
    assert len(ret) == len(sequence)
    return ret

def pad_sequence(sequence, source_sequence, length=SEGMENT_LENGTH):
    assert len(sequence) < length
    assert sequence == source_sequence or len(source_sequence) > length
    if len(source_sequence) > length:
        ret = source_sequence[-length:]
    else:
        ret = (source_sequence * ((length // len(sequence)) + 1))[:length]
    assert len(ret) == length
    return ret

def encode_sequence(sequence, nuc_order=None):
    if nuc_order is None:
        nuc_order = DEFAULT_NUC_ORDER
    seq = sequence.upper()[:SEGMENT_LENGTH]
    assert re.match('^[ATCG]+$', seq)
    encoded = [(nuc_order[x] + 1) for x in seq]
    return np.array([encoded])

def encode_sequences(sequences, nuc_order=None, segment_length=SEGMENT_LENGTH):
    encoded = [encode_sequence(s, nuc_order)[0] for s in sequences]
    return np.array(pad_sequences(encoded, maxlen=segment_length, padding='post'))

def save_model(model, model_path):
    model.save(model_path + '.h5')

def load_model_keras(model_path):
    return load_model(model_path)

# Assembly functions --------------------------------------------------------------------------------------------------
def assemble_right(read, read_list, score_list, score_read=1, sc_thr=0.5, runs=5000, sublen=SEARCHSUBLEN):
    # ... original implementation ...
    pass  # placeholder

def assemble_left(read, read_list, score_list, scores_read, sc_thr=0.5, runs=5000, sublen=SEARCHSUBLEN):
    # ... original implementation ...
    pass  # placeholder

def assemble_read(read, read_list, score_list, score_read):
    # ... original implementation ...
    pass  # placeholder

def assemble_read_loop(readsv, reads0, scores0, scoresv, filen, lenthr=48):
    # ... original implementation ...
    pass  # placeholder

def load_virus_model(model_path):
    return load_model_keras(model_path)

def filter_sequences(seqs):
    bad = ['A'*SEARCHSUBLEN, 'C'*SEARCHSUBLEN, 'G'*SEARCHSUBLEN, 'T'*SEARCHSUBLEN]
    return [s for s in seqs if all(b not in s for b in bad)]

def proc_fastq(infile):
    with open(infile) as f:
        lines = f.readlines()
    seqs = [handle_non_ATGC(lines[i].strip()) for i in range(1, len(lines), 4)]
    seqs = list(np.unique(seqs))
    seqs = filter_sequences(seqs)
    med = int(np.median([len(s) for s in seqs]))
    seqs = [s[:med] for s in seqs]
    return encode_sequences(seqs), seqs

def make_clist(lst):
    return (c_char_p * len(lst))(*[x.encode() for x in lst])

def assemble_read_call_c(readsv, reads0, scores0, scoresv, filen):
    # ... original implementation ...
    pass  # placeholder

# MAIN FUNCTION WITH SEED DUMP -----------------------------------------------------------------------------------------
def extract_contigs(invars, large_file_thr=1000000):
    inpath, outpath, fastmode, model_path = invars

    # Prepare output filenames
    base = os.path.basename(inpath)
    contig_fn = os.path.join(outpath, base.replace('_unmapped.fastq', '_contigs.txt'))
    seed_fn  = os.path.join(outpath, base.replace('_unmapped.fastq', '_seeds.txt'))

    # Skip if already processed
    if exists(contig_fn) or exists(seed_fn):
        return 0

    # Load and score reads
    model = load_virus_model(model_path)
    encoded_c, seqs = proc_fastq(inpath)
    scores = list(model.predict(encoded_c))

    # Select seeds
    threshold = 0.7
    seeds = [(seqs[i], scores[i]) for i in range(len(scores)) if scores[i] > threshold]

    # Dump seeds to file
    with open(seed_fn, 'w') as sf:
        for seq, score in seeds:
            sf.write(f">{score}\n{seq}\n")

    print(f"Seed reads dumped to {seed_fn}. Skipping contig assembly here.")
    return 1


def run_virna_pred(inpath, outpath, fastmode, multi_proc, model_path, num_threads):
    infastq = list(set(glob.glob(os.path.join(inpath, '*.fastq'))))
    outs = glob.glob(os.path.join(outpath, '*.txt'))
    processed = {os.path.basename(o).replace('_contigs.txt','').replace('_seeds.txt','') for o in outs}
    to_process = [f for f in infastq if os.path.basename(f).replace('_unmapped.fastq','') not in processed]

    print('starting_prediction...')
    if multi_proc:
        freeze_support()
        pool = mp.Pool(processes=num_threads)
        pool.map(extract_contigs, [[f, outpath, fastmode, model_path] for f in to_process])
    else:
        for f in to_process:
            extract_contigs([f, outpath, fastmode, model_path])

    print('Done processing')
