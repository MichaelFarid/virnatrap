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
import random
import tensorflow as tf
import ctypes
from ctypes import c_char_p, c_int, CDLL
from tensorflow import get_logger, autograph
from tensorflow.keras.models import load_model
from tensorflow.keras.preprocessing.sequence import pad_sequences
import glob
from os.path import exists, join, basename

# Configuration -------------------------------------------------------------------------------------------------------
# Default threshold; override via command-line wrapper by setting virnatrap.THRESHOLD
THRESHOLD = 0.7

# Suppress TensorFlow logs
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
get_logger().setLevel('ERROR')
autograph.set_verbosity(0)

# Globals -------------------------------------------------------------------------------------------------------------
DEFAULT_NUC_ORDER = {y: x for x, y in enumerate(["A", "T", "C", "G"])}
NUCLEOTIDES = list(DEFAULT_NUC_ORDER.keys())
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

def load_virus_model(model_path):
    return load_model(model_path)

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

# Placeholder assembly functions (original implementation should be retained)
def assemble_right(*args, **kwargs): pass
def assemble_left(*args, **kwargs): pass
def assemble_read(*args, **kwargs): pass
def assemble_read_loop(*args, **kwargs): pass

def make_clist(lst):
    return (c_char_p * len(lst))(*[x.encode() for x in lst])

def assemble_read_call_c(readsv, reads0, scores0, scoresv, filen):
    librd = CDLL(join(PWD, "src/assemble_read_c.so"))
    librd.connect()
    arr_f = (ctypes.c_float * len(scores0))(*scores0)
    arr_fv = (ctypes.c_float * len(scoresv))(*scoresv)
    filen_c = c_char_p(filen.encode())
    m = c_int(len(reads0))
    n = c_int(len(readsv[0]))
    arr_ch = make_clist(reads0)
    arr_chv = make_clist(readsv)
    result = librd.assemble_read_loop(arr_f, arr_fv, arr_ch, arr_chv, n, m, len(scoresv), filen_c)
    return result

# MAIN FUNCTION WITH SEED DUMP -----------------------------------------------------------------------------------------
def extract_contigs(invars, large_file_thr=1000000):
    inpath, outpath, fastmode, model_path = invars

    base = basename(inpath)
    contig_fn = join(outpath, base.replace('_unmapped.fastq', '_contigs.txt'))
    seed_fn   = join(outpath, base.replace('_unmapped.fastq', '_seeds.txt'))

    if exists(contig_fn) or exists(seed_fn):
        return 0

    model = load_virus_model(model_path)
    encoded_c, seqs = proc_fastq(inpath)
    scores = list(model.predict(encoded_c))

    seeds = [(seqs[i], scores[i]) for i in range(len(scores)) if scores[i] > THRESHOLD]

    with open(seed_fn, 'w') as sf:
        for seq, score in seeds:
            sf.write(f">{score}\n{seq}\n")

    print(f"Seed reads dumped to {seed_fn} with threshold {THRESHOLD}. Skipping contig assembly here.")
    return 1


def run_virna_pred(inpath, outpath, fastmode, multi_proc, model_path, num_threads):
    infastq = glob.glob(join(inpath, '*.fastq'))
    outs = glob.glob(join(outpath, '*.txt'))
    processed = {basename(o).replace('_contigs.txt','').replace('_seeds.txt','') for o in outs}
    to_process = [f for f in infastq if basename(f).replace('_unmapped.fastq','') not in processed]

    print('starting_prediction...')
    if multi_proc:
        freeze_support()
        pool = mp.Pool(processes=num_threads)
        pool.map(extract_contigs, [[f, outpath, fastmode, model_path] for f in to_process])
    else:
        for f in to_process:
            extract_contigs([f, outpath, fastmode, model_path])

    print('Done processing')
