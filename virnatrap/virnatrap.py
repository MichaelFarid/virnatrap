"""
Model classes and functions to identify viral reads and dump seed windows
"""
# Imports --------------------------------------------------------------------------------------------------------------
import os
import re
import sys
import glob
from multiprocessing import freeze_support
import multiprocessing as mp
import numpy as np
import random
import tensorflow as tf
import ctypes
from ctypes import c_char_p, c_int, CDLL
from tensorflow import get_logger, autograph
from tensorflow.keras.models import load_model
from tensorflow.keras.preprocessing.sequence import pad_sequences
from os.path import exists, join, basename

# Default threshold, overridden at runtime by run_virna_pred
THRESHOLD = 0.7

# Suppress TensorFlow logs
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
get_logger().setLevel('ERROR')
autograph.set_verbosity(0)

# Globals -------------------------------------------------------------------------------------------------------------
DEFAULT_NUC_ORDER = {y: x for x, y in enumerate(["A", "T", "C", "G"])}
NUCLEOTIDES = list(DEFAULT_NUC_ORDER.keys())
SEGMENT_LENGTH = 48  # window size for seed extraction
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
    assert re.match('^[ATCG]+$', seq)
    encoded = [(nuc_order[x] + 1) for x in seq]
    return np.array([encoded])

def encode_sequences(sequences, nuc_order=None, segment_length=SEGMENT_LENGTH):
    encoded = [encode_sequence(s, nuc_order)[0] for s in sequences]
    return np.array(pad_sequences(encoded, maxlen=segment_length, padding='post'))

def load_virus_model(model_path):
    return load_model(model_path)

def filter_sequences(seqs):
    bad = ["A"*SEGMENT_LENGTH, "C"*SEGMENT_LENGTH, "G"*SEGMENT_LENGTH, "T"*SEGMENT_LENGTH]
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

# Placeholder assembly functions --------------------------------------------------------------------------------------
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
    return librd.assemble_read_loop(arr_f, arr_fv, arr_ch, arr_chv, n, m, len(scoresv), filen_c)

# Main functions ------------------------------------------------------------------------------------------------------
def extract_contigs(invars, large_file_thr=1000000):
    inpath, outpath, fastmode, model_path = invars
    base = basename(inpath)
    sample = base.replace('_unmapped.fastq', '')
    seed_fasta = join(outpath, f"{sample}_seeds.txt")

    # Skip if already dumped
    if exists(seed_fasta):
        return seed_fasta

    # Load model and process reads
    model = load_virus_model(model_path)
    encodings, seqs = proc_fastq(inpath)
    # Flatten predictions to Python floats
    scores = model.predict(encodings).flatten().tolist()

    # Select seeds above threshold
    seeds = [(seqs[i], scores[i]) for i in range(len(scores)) if scores[i] > THRESHOLD]

    # Dump seed windows in FASTA
    with open(seed_fasta, 'w') as sf:
        for idx, (seq, score) in enumerate(seeds):
            first48 = seq[:SEGMENT_LENGTH]
            last48  = seq[-SEGMENT_LENGTH:]
            sf.write(f">{sample}_seed{idx}_1_thr{score:.4f}\n{first48}\n")
            sf.write(f">{sample}_seed{idx}_2_thr{score:.4f}\n{last48}\n")

    print(f"Seeds dumped to {seed_fasta} with threshold {THRESHOLD}.")
    return seed_fasta


def run_virna_pred(inpath, outpath, fastmode, multi_proc, model_path, num_threads, threshold=0.7):
    global THRESHOLD
    THRESHOLD = threshold

    fastq_files = glob.glob(join(inpath, '*.fastq'))
    processed = {basename(f).split('_seeds.txt')[0] for f in glob.glob(join(outpath, '*_seeds.txt'))}

    print('starting seed extraction...')
    if multi_proc:
        freeze_support()
        pool = mp.Pool(processes=num_threads)
        pool.map(extract_contigs, [[f, outpath, fastmode, model_path] for f in fastq_files if basename(f).replace('_unmapped.fastq','') not in processed])
    else:
        for f in fastq_files:
            sample = basename(f).replace('_unmapped.fastq','')
            if sample in processed:
                continue
            extract_contigs([f, outpath, fastmode, model_path])

    print('Seed extraction complete.')
