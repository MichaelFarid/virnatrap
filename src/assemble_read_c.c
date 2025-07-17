#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <omp.h>

#define SCORE_THR 0.5f
#define RUNS 100
#define SUBLEN 24
#define SEGMENT_LENGTH 48
#define MAX_CONTIG_LEN (MAX_READS * SEGMENT_LENGTH * 2)
#define TABLE_SIZE 2000003  // prime for k-mer hash table
#define MAX_READS 10000000  // adjust according to your dataset

// Simple chained hash table for k-mer → list of read indices
typedef struct kmer_entry {
    uint64_t key;
    int *idxs;
    int count;
    struct kmer_entry *next;
} kmer_entry;
static kmer_entry *hashtable[TABLE_SIZE];

// Encode 24-mer to 48-bit integer (2 bits per base)
static inline uint64_t encode_kmer(const char *kmer) {
    uint64_t x = 0;
    for (int i = 0; i < SUBLEN; ++i) {
        x <<= 2;
        switch (kmer[i]) {
            case 'A': x |= 0; break;
            case 'C': x |= 1; break;
            case 'G': x |= 2; break;
            case 'T': x |= 3; break;
            default: x |= 0;      // treat N as A
        }
    }
    return x;
}

// Add occurrence of k-mer
void add_kmer_entry(uint64_t key, int idx) {
    uint32_t h = key % TABLE_SIZE;
    kmer_entry *e = hashtable[h];
    while (e) {
        if (e->key == key) {
            e->idxs = realloc(e->idxs, (e->count + 1) * sizeof(int));
            e->idxs[e->count++] = idx;
            return;
        }
        e = e->next;
    }
    e = malloc(sizeof(*e));
    e->key = key;
    e->count = 1;
    e->idxs = malloc(sizeof(int));
    e->idxs[0] = idx;
    e->next = hashtable[h];
    hashtable[h] = e;
}

// Retrieve entry for k-mer
kmer_entry *get_kmer_entry(uint64_t key) {
    uint32_t h = key % TABLE_SIZE;
    kmer_entry *e = hashtable[h];
    while (e) {
        if (e->key == key) return e;
        e = e->next;
    }
    return NULL;
}

// Build k-mer index from all reads
void build_kmer_index(char **reads, int num_reads) {
    char buf[SUBLEN+1];
    for (int i = 0; i < num_reads; ++i) {
        for (int j = 0; j <= SEGMENT_LENGTH - SUBLEN; ++j) {
            memcpy(buf, reads[i] + j, SUBLEN);
            buf[SUBLEN] = '\0';
            uint64_t key = encode_kmer(buf);
            add_kmer_entry(key, i);
        }
    }
}

// Find best extension to the right using hash
bool find_sub_right(int *used, char *sb0, char **reads, int num_reads, int *best_idx, int *best_pos) {
    uint64_t key = encode_kmer(sb0);
    kmer_entry *e = get_kmer_entry(key);
    if (!e) return false;
    int min_pos = SEGMENT_LENGTH;
    int idx = -1;
    for (int i = 0; i < e->count; ++i) {
        int rid = e->idxs[i];
        if (!used[rid]) {
            char *ptr = strstr(reads[rid], sb0);
            if (ptr) {
                int pos = ptr - reads[rid];
                if (pos < min_pos && pos >= 0) {
                    min_pos = pos;
                    idx = rid;
                }
            }
        }
    }
    if (idx >= 0) {
        *best_idx = idx;
        *best_pos = min_pos;
        return true;
    }
    return false;
}

// Find best extension to the left
bool find_sub_left(int *used, char *sb0, char **reads, int num_reads, int *best_idx, int *best_pos) {
    uint64_t key = encode_kmer(sb0);
    kmer_entry *e = get_kmer_entry(key);
    if (!e) return false;
    int max_pos = -1;
    int idx = -1;
    for (int i = 0; i < e->count; ++i) {
        int rid = e->idxs[i];
        if (!used[rid]) {
            char *ptr = strstr(reads[rid], sb0);
            if (ptr) {
                int pos = ptr - reads[rid];
                if (pos > max_pos) {
                    max_pos = pos;
                    idx = rid;
                }
            }
        }
    }
    if (idx >= 0) {
        *best_idx = idx;
        *best_pos = max_pos;
        return true;
    }
    return false;
}

// Contig assembly result struct
struct ret { char *c; int len; int numel; float totsc; };

// Extend right
struct ret assemble_right(const char *seed, char **reads, float *scores, int num_reads, int *used) {
    char *contig = malloc(MAX_CONTIG_LEN);
    memcpy(contig, seed, SEGMENT_LENGTH);
    int clen = SEGMENT_LENGTH;
    float totsc = scores[0];
    int numel = 1;
    char sb0[SUBLEN+1];
    memcpy(sb0, contig + clen - SUBLEN, SUBLEN);
    sb0[SUBLEN] = '\0';
    int cnt = 0;
    while (cnt < RUNS && totsc/numel > SCORE_THR) {
        int idx, pos;
        if (!find_sub_right(used, sb0, reads, num_reads, &idx, &pos)) break;
        used[idx] = 1;
        memcpy(contig + clen, reads[idx] + pos + SUBLEN, SEGMENT_LENGTH);
        clen += SEGMENT_LENGTH;
        memcpy(sb0, contig + clen - SUBLEN, SUBLEN);
        totsc += scores[idx]; numel++; cnt++;
    }
    contig[clen] = '\0';
    return (struct ret){contig, clen, numel, totsc};
}

// Extend left
struct ret assemble_left(const char *seed, char **reads, float *scores, int num_reads, int *used) {
    char *contig = malloc(MAX_CONTIG_LEN);
    int start = MAX_CONTIG_LEN - SEGMENT_LENGTH;
    memcpy(contig + start, seed, SEGMENT_LENGTH);
    int clen = SEGMENT_LENGTH;
    float totsc = scores[0];
    int numel = 1;
    char sb0[SUBLEN+1];
    memcpy(sb0, seed, SUBLEN);
    sb0[SUBLEN] = '\0';
    int cnt = 0;
    while (cnt < RUNS && totsc/numel > SCORE_THR) {
        int idx, pos;
        if (!find_sub_left(used, sb0, reads, num_reads, &idx, &pos)) break;
        used[idx] = 1;
        int addlen = pos;
        memcpy(contig + start - addlen, reads[idx], addlen);
        start -= addlen; clen += addlen;
        memcpy(sb0, reads[idx], SUBLEN);
        totsc += scores[idx]; numel++; cnt++;
    }
    char *out = malloc(clen + 1);
    memcpy(out, contig + start, clen);
    out[clen] = '\0';
    free(contig);
    return (struct ret){out, clen, numel, totsc};
}

// Parallel loop over seeds
int assemble_read_loop(float *f_arr, float *f_arr2, char *ch_arr[], char *ch_arr2[], int num_reads, int nvr, const char *fname) {
    build_kmer_index(ch_arr, num_reads);
    omp_lock_t lock;
    omp_init_lock(&lock);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nvr; ++i) {
        int *used = calloc(num_reads, sizeof(int));
        struct ret rr = assemble_right(ch_arr2[i], ch_arr, f_arr2, num_reads, used);
        struct ret rl = assemble_left(ch_arr2[i], ch_arr, f_arr2, num_reads, used);
        int total_len = rl.len + rr.len;
        if (total_len >= SEGMENT_LENGTH && ((rl.totsc/rl.numel + rr.totsc/rr.numel)/2) > SCORE_THR) {
            char *full = malloc(total_len + 1);
            memcpy(full, rl.c, rl.len);
            memcpy(full + rl.len, rr.c, rr.len);
            full[total_len] = '\0';
            omp_set_lock(&lock);
            FILE *fp = fopen(fname, "a");
            fprintf(fp, ">contig_%d[]\n%s\n", i, full);
            fclose(fp);
            omp_unset_lock(&lock);
            free(full);
        }
        free(rr.c);
        free(rl.c);
        free(used);
    }
    omp_destroy_lock(&lock);
    // free hash table
    for (int h = 0; h < TABLE_SIZE; ++h) {
        kmer_entry *e = hashtable[h];
        while (e) {
            kmer_entry *tmp = e->next;
            free(e->idxs);
            free(e);
            e = tmp;
        }
    }
    return 1;
}
