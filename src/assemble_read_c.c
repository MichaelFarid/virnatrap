#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <omp.h>
#include "uthash.h"

#define SCORE_THR 0.5f
#define RUNS 100
#define SUBLEN 24
#define SEGMENT_LENGTH 48
#define MAX_CONTIG_LEN 1000000

// K-mer hash entry
typedef struct {
    char kmer[SUBLEN+1];
    int *idxs;
    int count;
    UT_hash_handle hh;
} kmer_entry;

static kmer_entry *kmer_index = NULL;

// Add a k-mer occurrence
void add_kmer(const char *kmer, int idx) {
    kmer_entry *e;
    HASH_FIND_STR(kmer_index, kmer, e);
    if (!e) {
        e = malloc(sizeof(*e));
        strcpy(e->kmer, kmer);
        e->idxs = malloc(sizeof(int));
        e->count = 1;
        e->idxs[0] = idx;
        HASH_ADD_STR(kmer_index, kmer, e);
    } else {
        e->idxs = realloc(e->idxs, (e->count + 1) * sizeof(int));
        e->idxs[e->count++] = idx;
    }
}

// Build k-mer index from all reads
void build_kmer_index(char **reads, int num_reads) {
    char buf[SUBLEN+1];
    for (int i = 0; i < num_reads; i++) {
        for (int j = 0; j <= SEGMENT_LENGTH - SUBLEN; j++) {
            memcpy(buf, reads[i] + j, SUBLEN);
            buf[SUBLEN] = '\0';
            add_kmer(buf, i);
        }
    }
}

// Find best extension to the right
bool find_sub_right(int *read_used, const char *sb0, int *best_idx, int *best_pos, char **reads) {
    kmer_entry *e;
    HASH_FIND_STR(kmer_index, sb0, e);
    if (!e) return false;
    int min_pos = SEGMENT_LENGTH;
    int idx = -1;
    for (int i = 0; i < e->count; i++) {
        int rid = e->idxs[i];
        if (!read_used[rid]) {
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
bool find_sub_left(int *read_used, const char *sb0, int *best_idx, int *best_pos, char **reads) {
    kmer_entry *e;
    HASH_FIND_STR(kmer_index, sb0, e);
    if (!e) return false;
    int max_pos = -1;
    int idx = -1;
    for (int i = 0; i < e->count; i++) {
        int rid = e->idxs[i];
        if (!read_used[rid]) {
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

// Contig assembly result
struct ret { char *c; int len; int numel; float totsc; };

// Extend contig to the right
struct ret assemble_right(const char *seed, char **reads, float *scores, int num_reads, int *read_used) {
    char *contig = malloc(MAX_CONTIG_LEN);
    int clen = SEGMENT_LENGTH;
    memcpy(contig, seed, SEGMENT_LENGTH);
    float totsc = scores[0];
    int numel = 1;
    char sb0[SUBLEN+1];
    memcpy(sb0, contig + clen - SUBLEN, SUBLEN);
    sb0[SUBLEN] = '\0';
    int cnt = 0;
    while (cnt < RUNS && totsc / numel > SCORE_THR) {
        int idx, pos;
        if (!find_sub_right(read_used, sb0, &idx, &pos, reads)) break;
        read_used[idx] = 1;
        memcpy(contig + clen, reads[idx] + pos + SUBLEN, SEGMENT_LENGTH);
        clen += SEGMENT_LENGTH;
        memcpy(sb0, contig + clen - SUBLEN, SUBLEN);
        totsc += scores[idx];
        numel++;
        cnt++;
    }
    contig[clen] = '\0';
    return (struct ret){contig, clen, numel, totsc};
}

// Extend contig to the left
struct ret assemble_left(const char *seed, char **reads, float *scores, int num_reads, int *read_used) {
    char *contig = malloc(MAX_CONTIG_LEN);
    int clen = 0;
    float totsc = scores[0];
    int numel = 1;
    // Copy seed at end of buffer
    memcpy(contig + MAX_CONTIG_LEN - SEGMENT_LENGTH, seed, SEGMENT_LENGTH);
    clen = SEGMENT_LENGTH;
    char sb0[SUBLEN+1];
    memcpy(sb0, seed, SUBLEN);
    sb0[SUBLEN] = '\0';
    int cnt = 0;
    while (cnt < RUNS && totsc / numel > SCORE_THR) {
        int idx, pos;
        if (!find_sub_left(read_used, sb0, &idx, &pos, reads)) break;
        read_used[idx] = 1;
        int add_len = pos;
        memcpy(contig + MAX_CONTIG_LEN - SEGMENT_LENGTH - add_len, reads[idx], add_len);
        clen += add_len;
        memcpy(sb0, reads[idx], SUBLEN);
        totsc += scores[idx];
        numel++;
        cnt++;
    }
    // shift to front
    char *out = malloc(clen + SEGMENT_LENGTH + 1);
    int left_start = MAX_CONTIG_LEN - SEGMENT_LENGTH - (clen - SEGMENT_LENGTH);
    memcpy(out, contig + left_start, clen + SEGMENT_LENGTH);
    out[clen + SEGMENT_LENGTH] = '\0';
    free(contig);
    return (struct ret){out, clen + SEGMENT_LENGTH, numel, totsc};
}

// Main loop: parallel over seeds
int assemble_read_loop(float *f_arr, float *f_arr2, char *ch_arr[], char *ch_arr2[], int num_reads, int nvr, const char *fname) {
    build_kmer_index(ch_arr, num_reads);
    omp_lock_t file_lock;
    omp_init_lock(&file_lock);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nvr; i++) {
        int *used1 = calloc(num_reads, sizeof(int));
        int *used2 = calloc(nvr, sizeof(int));
        struct ret rr = assemble_right(ch_arr2[i], ch_arr, f_arr2, num_reads, used1);
        struct ret rl = assemble_left(ch_arr2[i], ch_arr, f_arr2, num_reads, used1);
        // merge left and right
        int total_len = rl.len + rr.len;
        char *full = malloc(total_len + 1);
        memcpy(full, rl.c, rl.len);
        memcpy(full + rl.len, rr.c, rr.len);
        full[total_len] = '\0';

        if (total_len > SEGMENT_LENGTH && (rl.totsc/rl.numel + rr.totsc/rr.numel)/2 > SCORE_THR) {
            omp_set_lock(&file_lock);
            FILE *fp = fopen(fname, "a");
            fprintf(fp, ">contig_%d[]\n%s\n", i, full);
            fclose(fp);
            omp_unset_lock(&file_lock);
        }

        free(rr.c);
        free(rl.c);
        free(full);
        free(used1);
        free(used2);
    }

    omp_destroy_lock(&file_lock);
    // free kmer index
    kmer_entry *e, *tmp;
    HASH_ITER(hh, kmer_index, e, tmp) {
        HASH_DEL(kmer_index, e);
        free(e->idxs);
        free(e);
    }
    return 1;
}
