#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <omp.h>

#define SCORE_THR       0.5f
#define RUNS            100
#define SUBLEN          24
#define SEGMENT_LENGTH  48
#define MAX_CONTIG_LEN  ((2*RUNS + 1) * SEGMENT_LENGTH)
#define TABLE_SIZE      2000003

// k-mer hash table entry
typedef struct kmer_entry {
    uint64_t key;
    int     *idxs;
    int      count;
    struct kmer_entry *next;
} kmer_entry;
static kmer_entry *hashtable[TABLE_SIZE] = {NULL};

// Pack a SUBLEN-length k-mer into a 2-bit integer
typedef uint64_t kmer_t;
static inline kmer_t encode_kmer(const char *kmer) {
    kmer_t x = 0;
    for (int i = 0; i < SUBLEN; ++i) {
        x <<= 2;
        switch (kmer[i]) {
            case 'C': x |= 1; break;
            case 'G': x |= 2; break;
            case 'T': x |= 3; break;
            default:  break;  // 'A' or 'N'
        }
    }
    return x;
}

// Insert one k-mer occurrence into hash table
static void add_kmer(kmer_t key, int idx) {
    uint32_t h = key % TABLE_SIZE;
    for (kmer_entry *e = hashtable[h]; e; e = e->next) {
        if (e->key == key) {
            e->idxs = realloc(e->idxs, (e->count + 1) * sizeof(int));
            e->idxs[e->count++] = idx;
            return;
        }
    }
    kmer_entry *e = malloc(sizeof(*e));
    e->key     = key;
    e->count   = 1;
    e->idxs    = malloc(sizeof(int));
    e->idxs[0] = idx;
    e->next    = hashtable[h];
    hashtable[h] = e;
}

// Retrieve bucket list for given k-mer, or NULL if not present
static kmer_entry *get_bucket(kmer_t key) {
    uint32_t h = key % TABLE_SIZE;
    for (kmer_entry *e = hashtable[h]; e; e = e->next) {
        if (e->key == key) return e;
    }
    return NULL;
}

// Build k-mer index for all reads (each length SEGMENT_LENGTH)
static void build_kmer_index(char **reads, int num_reads) {
    for (int i = 0; i < num_reads; ++i) {
        for (int j = 0; j <= SEGMENT_LENGTH - SUBLEN; ++j) {
            kmer_t key = encode_kmer(reads[i] + j);
            add_kmer(key, i);
        }
    }
}

// Try to extend seed to the right; return true if extension found
static bool find_sub_right(int *used, const char *sb0,
                           char **reads, int num_reads,
                           int *best_idx, int *best_pos) {
    kmer_entry *e = get_bucket(encode_kmer(sb0));
    if (!e) return false;
    int min_pos = SEGMENT_LENGTH;
    int idx = -1;
    for (int k = 0; k < e->count; ++k) {
        int rid = e->idxs[k];
        if (used[rid]) continue;
        char *ptr = strstr(reads[rid], sb0);
        if (ptr) {
            int pos = ptr - reads[rid];
            if (pos >= 0 && pos < min_pos) {
                min_pos = pos;
                idx = rid;
            }
        }
    }
    if (idx < 0) return false;
    *best_idx = idx;
    *best_pos = min_pos;
    return true;
}

// Try to extend seed to the left; return true if extension found
static bool find_sub_left(int *used, const char *sb0,
                          char **reads, int num_reads,
                          int *best_idx, int *best_pos) {
    kmer_entry *e = get_bucket(encode_kmer(sb0));
    if (!e) return false;
    int max_pos = -1;
    int idx = -1;
    for (int k = 0; k < e->count; ++k) {
        int rid = e->idxs[k];
        if (used[rid]) continue;
        char *ptr = strstr(reads[rid], sb0);
        if (ptr) {
            int pos = ptr - reads[rid];
            if (pos > max_pos) {
                max_pos = pos;
                idx = rid;
            }
        }
    }
    if (idx < 0) return false;
    *best_idx = idx;
    *best_pos = max_pos;
    return true;
}

// Contig assembly result
struct ret { char *c; int len; int numel; float totsc; };

// Extend to the right
static struct ret assemble_right(const char *seed, char **reads,
                                 float *scores, int num_reads,
                                 int *used) {
    char *contig = malloc(MAX_CONTIG_LEN);
    memcpy(contig, seed, SEGMENT_LENGTH);
    int clen = SEGMENT_LENGTH;
    float totsc = scores[0];
    int numel = 1;
    char sb0[SUBLEN+1];
    memcpy(sb0, contig + clen - SUBLEN, SUBLEN);
    sb0[SUBLEN] = '\0';
    for (int it = 0; it < RUNS && totsc/numel > SCORE_THR; ++it) {
        int idx, pos;
        if (!find_sub_right(used, sb0, reads, num_reads, &idx, &pos)) break;
        used[idx] = 1;
        memcpy(contig + clen, reads[idx] + pos + SUBLEN, SEGMENT_LENGTH);
        clen += SEGMENT_LENGTH;
        memcpy(sb0, contig + clen - SUBLEN, SUBLEN);
        totsc += scores[idx];
        numel++;
    }
    contig[clen] = '\0';
    return (struct ret){contig, clen, numel, totsc};
}

// Extend to the left
static struct ret assemble_left(const char *seed, char **reads,
                                float *scores, int num_reads,
                                int *used) {
    char *contig = malloc(MAX_CONTIG_LEN);
    int start = MAX_CONTIG_LEN - SEGMENT_LENGTH;
    memcpy(contig + start, seed, SEGMENT_LENGTH);
    int clen = SEGMENT_LENGTH;
    float totsc = scores[0];
    int numel = 1;
    char sb0[SUBLEN+1];
    memcpy(sb0, seed, SUBLEN);
    sb0[SUBLEN] = '\0';
    for (int it = 0; it < RUNS && totsc/numel > SCORE_THR; ++it) {
        int idx, pos;
        if (!find_sub_left(used, sb0, reads, num_reads, &idx, &pos)) break;
        used[idx] = 1;
        memcpy(contig + start - pos, reads[idx], pos);
        start -= pos;
        clen += pos;
        memcpy(sb0, reads[idx], SUBLEN);
        totsc += scores[idx];
        numel++;
    }
    char *out = malloc(clen + 1);
    memcpy(out, contig + start, clen);
    out[clen] = '\0';
    free(contig);
    return (struct ret){out, clen, numel, totsc};
}

// Main entry: matches Python ctypes signature
int assemble_read_loop(float *f_arr, float *f_arr2,
                       char *ch_arr[], char *ch_arr2[],
                       int seg_len, int num_reads, int nvr,
                       const char *fname) {
    char outbuf[1024];
    strncpy(outbuf, fname, sizeof(outbuf)-1);
    outbuf[sizeof(outbuf)-1] = '\0';

    build_kmer_index(ch_arr, num_reads);

    FILE *fp = fopen(outbuf, "w");
    if (!fp) { perror("fopen"); return -1; }

    char **printed = calloc(nvr, sizeof(char*));
    int printed_count = 0;
    bool *seed_used = calloc(nvr, sizeof(bool));

    for (int i = 0; i < nvr; ++i) {
        if (seed_used[i]) continue;
        int *used_left  = calloc(num_reads, sizeof(int));
        int *used_right = calloc(num_reads, sizeof(int));
        // find seed index
        int seed_idx = -1;
        for (int j = 0; j < num_reads; ++j) {
            if (strcmp(ch_arr[j], ch_arr2[i]) == 0) { seed_idx = j; break; }
        }
        if (seed_idx >= 0) {
            used_left[seed_idx] = 1;
            used_right[seed_idx] = 1;
        }
        struct ret rr = assemble_right(ch_arr2[i], ch_arr, f_arr2, num_reads, used_right);
        struct ret rl = assemble_left(ch_arr2[i], ch_arr, f_arr2, num_reads, used_left);
        int left_ext_len = rl.len - SEGMENT_LENGTH;
        int right_len    = rr.len;
        int full_len     = left_ext_len + right_len;
        float avg_sc     = (rl.totsc/rl.numel + rr.totsc/rr.numel) * 0.5f;
        if (full_len >= SEGMENT_LENGTH && avg_sc > SCORE_THR) {
            char *full = malloc(full_len + 1);
            memcpy(full, rl.c, left_ext_len);
            memcpy(full + left_ext_len, rr.c, right_len);
            full[full_len] = '\0';
            bool is_sub = false;
            for (int k = 0; k < printed_count; ++k) {
                if (strstr(printed[k], full)) { is_sub = true; break; }
            }
            if (!is_sub) {
                printed[printed_count++] = full;
                fprintf(fp, ">contig_%d[]\n%s\n", i, full);
                for (int j = 0; j < nvr; ++j) {
                    if (!seed_used[j] && strstr(full, ch_arr2[j])) {
                        seed_used[j] = true;
                    }
                }
            } else {
                free(full);
            }
        }
        free(rr.c);
        free(rl.c);
        free(used_left);
        free(used_right);
    }

    for (int k = 0; k < printed_count; ++k) free(printed[k]);
    free(printed);
    free(seed_used);
    fclose(fp);

    // Cleanup hash table
    for (int h = 0; h < TABLE_SIZE; ++h) {
        for (kmer_entry *e = hashtable[h]; e; ) {
            kmer_entry *tmp = e->next;
            free(e->idxs);
            free(e);
            e = tmp;
        }
    }

    return 1;
}

// Stub for Python ctypes lookup
char* assemble_read() { return NULL; }
