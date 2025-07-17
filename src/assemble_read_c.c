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
    int *idxs;
    int count;
    struct kmer_entry *next;
} kmer_entry;
static kmer_entry *hashtable[TABLE_SIZE] = {NULL};

typedef uint64_t kmer_t;

// Encode a SUBLEN-mer into a 2-bit integer
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

// Insert one k-mer occurrence into the hash table
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
    e->key = key;
    e->count = 1;
    e->idxs = malloc(sizeof(int));
    e->idxs[0] = idx;
    e->next = hashtable[h];
    hashtable[h] = e;
}

// Retrieve bucket for a k-mer key
static kmer_entry *get_bucket(kmer_t key) {
    uint32_t h = key % TABLE_SIZE;
    for (kmer_entry *e = hashtable[h]; e; e = e->next) {
        if (e->key == key) return e;
    }
    return NULL;
}

// Build k-mer index from all reads
static void build_kmer_index(char **reads, int num_reads) {
    for (int i = 0; i < num_reads; ++i) {
        for (int j = 0; j <= SEGMENT_LENGTH - SUBLEN; ++j) {
            add_kmer(encode_kmer(reads[i] + j), i);
        }
    }
}

// Attempt right extension
static bool find_sub_right(int *used, const char *sb0,
                           char **reads, int num_reads,
                           int *best_idx, int *best_pos) {
    kmer_entry *e = get_bucket(encode_kmer(sb0));
    if (!e) return false;
    int min_pos = SEGMENT_LENGTH;
    int idx = -1;
    for (int i = 0; i < e->count; ++i) {
        int rid = e->idxs[i];
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

// Attempt left extension
static bool find_sub_left(int *used, const char *sb0,
                          char **reads, int num_reads,
                          int *best_idx, int *best_pos) {
    kmer_entry *e = get_bucket(encode_kmer(sb0));
    if (!e) return false;
    int max_pos = -1;
    int idx = -1;
    for (int i = 0; i < e->count; ++i) {
        int rid = e->idxs[i];
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

// Struct for contig assembly result
struct ret { char *c; int len; int numel; float totsc; };

// Extend contig to the right
static struct ret assemble_right(const char *seed, char **reads,
                                 float *scores, int num_reads,
                                 int *used) {
    int seed_len = strlen(seed);
    char *contig = malloc(MAX_CONTIG_LEN);
    memcpy(contig, seed, seed_len);
    int clen = seed_len;
    float totsc = scores[0];
    int numel = 1;
    char sb0[SUBLEN+1];
    // Init sb0 to last SUBLEN bases of seed
    if (seed_len >= SUBLEN) {
        memcpy(sb0, seed + seed_len - SUBLEN, SUBLEN);
    } else {
        memset(sb0, 'A', SUBLEN - seed_len);
        memcpy(sb0 + (SUBLEN - seed_len), seed, seed_len);
    }
    sb0[SUBLEN] = '\0';
    for (int it = 0; it < RUNS && totsc/numel > SCORE_THR; ++it) {
        int idx, pos;
        if (!find_sub_right(used, sb0, reads, num_reads, &idx, &pos)) break;
        used[idx] = 1;
        int append_len = strlen(reads[idx]) - (pos + SUBLEN);
        memcpy(contig + clen, reads[idx] + pos + SUBLEN, append_len);
        clen += append_len;
        if (clen >= SUBLEN) memcpy(sb0, contig + clen - SUBLEN, SUBLEN);
        sb0[SUBLEN] = '\0';
        totsc += scores[idx];
        numel++;
    }
    contig[clen] = '\0';
    return (struct ret){contig, clen, numel, totsc};
}

// Extend contig to the left
static struct ret assemble_left(const char *seed, char **reads,
                                float *scores, int num_reads,
                                int *used) {
    int seed_len = strlen(seed);
    char *contig = malloc(MAX_CONTIG_LEN);
    int start = MAX_CONTIG_LEN - seed_len;
    memcpy(contig + start, seed, seed_len);
    int clen = seed_len;
    float totsc = scores[0];
    int numel = 1;
    char sb0[SUBLEN+1];
    // Init sb0 to first SUBLEN bases of seed
    if (seed_len >= SUBLEN) {
        memcpy(sb0, seed, SUBLEN);
    } else {
        memcpy(sb0, seed, seed_len);
        memset(sb0 + seed_len, 'A', SUBLEN - seed_len);
    }
    sb0[SUBLEN] = '\0';
    for (int it = 0; it < RUNS && totsc/numel > SCORE_THR; ++it) {
        int idx, pos;
        if (!find_sub_left(used, sb0, reads, num_reads, &idx, &pos)) break;
        used[idx] = 1;
        memcpy(contig + start - pos, reads[idx], pos);
        start -= pos;
        clen += pos;
        if (pos >= SUBLEN) memcpy(sb0, reads[idx], SUBLEN);
        else {
            memcpy(sb0, reads[idx], pos);
            memset(sb0 + pos, 'A', SUBLEN - pos);
        }
        sb0[SUBLEN] = '\0';
        totsc += scores[idx];
        numel++;
    }
    char *out = malloc(clen + 1);
    memcpy(out, contig + start, clen);
    out[clen] = '\0';
    free(contig);
    return (struct ret){out, clen, numel, totsc};
}

// Main entry for Python ctypes
int assemble_read_loop(float *f_arr, float *f_arr2,
                       char *ch_arr[], char *ch_arr2[],
                       int seg_len, int num_reads, int nvr,
                       const char *fname) {
    build_kmer_index(ch_arr, num_reads);
    char outbuf[1024];
    strncpy(outbuf, fname, sizeof(outbuf)-1);
    outbuf[sizeof(outbuf)-1] = '\0';
    FILE *fp = fopen(outbuf, "w");
    if (!fp) { perror("fopen"); return -1; }
    char **printed = calloc(nvr, sizeof(char*));
    int printed_count = 0;
    for (int i = 0; i < nvr; ++i) {
        bool skip = false;
        for (int k = 0; k < printed_count; ++k) {
            if (strstr(printed[k], ch_arr2[i])) { skip = true; break; }
        }
        if (skip) continue;
        // Determine seed index in full reads
        int seed_idx = -1;
        for (int j = 0; j < num_reads; ++j) {
            if (strcmp(ch_arr[j], ch_arr2[i]) == 0) {
                seed_idx = j;
                break;
            }
        }
        // Prepare separate used arrays for each direction
        int *used_right = calloc(num_reads, sizeof(int));
        int *used_left  = calloc(num_reads, sizeof(int));
        if (seed_idx >= 0) {
            used_right[seed_idx] = 1;
            used_left[seed_idx]  = 1;
        }
        // Assemble right (seed + right)
        struct ret rr = assemble_right(ch_arr2[i], ch_arr, f_arr2, num_reads, used_right);
        // Assemble left (left + seed)
        struct ret rl = assemble_left(ch_arr2[i], ch_arr, f_arr2, num_reads, used_left);
        // Clean up used arrays
        free(used_right);
        free(used_left);
    }
    for (int k = 0; k < printed_count; ++k) free(printed[k]);
    free(printed);
    fclose(fp);
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

// Stub for ctypes
char* assemble_read() { return NULL; }
