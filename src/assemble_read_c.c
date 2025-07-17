#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

#define SCORE_THR       0.5f
#define RUNS            100
#define SUBLEN          24
#define SEGMENT_LENGTH  48
#define MAX_CONTIG_LEN  ((2*RUNS + 1) * SEGMENT_LENGTH)
#define TABLE_SIZE      2000003

// Hash table entry for k-mers
typedef struct kmer_entry {
    uint64_t key;
    int *idxs;
    int count;
    struct kmer_entry *next;
} kmer_entry;
static kmer_entry *hashtable[TABLE_SIZE] = {NULL};

typedef uint64_t kmer_t;

// Encode 24-mer into 2-bit packed integer
static inline kmer_t encode_kmer(const char *k) {
    kmer_t x = 0;
    for (int i = 0; i < SUBLEN; ++i) {
        x <<= 2;
        switch (k[i]) {
            case 'C': x |= 1; break;
            case 'G': x |= 2; break;
            case 'T': x |= 3; break;
            default: /* A or N */ break;
        }
    }
    return x;
}

// Insert k-mer occurrence
static void add_kmer(kmer_t key, int idx) {
    uint32_t h = key % TABLE_SIZE;
    for (kmer_entry *e = hashtable[h]; e; e = e->next) {
        if (e->key == key) {
            e->idxs = realloc(e->idxs, (e->count+1)*sizeof(int));
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

// Retrieve bucket
static kmer_entry *get_bucket(kmer_t key) {
    uint32_t h = key % TABLE_SIZE;
    for (kmer_entry *e = hashtable[h]; e; e = e->next)
        if (e->key == key) return e;
    return NULL;
}

// Build index from reads
static void build_kmer_index(char **reads, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= SEGMENT_LENGTH - SUBLEN; ++j)
            add_kmer(encode_kmer(reads[i]+j), i);
}

// Extend contig to the right
typedef struct { char *c; int len; int numel; float totsc; } ret_t;
static ret_t assemble_right(const char *seed, char **reads, float *scores, int n, int *used) {
    int sl = strlen(seed);
    char *contig = malloc(MAX_CONTIG_LEN);
    memcpy(contig, seed, sl);
    int clen = sl;
    float totsc = scores[0]; int numel = 1;
    char sb[SUBLEN+1];
    memcpy(sb, seed + sl - SUBLEN, SUBLEN); sb[SUBLEN] = '\0';
    for (int it = 0; it < RUNS && totsc/numel > SCORE_THR; ++it) {
        int best, pos;
        if (!find_sub_right(used, sb, reads, n, &best, &pos)) break;
        used[best] = 1;
        int alen = strlen(reads[best]) - (pos+SUBLEN);
        memcpy(contig+clen, reads[best]+pos+SUBLEN, alen);
        clen += alen;
        memcpy(sb, contig+clen-SUBLEN, SUBLEN); sb[SUBLEN]='\0';
        totsc += scores[best]; numel++;
    }
    contig[clen] = '\0';
    return (ret_t){contig, clen, numel, totsc};
}

// Extend contig to the left
static ret_t assemble_left(const char *seed, char **reads, float *scores, int n, int *used) {
    int sl = strlen(seed);
    char *contig = malloc(MAX_CONTIG_LEN);
    int start = MAX_CONTIG_LEN - sl;
    memcpy(contig+start, seed, sl);
    int clen=sl; float totsc=scores[0]; int numel=1;
    char sb[SUBLEN+1];
    memcpy(sb, seed, SUBLEN); sb[SUBLEN]='\0';
    for (int it=0; it<RUNS && totsc/numel>SCORE_THR; ++it) {
        int best,pos;
        if (!find_sub_left(used,sb,reads,n,&best,&pos)) break;
        used[best]=1;
        memcpy(contig+start-pos,reads[best],pos);
        start-=pos; clen+=pos;
        memcpy(sb,reads[best],SUBLEN); sb[SUBLEN]='\0';
        totsc+=scores[best]; numel++;
    }
    char *out=malloc(clen+1);
    memcpy(out,contig+start,clen); out[clen]='\0'; free(contig);
    return (ret_t){out,clen,numel,totsc};
}

// Main entry
int assemble_read_loop(float *f_arr, float *f_arr2, char *ch_arr[], char *ch_arr2[], int seg_len, int n, int nvr, const char *fname) {
    build_kmer_index(ch_arr,n);
    char out[1024]; strncpy(out,fname,1023); out[1023]='\0';
    FILE *fp=fopen(out,"w"); if(!fp){perror("fopen");return-1;}
    for(int i=0;i<nvr;i++){
        int *usedr=calloc(n,sizeof(int)), *usedl=calloc(n,sizeof(int));
        for(int j=0;j<n;j++) if(!strcmp(ch_arr[j],ch_arr2[i])) usedr[j]=usedl[j]=1;
        ret_t rr=assemble_right(ch_arr2[i],ch_arr,f_arr2,n,usedr);
        ret_t rl=assemble_left(ch_arr2[i],ch_arr,f_arr2,n,usedl);
        free(usedr); free(usedl);
        int sl=strlen(ch_arr2[i]); int lel=rl.len-sl; int rlen=rr.len;
        int fl=lel+rlen; float sc=(rl.totsc/rl.numel+rr.totsc/rr.numel)*.5f;
        if(fl>=SEGMENT_LENGTH&&sc>SCORE_THR) {
            char *full=malloc(fl+1);
            memcpy(full,rl.c,lel);
            memcpy(full+lel,rr.c,rlen);
            full[fl]='\0';
            fprintf(fp,">contig_%d[]\n%s\n",i,full);
            free(full);
        }
        free(rr.c); free(rl.c);
    }
    fclose(fp);
    for(int h=0;h<TABLE_SIZE;h++)for(kmer_entry*e=hashtable[h];e;){kmer_entry*t=e->next;free(e->idxs);free(e);e=t;}
    return 1;
}
char* assemble_read(){return NULL;}
