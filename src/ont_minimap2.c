/*
 * ont_minimap2.c
 *
 *  Created on: August 10, 2017
 *  Proprietary and confidential information of Oxford Nanopore Technologies, Limited
 *  All rights reserved; (c)2017: Oxford Nanopore Technologies, Limited
 */

#include <stdio.h>
#include <zlib.h>
#include "minimap.h"
#include "bseq.h"
#include "kseq.h"
#include "mmpriv.h"
KSEQ_INIT(gzFile, gzread)

#include "ont_minimap2.h"

// These are prototypes of methods in minimap2's index.c
int64_t mm_idx_is_idx(const char *fn);
mm_idx_t *mm_idx_gen(mm_bseq_file_t *fp, int w, int k, int b,
                     int flag, int mini_batch_size, int n_threads,
                     uint64_t batch_size);
mm_idx_t *mm_idx_load(FILE *fp);


static inline void str_enlarge(kstring_t *s, int l) {
    // copy from format.c:7
    if (s->l + l + 1 > s->m) {
        s->m = s->l + l + 1;
        kroundup32(s->m);
        s->s = (char*)realloc(s->s, s->m);
    }
}


static inline void str_copy(kstring_t *s, const char *st, const char *en) {
    // copy from format.c:16
    str_enlarge(s, en - st);
    memcpy(&s->s[s->l], st, en - st);
    s->l += en - st;
}


mm_idx_t *ontmm_load_index(const char *index_filename) {
    // disable printing to stdout
    mm_verbose = 0;

    // see main.c:225sq.
    int64_t is_idx = mm_idx_is_idx(index_filename);
    if (is_idx < 0) return NULL;

    FILE *fpr = NULL;
    mm_bseq_file_t *fp = NULL;
    if (is_idx && !(fpr = fopen(index_filename, "rb"))) return NULL;
    if (!is_idx) fp = mm_bseq_open(index_filename);

    // parameters for index generation
    int k = 15, bucket_bits = 14, n_threads = 3, is_hpc = 0;
    int w = (int)(.6666667 * k + .499);
    int minibatch_size = 200000000;
    // Always build a single unsharded index
    uint64_t batch_size = 0x7fffffffffffffffULL;

    mm_idx_t *mi = NULL;
    if (fpr) mi = mm_idx_load(fpr);
    else if (fp) mi = mm_idx_gen(fp, w, k, bucket_bits, is_hpc, minibatch_size,
                                 n_threads, batch_size);

    if (fpr) fclose(fpr);
    if (fp) mm_bseq_close(fp);
    return mi;
}


int32_t ontmm_cache_idx_occ_thres(const mm_idx_t *index) {
    // opt.mid_occ_frac=2e-4f is the only initialisation value required from mm_mapopt_init
    return mm_idx_cal_max_occ(index, 2e-4f);
}


char *ontmm_align(mm_bseq1_t query, const mm_idx_t *index, const int32_t idx_mid_occ, int type) {
    // disable printing to stdout
    mm_verbose = 0;

    // see example.c:22
    mm_mapopt_t opt;
    mm_mapopt_init(&opt); // initialize mapping parameters
    opt.mid_occ = idx_mid_occ; // cached value as slow to calculate and invariant
    mm_mapopt_update(&opt, index); // this sets the maximum minimizer occurrence
    if (type == ONTMM_FULL_ALIGNMENT) {
        if (index->flag & MM_I_NO_SEQ) {
            return strdup("ERROR:full alignment requested, but index does not support full alignment.");
        }
        opt.flag |= MM_F_CIGAR;
    }
    else if (type != ONTMM_COARSE_ALIGNMENT) {
        return strdup("ERROR:invalid value for type argument.");
    }
    mm_tbuf_t *tbuf = mm_tbuf_init();

    // get all hits for the query
    const mm_reg1_t *reg;
    int j, n_reg;
    reg = mm_map(index, query.l_seq, query.seq, &n_reg, tbuf, &opt, 0);

    // write alignment string
    kstring_t alignment_string;
    alignment_string.l = alignment_string.m = 0, alignment_string.s = 0;
    if (!n_reg) {
        if (opt.flag & MM_F_CIGAR) {
            char *empty_sam = "\t4\t*\t0\t0\t*\t*\t0\t0\n";
            str_copy(&alignment_string, query.name, query.name + strlen(query.name));
            str_copy(&alignment_string, empty_sam, empty_sam + strlen(empty_sam));
        }
        else {
            char *empty_paf = "\t0\t0\t0\t+\t*\t0\t0\t0\t0\t0\t0\n";
            str_copy(&alignment_string, query.name, query.name + strlen(query.name));
            str_copy(&alignment_string, empty_paf, empty_paf + strlen(empty_paf));
        }
    }
    else {
        // traverse hits and print them out
        for (j = 0; j < n_reg; ++j) {
            const mm_reg1_t *r = &reg[j];
            #ifndef NDEBUG
            printf("[ont_minimap2 lib] Score: %d\n", r->score);
            #endif
            if (!r->p && (opt.flag & MM_F_CIGAR)) {
                return strdup("ERROR:alignment has not returned cigar string.");
            }
            kstring_t alignment_line;
            alignment_line.l = alignment_line.m = 0, alignment_line.s = 0;
            // Get the reg_idx for the updated 'mm_write_sam2' interface
            int i;
            for (i = 0; i < n_reg; ++i) {
                if (r == &reg[i]) break;
            }
            if (opt.flag & MM_F_CIGAR) {
                mm_write_sam2(&alignment_line, index, &query, 0, i, 1, &n_reg, &reg, NULL, MM_F_OUT_MD);
            }
            else {
                mm_write_paf(&alignment_line, index, &query, reg, NULL, MM_F_OUT_MD);
            }
            #ifndef NDEBUG
            printf("[ont_minimap2 lib] Alignment line: %s\n", alignment_line.s);
            #endif
            str_copy(&alignment_string, alignment_line.s, alignment_line.s + alignment_line.l);
            char *newline = "\n";
            str_copy(&alignment_string, newline, newline + 1);
            // free the alignment line.
            free(alignment_line.s);
        }
        // Free the mm_reg1_t structures.
        for (j = 0; j < n_reg; ++j) {
            mm_reg1_t *r = &reg[j];
            free(r->p);
        }
        free(reg);
    }
    char *endstring = "\0";
    str_copy(&alignment_string, endstring, endstring + 1);
    mm_tbuf_destroy(tbuf);
    return alignment_string.s;
}


void ontmm_unload_index(mm_idx_t *index) {
    mm_idx_destroy(index);
}
