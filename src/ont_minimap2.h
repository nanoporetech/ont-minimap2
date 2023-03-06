/*
 * ont_minimap2.h
 *
 *  Created on: August 10, 2017
 *  Proprietary and confidential information of Oxford Nanopore Technologies, Limited
 *  All rights reserved; (c)2017: Oxford Nanopore Technologies, Limited
 */
#pragma once

#include "minimap.h"
#include "bseq.h"

#if defined(_MSC_VER)
//  Microsoft
#define EXPORT __declspec(dllexport)
#define IMPORT __declspec(dllimport)
#elif defined(__GNUC__)
//  GCC
#define EXPORT __attribute__((visibility("default")))
#define IMPORT
#else
//  do nothing and hope for the best?
#define EXPORT
#define IMPORT
#pragma warning Unknown dynamic link import/export semantics.
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define ONTMM_FULL_ALIGNMENT 1
#define ONTMM_COARSE_ALIGNMENT 0

    /** Loads index and returns a pointer to it. Call ontmm_unload_index() to free  the return value.
     */
    EXPORT mm_idx_t *ontmm_load_index(const char *index_filename);

    /** Do the expensive idx_cal_max_occ operation once and cache the result
     *  This avoids recalculating when initialising opt on each call to ontmm_align
     */
    EXPORT int32_t ontmm_cache_idx_occ_thres(const mm_idx_t *index);
  
    /** Call minimap2 alignment and return an alignment SAM or PAF string (without header!).
     *
     *  @param query The query to be aligned.
     *  @param index Pointer to index to align the query to.
     *  @param idx_mid_occ Precomputed parameter returned by prior call to ontmm_cache_idx_occ_thres().
     *  @param type Specified full or coarse alignment.
     *  @return Pointer to string representation of the alignment.
     *
     *  If anything goes wrong the returned string will begin with "ERROR:" followed by an error message.
     *  Regardless, you must call free() on the returned pointer to free the memory it is using.
     *
     *  If a full alignment is requested then a SAM string will be returned. Otherwise a PAF string will
     *  be returned. If you request a full alignment and the index does not support it (the MM_I_NO_SEQ
     *  flag is set), then an error will occur.
     */
    EXPORT char *ontmm_align(mm_bseq1_t query, const mm_idx_t *index, const int32_t idx_mid_occ, const int type);

    /** Frees any memory being used by the index.
     */
    EXPORT void ontmm_unload_index(mm_idx_t *index);

    EXPORT char *ontmm_test() { return strdup("Hello World."); }

#ifdef __cplusplus
}
#endif
