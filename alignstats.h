#ifndef _ALIGNSTATS_H
#define _ALIGNSTATS_H

#include "align.h"
#include "alignlen.h"
#include "bed.h"
#include "coverage.h"
#include "filter.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "insertsize.h"
#include "pairstats.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#ifdef HGSC_BUILD
#define USE_PTHREAD
#else
#include "config.h"
#endif

#ifdef USE_PTHREAD
#include <pthread.h>
#include <semaphore.h>
#endif

#define ALIGNSTATS_VERSION "0.10"
#define RECORD_BUFFER_SIZE 0x4000
#define RECORD_BUFFER_SECTIONS 0x10

/* TODO split or make hierarchical */
struct args {
    /* Pthread stuff */
#ifdef USE_PTHREAD
    sem_t read_sem, proc_sem;
#endif
    bam1_t **rec_buff_arr[RECORD_BUFFER_SECTIONS];
    bam1_t **read_buff;
    bam1_t **curr_buff;
    uint32_t read_buff_size[RECORD_BUFFER_SECTIONS];
    uint32_t reads_per_buffer;
    size_t read_section_idx;
    size_t process_section_idx;

    /* Options */
    bool verbose;
    bool do_alignment;
    bool do_capture;
    bool do_wgs;
    bool do_cov_mask;
    bool do_pthread;
    bool remove_dups;
    bool remove_overlaps;
    bool remove_overlaps_mc;
    bool process_unmapped;
    bool process_unmapped_done;

    /* Read processing */
    char *prev_chrom_name;
    char *curr_chrom_name;
    uint8_t min_qual;
    uint8_t min_base_qual;
    int32_t prev_chrom_idx;
    int32_t prev_mapped_chrom_idx;
    int32_t curr_chrom_idx;
    int32_t prev_chrom_len;
    int32_t curr_chrom_len;
    uint64_t num_records_processed;
    uint64_t interval;
    bool new_chrom;
    bool process_cigar;
    bool order_warn;

    /* Files and alignment */
    char *input_fn;
    hts_itr_t *iter;
    hts_idx_t *index;
    bed_t *regions;
    size_t regions_curr_chrom_idx;
    size_t regions_curr_target_idx;
    uint32_t (*read_bam_func)(struct args *);

    /* File input/output */
    FILE *output_fp;
    FILE *cov_mask_fp;
    FILE *target_fp;
    samFile *input_sf;
    bam_hdr_t *hdr;

    /* Capture/coverage data structures */
    uint32_t *coverage;
    uint8_t *target_cov;
    capture_metrics_t *cm;
    capture_metrics_t *cm_wgs;
    coverage_info_t *ci;
    coverage_info_t *ci_wgs;
    bed_t *ti;
    bed_t *cov_mask_ti;
    overlap_handler_t *olh;

    /* Calculators */
    align_metrics_t *am_all;
    align_metrics_t *am_read1;
    align_metrics_t *am_read2;
    /*align_metrics_t *am_fragment; */
    align_len_metrics_t *alm_all;
    align_len_metrics_t *alm_read1;
    align_len_metrics_t *alm_read2;
    /*align_len_metrics_t *alm_fragment; */
    pair_stats_metrics_t *psm;
    insert_size_metrics_t *ism;
    filter_counter_t *fc;

    /* Reports */
    report_t *report;
};
typedef struct args args_t;

#endif /* _ALIGNSTATS_H */
