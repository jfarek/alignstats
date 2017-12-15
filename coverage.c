#include "coverage.h"
#include "err.h"
#include "logging.h"
#include "print.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Capture and whole genome statistics
 */

/* Coverage info structure */

/**
 * Create and return new *coverage_info_t.
 */
coverage_info_t *coverage_info_init()
{
    coverage_info_t *ci = calloc(1, sizeof(coverage_info_t));
    die_on_alloc_fail(ci);

    ci->cov_histo_len = 0x40000; /* == 2^18 == 262144 */
    ci->cov_histo = calloc(ci->cov_histo_len, sizeof(uint64_t));
    die_on_alloc_fail(ci->cov_histo);

    return ci;
}

/**
 * Free *coverage_info_t from memory.
 */
void coverage_info_destroy(coverage_info_t *ci)
{
    if (ci != NULL) {
        free(ci->cov_histo);
        free(ci);
    }
}

/* Capture metrics structure */

/**
 * Create and return new *capture_metrics_t.
 */
capture_metrics_t *capture_metrics_init()
{
    capture_metrics_t *cm = calloc(1, sizeof(capture_metrics_t));
    die_on_alloc_fail(cm);

    return cm;
}

/**
 * Finalize capture metrics once all records are processed.
 * Calculates median coverage.
 */
void capture_metrics_finalize(capture_metrics_t *cm, coverage_info_t *ci, bed_t *ti)
{
    uint64_t k = 0, sum = 0;
    uint64_t mid = (ti == NULL ? cm->b_total : cm->b_targeted) / 2;
    bool median_set = false;

    /* Set median coverage */
    for (size_t i = 0; i < ci->cov_histo_len; ++i) {
        if (ci->cov_histo[i] > 0) {
            k += ci->cov_histo[i];
            sum += i * ci->cov_histo[i];
            if (!median_set && sum >= mid) {
                cm->c_median = i;
                median_set = true;
            }
        }
    }

    cm->c_std_dev = k > 1
        ? sqrt(((double)cm->c_sum_sq - (double)(sum * sum) / (double)k) /
               (double)(k - 1))
        : 0.0;

    /* Masked bases */
    cm->b_total -= cm->b_masked;
}

/**
 * Free *capture_metrics_t from memory.
 */
void capture_metrics_destroy(capture_metrics_t *cm)
{
    free(cm);
}

/* Capture target and coverage statistics calculation */

/**
 * Increment ci->cov_histo[cov], realloc if needed
 */
void incr_cov_histo(coverage_info_t *ci, uint32_t cov)
{
    uint64_t *tmp_cov_histo;

    if (cov + 1 > ci->cov_histo_len) {
        /* Buffer an additional 256 to reduce # of reallocs for slowly increasing coverage */
        ci->cov_histo_len = (size_t)(cov + 1 + 256);
        tmp_cov_histo = realloc(ci->cov_histo, ci->cov_histo_len * sizeof(uint64_t));

        if (tmp_cov_histo != NULL) {
            ci->cov_histo = tmp_cov_histo;
        } else {
            // kill program
            free(ci->cov_histo);
            die_on_alloc_fail(tmp_cov_histo);
            return;
        }
    }

    ++ci->cov_histo[cov];
}

/**
 * Record whole genome coverage metrics for chromosome cname.
 */
void handle_wgs_coverage(const uint32_t *coverage, capture_metrics_t *cm,
                         coverage_info_t *ci, int32_t chrom_len)
{
    uint32_t cov;

    /* for each base position in chromosome */
    for (int32_t i = 0; i < chrom_len; ++i) {
        cov = get_coverage(coverage[i]);

        /* Bases with coverage of at least 1, 10, 20, etc. */
        if (cov < 30) {
            if (cov < 10) {
                if (cov < 1) {
                    goto cov0;
                } else {
                    goto cov1;
                }
            } else {
                if (cov < 20) {
                    goto cov10;
                } else {
                    goto cov20;
                }
            }
        } else {
            if (cov < 50) {
                if (cov < 40) {
                    goto cov30;
                } else {
                    goto cov40;
                }
            } else {
                if (cov < 100) {
                    goto cov50;
                } else {
                    if (cov < 1000) {
                        goto cov100;
                    } else {
                        goto cov1000;
                    }
                }
            }
        }

cov1000:++cm->b_1000_plus_hits;
cov100: ++cm->b_100_plus_hits;
cov50:  ++cm->b_50_plus_hits;
cov40:  ++cm->b_40_plus_hits;
cov30:  ++cm->b_30_plus_hits;
cov20:  ++cm->b_20_plus_hits;
cov10:  ++cm->b_10_plus_hits;
cov1:   ++cm->b_1_plus_hits;
        cm->c_total += (uint64_t)cov;
        cm->c_sum_sq += (uint64_t)cov * (uint64_t)cov;
cov0:   incr_cov_histo(ci, cov);
    }
}

/**
 * Record capture coverage metrics for chromosome cname.
 */
void handle_target_coverage(const uint32_t *coverage, capture_metrics_t *cm,
                            coverage_info_t *ci, bed_t *ti, int32_t chrom_idx,
                            const char *chrom, int32_t chrom_len)
{
    /* Target coverage statistics */
    bool target_hit, buffer_hit;
    int32_t start, end, j, buffer_end;
    uint32_t cov;
    bed_chrom_t *tic = ti->chroms[chrom_idx];

    /* for each target */
    for (size_t i = 0; i < tic->num_targets; ++i) {
        start = tic->start_pos[i];
        end = tic->end_pos[i];
        target_hit = false;

        /* for each base position */
        for (int32_t j = start; j <= end; ++j) {
            cov = get_coverage(coverage[j]);

            /* Bases with coverage of at least 1, 10, 20, etc. */
            if (cov < 30) {
                if (cov < 10) {
                    if (cov < 1) {
                        goto tgtcov0;
                    } else {
                        goto tgtcov1;
                    }
                } else {
                    if (cov < 20) {
                        goto tgtcov10;
                    } else {
                        goto tgtcov20;
                    }
                }
            } else {
                if (cov < 50) {
                    if (cov < 40) {
                        goto tgtcov30;
                    } else {
                        goto tgtcov40;
                    }
                } else {
                    if (cov < 100) {
                        goto tgtcov50;
                    } else {
                        if (cov < 1000) {
                            goto tgtcov100;
                        } else {
                            goto tgtcov1000;
                        }
                    }
                }
            }

tgtcov1000: ++cm->b_1000_plus_hits;
tgtcov100:  ++cm->b_100_plus_hits;
tgtcov50:   ++cm->b_50_plus_hits;
tgtcov40:   ++cm->b_40_plus_hits;
tgtcov30:   ++cm->b_30_plus_hits;
tgtcov20:   ++cm->b_20_plus_hits;
tgtcov10:   ++cm->b_10_plus_hits;
tgtcov1:    ++cm->b_1_plus_hits;
            cm->c_total += (uint64_t)cov;
            cm->c_sum_sq += (uint64_t)cov * (uint64_t)cov;
            target_hit = true;
tgtcov0:    incr_cov_histo(ci, cov);
        }

        if (target_hit) {
            ++cm->t_hit;
        } else {
            /* Check if buffers were hit */
            buffer_hit = false;

            /* Left buffer */
            j = (start > BUFFER) ? start - BUFFER : 0;
            buffer_end = (start < chrom_len) ? start : chrom_len - 1;

            while (j < buffer_end) {
                cov = get_coverage(coverage[j]);
                if (cov > 0) {
                    buffer_hit = true;
                    ++cm->t_buffers_hit;
                    break;
                }
                ++j;
            }

            if (!buffer_hit) {
                /* Right buffer */
                ++end;
                j = (end > 0) ? end : 0;
                buffer_end = (end + BUFFER < chrom_len) ? end + BUFFER : chrom_len - 1;

                while (j < buffer_end) {
                    cov = get_coverage(coverage[j]);
                    if (cov > 0) {
                        ++cm->t_buffers_hit;
                        break;
                    }
                    ++j;
                }
            }
        }
    }
}

/**
 * Set values from start to end in coverage to 0.
 */
void clear_coverage(uint32_t *coverage, int32_t start, int32_t end, int32_t chrom_len)
{
    if (start < 0) {
        start = 0;
    }
    if (end >= chrom_len) {
        end = chrom_len - 1;
    }

    /* Clear target and buffer bases */
    memset(coverage + start, 0, (end - start + 1) * sizeof(uint32_t));
}

/**
 * Record contiguous regions of >= 20X coverage outside of target and buffer
 * regions.
 */
void handle_miss_reads(uint32_t *coverage, capture_metrics_t *cm, bed_t *ti,
                       int32_t chrom_idx, int32_t chrom_len)
{
    if (ti->num_targets > 0) {
        bed_chrom_t *tic = ti->chroms[chrom_idx];

        int32_t start = 0, end;

        for (size_t i = 0; i < tic->num_targets; ++i) {
            end = tic->start_pos[i] - 1 - MISS_BUFFER;
            if (end >= start) {
                if (start < 0) {
                    start = 0;
                }
                if (end >= chrom_len) {
                    end = chrom_len - 1;
                }

                for (int32_t j = start; j <= end; ++j) {
                    if (coverage[j] >= 20) {
                        ++cm->t_non_target_good_hits;
                        while (j <= end && coverage[j] > 0) {
                            ++j;
                        }
                    }
                }
            }
            start = tic->end_pos[i] + 1 + MISS_BUFFER;
        }
        end = chrom_len;

        for (int32_t j = start; j <= end; ++j) {
            if (coverage[j] >= 20) {
                ++cm->t_non_target_good_hits;
                while (j <= end && coverage[j] > 0) {
                    ++j;
                }
            }
        }
    }
}

/**
 * Set coverage values for regions of known N bases in reference to 0.
 * Regions are defined as targets in cov_mask_ti.
 */
void handle_coverage_mask(uint32_t *coverage, bed_t *cov_mask_ti,
                          int32_t chrom_idx, int32_t chrom_len)
{
    if (cov_mask_ti->num_targets > 0) {
        bed_chrom_t *tic = cov_mask_ti->chroms[chrom_idx];

        /* for each target set coverage[start:end] values to 0 */
        for (size_t i = 0; i < tic->num_targets; ++i) {
            clear_coverage(coverage, tic->start_pos[i], tic->end_pos[i], chrom_len);
        }
    }
}

/**
 * Erase target regions overlapping masked regions
 */
/*void handle_coverage_mask_target(uint8_t *target_cov, capture_metrics_t *cm,*/
void handle_coverage_mask_target(uint32_t *coverage, capture_metrics_t *cm,
                                 bed_t *cov_mask_ti, int32_t chrom_idx,
                                 int32_t chrom_len)
{
    int32_t start, end;
    uint32_t *curr_pos, *end_pos;
    bed_chrom_t *tic;

    if (cov_mask_ti->num_targets > 0) {
        tic = cov_mask_ti->chroms[chrom_idx];

        /* for each target */
        for (size_t i = 0; i < tic->num_targets; ++i) {
            start = tic->start_pos[i];
            end = tic->end_pos[i];

            for (curr_pos = coverage + start, end_pos = coverage + end;
                 curr_pos <= end_pos;
                 ++curr_pos)
            {
                switch (get_target(*curr_pos)) {
                case TARGET_IN:
                    --cm->b_targeted;
                    break;
                case TARGET_BUFFER:
                    --cm->b_buffer;
                    break;
                default:
                    break;
                }
            }
        }
    }
}

/**
 * Set target positions in target_cov.
 */
void set_target_cov(uint32_t *coverage, capture_metrics_t *cm, bed_t *ti,
                    int32_t chrom_idx, int32_t chrom_len)
{
    uint32_t *curr_pos, *end_pos;
    int32_t start, end, start_, end_;
    bed_chrom_t *tic;

    if (ti->num_targets > 0) {
        tic = ti->chroms[chrom_idx];

        /* for each target */
        for (size_t j = 0; j < tic->num_targets; ++j) {
            start = tic->start_pos[j];
            end = tic->end_pos[j];

            /* Buffers (and Target) */
            start_ = (start > BUFFER) ? start - BUFFER : 0;
            end_ = (end + 1 + BUFFER < chrom_len) ? end + 1 + BUFFER : chrom_len;
            for (size_t pos = start_; pos <= end_; ++pos) {
                coverage[pos] |= TARGET_BUFFER << TARGET_SHIFT;
            }
        }

        for (size_t j = 0; j < tic->num_targets; ++j) {
            start = tic->start_pos[j];
            end = tic->end_pos[j];

            /* Target */
            start_ = (start > 0) ? start : 0;
            end_ = (end < chrom_len) ? end + 1 : chrom_len;
            for (size_t pos = start_; pos <= end_; ++pos) {
                coverage[pos] |= TARGET_IN << TARGET_SHIFT;
            }
        }

        /* Record number of buffer and targeted bases in target_cov */
        for (curr_pos = coverage, end_pos = coverage + chrom_len;
             curr_pos < end_pos;
             ++curr_pos)
        {
            switch (get_target(*curr_pos)) {
            case TARGET_IN:
                ++cm->b_targeted;
                break;
            case TARGET_BUFFER:
                ++cm->b_buffer;
                break;
            default:
                break;
            }
        }
    }
}

void _capture_process_record1(bam1_t *rec, capture_metrics_t *cm)
{
    ++cm->r_aligned;

    /* if paired read */
    if (rec->core.flag & BAM_FPAIRED) {
        ++cm->r_paired;

        /* if mate is mapped */
        if (!(rec->core.flag & BAM_FMUNMAP)) {
            ++cm->r_paired_w_mate;
        }
    }

    /* if duplicate read */
    if (rec->core.flag & BAM_FDUP) {
        ++cm->r_dup;
    }
}

/**
 * Process coverage metrics specific to whole genome or capture.
 */
void _capture_process_record2(bam1_t *rec, capture_metrics_t *cm,
                              target_state_t target_status)
{
    cm->b_aligned += rec->core.l_qseq;

    /* Record read as out of target, in target buffer, or in target */
    switch (target_status) {
    case TARGET_IN:
        ++cm->r_in_target;
        break;
    case TARGET_BUFFER:
        ++cm->r_in_buffer;
        break;
    default:
        ++cm->r_out_target;
        break;
    }
}

/**
 * Process record for target and coverage info.
 * cm_wgs: Capture metrics for whole genome statistics
 * cm_cap: Capture metrics for capture statistics
 */
void capture_process_record(bam1_t *rec, uint32_t *coverage,
                            capture_metrics_t *cm_wgs,
                            capture_metrics_t *cm_cap, int32_t chrom_len,
                            bool remove_dups)
{
    bool in_target, in_buffer;
    uint8_t *qual, *bqual;
    int32_t pos, start, start_pos, end_pos, ref_pos;
    uint32_t target, oplen, *cigar;
    const uint16_t FILTER = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL;
    const uint32_t COV_MAX = (UINT32_MAX - 1) >> 2;


    if (cm_wgs != NULL) {
        ++cm_wgs->r_total;
    }
    if (cm_cap != NULL) {
        ++cm_cap->r_total;
    }

    /* Filter out reads with filtered flags */
    if (rec->core.flag & FILTER) {
        return;
    }

    if (cm_wgs != NULL) {
        _capture_process_record1(rec, cm_wgs);
    }
    if (cm_cap != NULL) {
        _capture_process_record1(rec, cm_cap);
    }

    /* remove duplicate reads */
    if (remove_dups && (rec->core.flag & BAM_FDUP)) {
        return;
    }

    /*
     * Record coverage and whether aligned bases hit target or buffer
     * Confirmed this now matches up with the correct ref pos and # CIGAR Ms
     */
    in_target = false;
    in_buffer = false;
    ref_pos = 0;
    start = rec->core.pos;
    cigar = bam_get_cigar(rec);
    bqual = bam_get_qual(rec);

    /* for each CIGAR op */
    for (uint32_t i = 0; i < rec->core.n_cigar; ++i) {
        oplen = bam_cigar_oplen(cigar[i]);

        switch (bam_cigar_op(cigar[i])) {
        case BAM_CHARD_CLIP:
            break;
        /*
         * M, =, X: record coverage for CIGAR M bases (matches + mismatches)
         * or = (matches) and X (mismatches)
         */
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
            start_pos = pos = start + ref_pos;
            end_pos = start_pos + oplen - 1;

            if (start_pos < 0) {
                start_pos = 0;
                pos = 0;
            }
            if (end_pos >= chrom_len) {
                end_pos = chrom_len - 1;
            }

            /* Record coverage */
            qual = bqual + ref_pos;
            while (pos <= end_pos) {
                if (get_coverage(coverage[pos]) < COV_MAX) {
                    ++coverage[pos];
                } else {
                    /* Y'know just in case */
                    log_warning("Coverage of greater than %u detected. "
                                "Coverage statistics may not be accurate.", COV_MAX);
                }
                ++pos;
                ++qual;
            }
            /* pos == end_pos */

            if (coverage != NULL) {
                while (--pos >= start_pos) {
                    target = get_target(coverage[pos]);
                    if (target == TARGET_BUFFER) {
                        in_buffer = true;
                    } else if (target == TARGET_IN) {
                        in_target = true;
                        break;
                    }
                }
            }
            ref_pos += oplen;
            break;
        /* D: advance ref_pos past deletion */
        case BAM_CDEL:
            ref_pos += oplen;
            break;
        default:
            break;
        }
    }

    if (cm_wgs != NULL) {
        _capture_process_record2(rec, cm_wgs, TARGET_OUT);
    }
    if (cm_cap != NULL) {
        _capture_process_record2(rec, cm_cap,
            in_target ? TARGET_IN : in_buffer ? TARGET_BUFFER : TARGET_OUT);
    }
}

/**
 * Write capture metrics to report.
 */
void capture_report(report_t *report, capture_metrics_t *cm, bed_t *ti)
{
    char *key_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(key_buffer);
    char *value_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(value_buffer);

    const char *prefix = (ti == NULL) ? "Wgs" : "Cap";
    size_t prefix_len = strlen(prefix);
    size_t copy_size = REPORT_BUFFER_SIZE - prefix_len;
    char *key_start = key_buffer + prefix_len;

    /* Capture (targets) or whole genome (no targets) percentage denominator */
    uint64_t denominator = (ti == NULL) ? cm->b_total : cm->b_targeted;

    copy_to_buffer(key_buffer, prefix, REPORT_BUFFER_SIZE);

    copy_to_buffer(key_start, "TotalReads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CovDuplicateReads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_dup);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CovDuplicateReadsPct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->r_dup, cm->r_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "AlignedReads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_aligned);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "AlignedReadsPct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->r_aligned, cm->r_total);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "ReadsPaired", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_paired);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "ReadsPairedWithMates", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_paired_w_mate);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageMean", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%f",
             (denominator != 0) ? (double)cm->c_total / (double)denominator
                                : 0.0);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageMedian", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->c_median);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageStandardDeviation", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%f", cm->c_std_dev);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "ExpectedAlignedReads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_aligned);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CalculatedAlignedReads", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_in_target + cm->r_in_buffer + cm->r_out_target);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases1", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_1_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases1Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_1_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases10", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_10_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases10Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_10_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases20", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_20_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases20Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_20_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases30", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_30_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases30Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_30_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases40", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_40_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases40Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_40_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases50", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_50_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases50Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_50_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases100", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_100_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases100Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_100_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases1000", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_1000_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases1000Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_1000_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    /* ti != NULL? capture stats: wgs */
    if (ti != NULL) {
        copy_to_buffer(key_start, "BufferAlignedReads", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_in_buffer);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "BufferAlignedReadsPct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->r_in_buffer, cm->r_aligned);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "TargetAlignedReads", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_in_target);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "TargetAlignedReadsPct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->r_in_target, cm->r_aligned);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "TargetsHit", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->t_hit);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "TargetsHitPct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->t_hit, cm->t_total);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "TargetBuffersHit", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->t_buffers_hit);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "TargetBuffersHitPct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->t_buffers_hit, cm->t_total);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "TotalTargets", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->t_total);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "HighCoverageNonTargetHits", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->t_non_target_good_hits);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "BasesOnTarget", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_targeted);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "BasesOnBuffer", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_buffer);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "ReadsOnTargetOrBuffer", copy_size);
        snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->r_in_target + cm->r_in_buffer);
        report_add_key_value(report, key_buffer, value_buffer);

        copy_to_buffer(key_start, "ReadsOnTargetOrBufferPct", copy_size);
        print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->r_in_target + cm->r_in_buffer, cm->r_aligned);
        report_add_key_value(report, key_buffer, value_buffer);
    }

    free(key_buffer);
    free(value_buffer);
}
