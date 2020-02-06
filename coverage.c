#include "coverage.h"
#include "err.h"
#include "logging.h"
#include "print.h"
#include <ctype.h>
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

/* Read pair overlapping pair structure */

overlap_handler_t *overlap_handler_init()
{
    overlap_handler_t *olh = calloc(1, sizeof(overlap_handler_t));
    die_on_alloc_fail(olh);

    overlap_buffer_init(olh);
    mc_buffer_init(olh);

    return olh;
}

void overlap_handler_add(overlap_handler_t *olh, int32_t start, int32_t end)
{
    /* if need to enlarge arr */
    if (olh->n_overlap_pairs == olh->overlap_buffer_len) {
        /* double size and realloc */
        uint32_t new_size = 2 * olh->overlap_buffer_len;
        overlap_pair_t *tmp_overlap_buffer = realloc(
            olh->overlap_buffer,
            new_size * sizeof(overlap_pair_t)
        );

        /* check realloc */
        if (tmp_overlap_buffer != NULL) {
            olh->overlap_buffer = tmp_overlap_buffer;
            olh->overlap_buffer_len = new_size;
        } else {
            free(olh->overlap_buffer);
            die_on_alloc_fail(tmp_overlap_buffer);

            return;
        }
    }

    /* set start and end of new pair */
    overlap_pair_t *pair = olh->overlap_buffer + olh->n_overlap_pairs;
    pair->start = start;
    pair->end = end;
    ++olh->n_overlap_pairs;
}

void overlap_handler_clear(overlap_handler_t *olh)
{
    memset(olh->overlap_buffer, 0, olh->overlap_buffer_len * sizeof(overlap_pair_t));
    olh->n_overlap_pairs = 0;
}

void _advance_to_coverage_cigar(uint32_t **cigar_ptr, uint32_t **cigar_end_ptr, int32_t *rpos_ptr)
{
    uint8_t cigar_op;

    while (*cigar_ptr != *cigar_end_ptr) {
        cigar_op = bam_cigar_op(**cigar_ptr);

        if (cigar_op == BAM_CMATCH ||
            cigar_op == BAM_CEQUAL ||
            cigar_op == BAM_CDIFF)
        {
            break;
        } else {
            if (bam_cigar_type(cigar_op) & 0x02) {
                *rpos_ptr += bam_cigar_oplen(**cigar_ptr);
            }

            ++(*cigar_ptr);
        }
    }
}

/* Generate overlap intervals for a read after loading mc_buffer with mc_buffer_load() */
void overlap_handler_generate(overlap_handler_t *olh, bam1_t *rec)
{
    /* ptrs to current cigar ops */
    uint32_t *read_cigar, *mate_cigar, *read_cigar_end, *mate_cigar_end;
    int32_t read_rpos, mate_rpos;

    /* if both read and mate have non-zero # of cigar ops */
    if (rec->core.n_cigar > 0 && olh->n_mc_cigar) {
        /* read and mate pos */
        read_rpos = rec->core.pos;
        mate_rpos = rec->core.mpos;

        /* set read and mate start and end */
        read_cigar = bam_get_cigar(rec);
        read_cigar_end = read_cigar + rec->core.n_cigar;
        mate_cigar = olh->mc_buffer;
        mate_cigar_end = olh->mc_buffer + olh->n_mc_cigar;

        /* advance to first coverage cigar */
        _advance_to_coverage_cigar(&read_cigar, &read_cigar_end, &read_rpos);
        _advance_to_coverage_cigar(&mate_cigar, &mate_cigar_end, &mate_rpos);

        int32_t read_endrpos, mate_endrpos, overlap_start, overlap_end;

        while (read_cigar != read_cigar_end && 
               mate_cigar != mate_cigar_end)
        {
            read_endrpos = read_rpos + bam_cigar_oplen(*read_cigar);
            mate_endrpos = mate_rpos + bam_cigar_oplen(*mate_cigar);

            /* find overlap, add to overlap buffer, and advance cigars */
            if (read_rpos < mate_rpos) {
                /* record overlap */
                if (mate_rpos < read_endrpos) {
                    overlap_start = mate_rpos;
                    overlap_end = (read_endrpos < mate_endrpos)? read_endrpos: mate_endrpos;
                    overlap_handler_add(olh, overlap_start, overlap_end);
                }

                /* advance read_cigar */
                if (bam_cigar_type(*read_cigar) & 0x02) {
                    read_rpos += bam_cigar_oplen(*read_cigar);
                }

                ++read_cigar;
                _advance_to_coverage_cigar(&read_cigar, &read_cigar_end, &read_rpos);
            } else {
                /* record overlap */
                if (read_rpos < mate_endrpos) {
                    overlap_start = read_rpos;
                    overlap_end = (mate_endrpos < read_endrpos)? mate_endrpos: read_endrpos;
                    overlap_handler_add(olh, overlap_start, overlap_end);
                }

                /* advance mate_cigar */
                if (bam_cigar_type(*mate_cigar) & 0x02) {
                    mate_rpos += bam_cigar_oplen(*mate_cigar);
                }

                ++mate_cigar;
                _advance_to_coverage_cigar(&mate_cigar, &mate_cigar_end, &mate_rpos);
            }
        }
    }
}

void overlap_handler_destroy(overlap_handler_t *olh)
{
    mc_buffer_destroy(olh);
    overlap_buffer_destroy(olh);
    free(olh);
}

/* Overlap buffer structure */

void overlap_buffer_init(overlap_handler_t *olh)
{
    int32_t initial_len = 32;

    olh->overlap_buffer = calloc(initial_len, sizeof(overlap_pair_t));
    die_on_alloc_fail(olh->overlap_buffer);
    olh->overlap_buffer_len = initial_len;
    olh->n_overlap_pairs = 0;
}

void overlap_buffer_destroy(overlap_handler_t *olh)
{
    if (olh != NULL) {
        free(olh->overlap_buffer);
        olh->overlap_buffer = NULL;
        olh->overlap_buffer_len = 0;
        olh->n_overlap_pairs = 0;
    }
}

/* Mate CIGAR buffer structure */

void mc_buffer_init(overlap_handler_t *olh)
{
    uint32_t initial_len = 32;

    olh->mc_buffer = calloc(initial_len, sizeof(uint32_t));
    die_on_alloc_fail(olh->mc_buffer);

    olh->mc_buffer_len = initial_len;
    olh->n_mc_cigar = 0;
}

void mc_buffer_load(overlap_handler_t *olh, uint32_t n_mc_cigar, char *mc_str)
{
    bool valid_cigar;
    char *c, *end, cigar_opchr;
    uint8_t cigar_op;
    uint32_t cigar_len, cigar_idx;

    cigar_idx = 0;

    if (mc_str != NULL && n_mc_cigar > 0) {
        /* if need to enlarge arr */
        if (n_mc_cigar > olh->mc_buffer_len) {
            /* size = n_mc_cigar + 4 and realloc */
            uint32_t new_size = n_mc_cigar + 4;
            uint32_t *tmp_mc_buffer = realloc(olh->mc_buffer, new_size * sizeof(uint32_t));

            /* check realloc */
            if (tmp_mc_buffer != NULL) {
                olh->mc_buffer = tmp_mc_buffer;
                olh->mc_buffer_len = new_size;
            } else {
                free(olh->mc_buffer);
                die_on_alloc_fail(tmp_mc_buffer);
            }
        }

        c = mc_str;
        cigar_len = 0;
        cigar_opchr = '\0';
        cigar_op = 0;
        valid_cigar = true;

        while (*c != '\0' && cigar_idx < n_mc_cigar) {
            /* parse cigar len */
            cigar_len = strtol(c, &end, 10);

            /* parse cigar opchr */
            if (end != c && *end != '\0') {
                cigar_opchr = *end;
                ++end;
                c = end;
            } else {
                break;
            }

            /* cigar op from opchr */
            switch (cigar_opchr) {
            case 'M': cigar_op = 0x00; break;
            case 'I': cigar_op = 0x01; break;
            case 'D': cigar_op = 0x02; break;
            case 'N': cigar_op = 0x03; break;
            case 'S': cigar_op = 0x04; break;
            case 'H': cigar_op = 0x05; break;
            case 'P': cigar_op = 0x06; break;
            case '=': cigar_op = 0x07; break;
            case 'X': cigar_op = 0x08; break;
            case 'B': cigar_op = 0x09; break;
            default: valid_cigar = false; break;
            }

            if (valid_cigar) {
                olh->mc_buffer[cigar_idx] = (cigar_len << 4) | cigar_op;
                ++cigar_idx;
            } else {
                break;
            }
        }
    }

    /*return cigar_idx;*/
    olh->n_mc_cigar = cigar_idx;
}

void mc_buffer_destroy(overlap_handler_t *olh)
{
    if (olh != NULL) {
        free(olh->mc_buffer);
        olh->mc_buffer = NULL;
        olh->n_mc_cigar = 0;
    }
}

uint32_t _count_cigar_ops(char *cigar_str)
{
    char *c = cigar_str;
    uint32_t count = 0;

    while (*c != '\0') {
        if (isalpha(*c)) {
            ++count;
        }
        ++c;
    }

    return count;
}

void process_record_overlap(overlap_handler_t *olh, bam1_t *rec, bool remove_overlaps, bool remove_overlaps_mc)
{
    char *mc_str;
    uint8_t *mc_aux;
    int32_t r_end_pos;

    /* clear out overlap buffer */
    overlap_handler_clear(olh);

    /* position-based overlap removal */
    if (remove_overlaps) {
        int32_t overlap_start, overlap_end;
        r_end_pos = bam_endpos(rec);

        if (!(rec->core.flag & BAM_FMUNMAP) && /* mate pair is mapped */
            rec->core.tid == rec->core.mtid && /* mate pair is on same chromosome */
            rec->core.pos < rec->core.mpos &&  /* read is first read in coord-sort order */
            rec->core.mpos < r_end_pos)        /* mate pair overlaps seq before ref-aligned end pos */
        {
            /* add single overlap pair (mate pos, read end pos) */
            overlap_start = rec->core.mpos;
            overlap_end = bam_endpos(rec);
            overlap_handler_add(olh, overlap_start, overlap_end);
        }
    /* MC tag-based overlap removal */
    } else if (remove_overlaps_mc) {
        r_end_pos = bam_endpos(rec);

        if (!(rec->core.flag & BAM_FMUNMAP) && /* mate pair is mapped */
            rec->core.tid == rec->core.mtid && /* mate pair is on same chromosome */
            rec->core.pos < rec->core.mpos &&  /* read is first read in coord-sort order */
            rec->core.mpos < r_end_pos)        /* mate pair overlaps seq before ref-aligned end pos */
            /*rec->core.flag & BAM_FREAD2)       / * read is read 2 */
        {
            /* if MC tag is present */
            if ((mc_aux = bam_aux_get(rec, "MC")) != NULL) {
                mc_str = bam_aux2Z(mc_aux);

                /* load mc buffer from MC tag */
                mc_buffer_load(olh, _count_cigar_ops(mc_str), mc_str);

                /* regenerate overlap pairs */
                overlap_handler_generate(olh, rec);
            }
        }
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
    uint64_t k = 0, median_sum = 0, sum = 0;
    uint64_t mid = (ti == NULL ? cm->b_total : cm->b_targeted) / 2;
    bool median_set = false;

    /* Set median coverage */
    for (size_t i = 0; i < ci->cov_histo_len; ++i) {
        if (ci->cov_histo[i] > 0) {
            k += ci->cov_histo[i];
            median_sum += ci->cov_histo[i];
            if (!median_set && median_sum >= mid) {
                cm->c_median = i;
                median_set = true;
            }
            sum += i * ci->cov_histo[i];
        }
    }

    cm->c_std_dev = k > 1
        ? sqrt(((double)cm->c_sum_sq - (double)sum * ((double)sum / (double)k)) /
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
    uint64_t *tmp_cov_histo, prev_histo_len;

    if ((size_t)cov >= ci->cov_histo_len) {
        /* Buffer an additional 256 to reduce # of reallocs for slowly increasing coverage */
        prev_histo_len = ci->cov_histo_len;
        ci->cov_histo_len = (size_t)(cov + 1 + 256);
        tmp_cov_histo = realloc(ci->cov_histo, ci->cov_histo_len * sizeof(uint64_t));

        /* check realloc */
        if (tmp_cov_histo != NULL) {
            ci->cov_histo = tmp_cov_histo;
            for (size_t i = (size_t)prev_histo_len; i < ci->cov_histo_len; ++i) {
                ci->cov_histo[i] = 0;
            }
        } else {
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
            if (cov == 0) {
                goto cov0;
            } else {
                if (cov < 15) {
                    if (cov < 10) {
                        goto cov1;
                    } else {
                        goto cov10;
                    }
                } else {
                    if (cov < 20) {
                        goto cov15;
                    } else {
                        goto cov20;
                    }
                }
            }
        } else {
            if (cov < 70) {
                if (cov < 50) {
                    if (cov < 40) {
                        goto cov30;
                    } else {
                        goto cov40;
                    }
                } else {
                    if (cov < 60) {
                        goto cov50;
                    } else {
                        goto cov60;
                    }
                }
            } else {
                if (cov < 500) {
                    if (cov < 100) {
                        goto cov70;
                    } else {
                        goto cov100;
                    }
                } else {
                    if (cov < 1000) {
                        goto cov500;
                    } else {
                        goto cov1000;
                    }
                }
            }
        }

cov1000:++cm->b_1000_plus_hits;
cov500: ++cm->b_500_plus_hits;
cov100: ++cm->b_100_plus_hits;
cov70:  ++cm->b_70_plus_hits;
cov60:  ++cm->b_60_plus_hits;
cov50:  ++cm->b_50_plus_hits;
cov40:  ++cm->b_40_plus_hits;
cov30:  ++cm->b_30_plus_hits;
cov20:  ++cm->b_20_plus_hits;
cov15:  ++cm->b_15_plus_hits;
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
                if (cov == 0) {
                    goto tgtcov0;
                } else {
                    if (cov < 15) {
                        if (cov < 10) {
                            goto tgtcov1;
                        } else {
                            goto tgtcov10;
                        }
                    } else {
                        if (cov < 20) {
                            goto tgtcov15;
                        } else {
                            goto tgtcov20;
                        }
                    }
                }
            } else {
                if (cov < 70) {
                    if (cov < 50) {
                        if (cov < 40) {
                            goto tgtcov30;
                        } else {
                            goto tgtcov40;
                        }
                    } else {
                        if (cov < 60) {
                            goto tgtcov50;
                        } else {
                            goto tgtcov60;
                        }
                    }
                } else {
                    if (cov < 500) {
                        if (cov < 100) {
                            goto tgtcov70;
                        } else {
                            goto tgtcov100;
                        }
                    } else {
                        if (cov < 1000) {
                            goto tgtcov500;
                        } else {
                            goto tgtcov1000;
                        }
                    }
                }
            }

tgtcov1000: ++cm->b_1000_plus_hits;
tgtcov500:  ++cm->b_500_plus_hits;
tgtcov100:  ++cm->b_100_plus_hits;
tgtcov70:   ++cm->b_70_plus_hits;
tgtcov60:   ++cm->b_60_plus_hits;
tgtcov50:   ++cm->b_50_plus_hits;
tgtcov40:   ++cm->b_40_plus_hits;
tgtcov30:   ++cm->b_30_plus_hits;
tgtcov20:   ++cm->b_20_plus_hits;
tgtcov15:   ++cm->b_15_plus_hits;
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
            end_ = (end < chrom_len) ? end : chrom_len;
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
                            capture_metrics_t *cm_wgs, capture_metrics_t *cm_cap,
                            overlap_handler_t *olh, int32_t chrom_len,
                            bool remove_dups, bool remove_overlaps, bool remove_overlaps_mc,
                            uint8_t min_base_qual)
{
    bool in_target, in_buffer, check_overlap;
    uint8_t *qual, *bqual;
    int32_t pos, start, start_pos, end_pos, ref_pos;
    uint32_t target, oplen, *cigar;
    uint64_t b_filt_overlap, b_filt_basequal;
    uint64_t b_filt_target_overlap, b_filt_target_basequal;
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

    /* process record for possible overlap with mate pair */
    process_record_overlap(olh, rec, remove_overlaps, remove_overlaps_mc);

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
    b_filt_overlap = 0;
    b_filt_basequal = 0;
    b_filt_target_overlap = 0;
    b_filt_target_basequal = 0;

    overlap_pair_t *curr_overlap_pair, *overlap_pair_end;

    check_overlap = (olh->n_overlap_pairs > 0);

    if (check_overlap) {
        curr_overlap_pair = olh->overlap_buffer;
        overlap_pair_end = curr_overlap_pair + olh->n_overlap_pairs;
    } else {
        curr_overlap_pair = NULL;
        overlap_pair_end = NULL;
    }

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
            /* Record coverage */
            if (coverage != NULL) {
                start_pos = pos = start + ref_pos;
                end_pos = start_pos + oplen - 1;

                if (start_pos < 0) {
                    start_pos = 0;
                    pos = 0;
                }
                if (end_pos >= chrom_len) {
                    end_pos = chrom_len - 1;
                }

                qual = bqual + ref_pos;

                while (pos <= end_pos) {
                    if (get_coverage(coverage[pos]) < COV_MAX) {
                        if (check_overlap &&
                            pos >= curr_overlap_pair->start &&
                            pos < curr_overlap_pair->end)
                        {
                            ++b_filt_overlap;
                            if (get_target(coverage[pos]) == TARGET_IN) {
                                ++b_filt_target_overlap;
                            }
                        } else if (*qual < min_base_qual) {
                            ++b_filt_basequal;
                            if (get_target(coverage[pos]) == TARGET_IN) {
                                ++b_filt_target_overlap;
                            }
                        } else {
                            ++coverage[pos];
                        }
                    } else {
                        /* Y'know just in case */
                        log_warning("Coverage of greater than %u detected. "
                                    "Coverage statistics may not be accurate.", COV_MAX);
                    }

                    ++pos;
                    ++qual;

                    /* if pos is past current overlap pair, advance to next */
                    if (check_overlap && pos >= curr_overlap_pair->end) {
                        ++curr_overlap_pair;

                        /* no more overlap pairs */
                        if (curr_overlap_pair == overlap_pair_end) {
                            curr_overlap_pair = NULL;
                            overlap_pair_end = NULL;
                            check_overlap = false;
                        }
                    }
                }

                /* pos == end_pos */
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
        cm_wgs->b_filt_overlap += b_filt_overlap;
        cm_wgs->b_filt_basequal += b_filt_basequal;
        _capture_process_record2(rec, cm_wgs, TARGET_OUT);
    }
    if (cm_cap != NULL) {
        cm_cap->b_filt_overlap += b_filt_target_overlap;
        cm_cap->b_filt_basequal += b_filt_target_basequal;
        _capture_process_record2(rec, cm_cap, in_target ? TARGET_IN : in_buffer ? TARGET_BUFFER : TARGET_OUT);
    }
}

/**
 * Write capture metrics to report.
 */
void capture_report(report_t *report, capture_metrics_t *cm, bed_t *ti, char *key_buffer, char *value_buffer)
{
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

    copy_to_buffer(key_start, "CoverageBases15", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_15_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases15Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_15_plus_hits, denominator);
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

    copy_to_buffer(key_start, "CoverageBases60", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_60_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases60Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_60_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases70", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_70_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases70Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_70_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases100", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_100_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases100Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_100_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases500", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_500_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases500Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_500_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases1000", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_1000_plus_hits);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "CoverageBases1000Pct", copy_size);
    print_pct(value_buffer, REPORT_BUFFER_SIZE, cm->b_1000_plus_hits, denominator);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "FilteredOverlapBases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_filt_overlap);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_start, "FilteredLowBaseQualityBases", copy_size);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", cm->b_filt_basequal);
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
}
