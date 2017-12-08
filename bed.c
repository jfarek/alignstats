#include "bed.h"
#include "err.h"
#include "logging.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * File structure for reading files in subset of BED format (first 3 fields)
 */

/**
 * Create and return new *bed_chrom_t.
 */
bed_chrom_t *bed_chrom_init(size_t count)
{
    size_t size;

    bed_chrom_t *tic = calloc(1, sizeof(bed_chrom_t));
    die_on_alloc_fail(tic);

    tic->start_pos = NULL;
    tic->end_pos = NULL;

    if ((tic->num_targets = count) > 0) {
        size = count * sizeof(uint32_t);
        tic->start_pos = malloc(size);
        die_on_alloc_fail(tic->start_pos);
        tic->end_pos = malloc(size);
        die_on_alloc_fail(tic->end_pos);
    }

    return tic;
}

/**
 * Free *bed_chrom_t from memory.
 */
void bed_chrom_destroy(bed_chrom_t *tic)
{
    if (tic != NULL) {
        free(tic->start_pos);
        free(tic->end_pos);
        free(tic);
    }
}

/**
 * Create and return new *bed_t.
 */
bed_t *bed_init()
{
    bed_t *ti = calloc(1, sizeof(bed_t));
    die_on_alloc_fail(ti);

    ti->chroms = NULL;

    return ti;
}

/**
 * Returns sum of bed region lengths.
 */
uint32_t bed_sum_bases(bed_t *ti)
{
    uint32_t sum = 0;

    if (ti != NULL) {
        sum = (uint32_t)ti->num_chroms;
        for (size_t i = 0; i < ti->num_chroms; ++i) {
            bed_chrom_t *tic = ti->chroms[i];
            for (size_t j = 0; j < tic->num_targets; ++j) {
                sum += tic->end_pos[j] - tic->start_pos[j];
            }
        }
    }

    return sum;
}

/**
 * Free *bed_t from memory.
 */
void bed_destroy(bed_t *ti)
{
    if (ti != NULL) {
        for (size_t i = 0; i < ti->num_chroms; ++i) {
            bed_chrom_destroy(ti->chroms[i]);
        }
        free(ti->chroms);
        free(ti);
    }
}

/**
 * Get index of chrom_buffer (chromosome name) in chrom_names of length
 * n_targets.
 */
int _get_chrom_idx(char **chrom_names, char *chrom_buffer, int32_t n_targets)
{
    int ret = -9;

    for (int32_t i = 0; i < n_targets && ret == -9; ++i) {
        if (strcmp(chrom_buffer, chrom_names[i]) == 0) {
            ret = i;
        }
    }

    return ret;
}

/**
 * Parse line from fp and set chrom_buffer, start and end.
 * Returns whether parsing was successful.
 */
bool _parse_bed(FILE *fp, char *line_buffer, char *chrom_buffer, int *start, int *end)
{
    bool ret = true;
    int num_read_chrom;
    char *c_start, *c_end;
    size_t strlen_browser, strlen_track;

    strlen_browser = strlen("browser");
    strlen_track = strlen("track");

    if (fgets(line_buffer, CHROM_BUFFER_SIZE, fp) == NULL) {
        ret = false;
    } else {
        num_read_chrom = sscanf(line_buffer, CHROM_BUFFER_SIZE_SCAN, chrom_buffer);

        /* Ignore header lines and improperly formatted lines. */
        if (num_read_chrom < 1 ||
            *chrom_buffer == '#' ||
            strncmp(chrom_buffer, "browser", strlen_browser) == 0 ||
            strncmp(chrom_buffer, "track", strlen_track) == 0)
        {
            ret = false;
        } else {
            c_start = line_buffer + strlen(chrom_buffer);
            *start = strtol(c_start, &c_end, 10);
            c_start = c_end;
            *end = strtol(c_start, &c_end, 10);
        }
    }

    return ret;
}

/**
 * Load targets file. Header lines and improperly formatted lines are ignored.
 * Returns number of targets loaded.
 */
size_t load_bed(FILE *fp, bed_t *ti, bam_hdr_t *hdr)
{
    char *chrom_buffer, *line_buffer;
    int chrom_idx;
    size_t target_size, counts[hdr->n_targets];
    int32_t start, end, chrom_len;

    ti->num_chroms = hdr->n_targets;
    ti->chroms = malloc(((size_t)hdr->n_targets * sizeof(bed_chrom_t *)));
    die_on_alloc_fail(ti->chroms);

    line_buffer = malloc(CHROM_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(line_buffer);
    chrom_buffer = malloc(CHROM_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(chrom_buffer);

    target_size = hdr->n_targets * sizeof(size_t);
    memset(counts, 0, target_size);

    /* for each line in target file */
    while (!feof(fp)) {
        /* count number of targets per chromosome */
        if (!_parse_bed(fp, line_buffer, chrom_buffer, &start, &end) ||
            (chrom_idx = _get_chrom_idx(hdr->target_name, chrom_buffer, hdr->n_targets)) < 0)
        {
            continue;
        }
        ++counts[chrom_idx];
    }

    for (int32_t i = 0; i < hdr->n_targets; ++i) {
        ti->chroms[i] = bed_chrom_init(counts[i]);
        ti->num_targets += counts[i];
    }

    rewind(fp);
    memset(counts, 0, target_size);

    /* for each line in target file */
    while (!feof(fp)) {
        /* record targets */
        if (!_parse_bed(fp, line_buffer, chrom_buffer, &start, &end) ||
            (chrom_idx = _get_chrom_idx(hdr->target_name, chrom_buffer, hdr->n_targets)) < 0)
        {
            continue;
        }

        if (start < 0) {
            log_warning("BED start coordinate less than minimum allowed: %d < 0", start);
            start = 0;
        }
        chrom_len = (int32_t)hdr->target_len[chrom_idx];
        if (chrom_len > INT32_MAX) {
            chrom_len = INT32_MAX;
        }
        if (end > chrom_len) {
            log_warning("BED end coordinate greater than maximum allowed: %d > %d", end, chrom_len);
            end = chrom_len;
        }

        if (ti->chroms[chrom_idx]->start_pos != NULL) {
            ti->chroms[chrom_idx]->start_pos[counts[chrom_idx]] = start;
        }
        if (ti->chroms[chrom_idx]->end_pos != NULL) {
            ti->chroms[chrom_idx]->end_pos[counts[chrom_idx]] = end - 1;
        }
        ++counts[chrom_idx];
    }

    free(line_buffer);
    free(chrom_buffer);

    return ti->num_targets;
}
