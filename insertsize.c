#include "insertsize.h"
#include "err.h"
#include "logging.h"
#include "print.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* Insert size metrics structure */

/**
 * Create new *insert_size_metrics_t and initialize values
 */
insert_size_metrics_t *insert_size_metrics_init()
{
    insert_size_metrics_t *ism = calloc(1, sizeof(insert_size_metrics_t));
    die_on_alloc_fail(ism);

    ism->insert_size_map = tree_map_init();
    ism->filter_incl = BAM_FPAIRED |
                       BAM_FPROPER_PAIR;
    ism->filter_excl = BAM_FREAD1 |
                       BAM_FSECONDARY |
                       BAM_FSUPPLEMENTARY |
                       BAM_FUNMAP |
                       BAM_FMUNMAP |
                       BAM_FDUP;

    return ism;
}

/**
 * Free *ism from memory.
 */
void insert_size_metrics_destroy(insert_size_metrics_t *ism)
{
    if (ism != NULL) {
        tree_map_destroy(ism->insert_size_map);
        free(ism);
    }
}

/* Insert size metrics calculation */

/**
 * Process a bam1_t record for insert size metrics.
 */
void insert_size_process_record(bam1_t *rec, insert_size_metrics_t *ism)
{
    int32_t isize;
    tree_node_t *node;

    if (((rec->core.flag & ism->filter_incl) == ism->filter_incl) && /* reads to filter in */
        !(rec->core.flag & ism->filter_excl) && /* reads to filter out */
        rec->core.tid == rec->core.mtid)        /* read and mate on same chr */
    {
        isize = abs(rec->core.isize);
        node = tree_map_get(ism->insert_size_map, isize);

        /* Increment insert size in map */
        tree_map_set(ism->insert_size_map, isize, (node == NULL) ? 1 : node->value + 1);

        ism->sum_sq += (uint64_t)isize * (uint64_t)isize;
    }
}

/**
 * Finalize insert size metrics once all records are processed.
 * Calculates mean, median, and mode insert sizes.
 */
void insert_size_finalize(insert_size_metrics_t *ism)
{
    tree_node_key_t *keyset;
    tree_node_t *node;
    uint64_t k, sum, mode_size, curr_size, median_idx;
    bool median_set = false;

    if (tree_map_set_keyset(&keyset, ism->insert_size_map)) {
        k = sum = mode_size = 0;

        /* Mean and mode alignment sizes */
        for (size_t i = 0; i < ism->insert_size_map->num_nodes; ++i) {
            if (!tree_map_set_node(&node, ism->insert_size_map, keyset[i])) {
                goto fail;
            }

            if ((curr_size = (uint64_t)node->value) > mode_size) {
                mode_size = curr_size;
                ism->mode = (uint64_t)keyset[i];
            }
            log_warning("%d\t%d", keyset[i], node->value);

            k += curr_size;
            sum += (uint64_t)keyset[i] * curr_size;
        }

        if (k != 0) {
            ism->mean = (double)sum / (double)k;

            /* Median alignment size */
            median_idx = k / 2;
            k = 0;
            for (size_t i = 0; i < ism->insert_size_map->num_nodes; ++i) {
                if (!tree_map_set_node(&node, ism->insert_size_map, keyset[i])) {
                    goto fail;
                }

                k += (uint64_t)node->value;
                if (!median_set && k >= median_idx) {
                    ism->median = (uint64_t)keyset[i];
                    median_set = true;
                }
            }

            /* Sample standard deviation insertion size */
            ism->std_dev = k > 1
                ? sqrt(((double)ism->sum_sq - ((double)sum * (double)sum) / (double)k) /
                       (double)(k - 1))
                : 0.0;
        } else {
            ism->mean = 0.0;
            ism->median = 0;
            ism->std_dev = 0.0;
        }

fail:
        free(keyset);
    }
}

/**
 * Write metrics in ism to report.
 */
void insert_size_report(report_t *report, insert_size_metrics_t *ism)
{
    char *key_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(key_buffer);
    char *value_buffer = malloc(REPORT_BUFFER_SIZE * sizeof(char));
    die_on_alloc_fail(value_buffer);

    copy_to_buffer(key_buffer, "InsertSizeMean", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%f", ism->mean);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "InsertSizeMedian", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", ism->median);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "InsertSizeMode", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%lu", ism->mode);
    report_add_key_value(report, key_buffer, value_buffer);

    copy_to_buffer(key_buffer, "InsertSizeStandardDeviation", REPORT_BUFFER_SIZE);
    snprintf(value_buffer, REPORT_BUFFER_SIZE, "%f", ism->std_dev);
    report_add_key_value(report, key_buffer, value_buffer);

    free(key_buffer);
    free(value_buffer);
}
