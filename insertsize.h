#ifndef _INSERTSIZE_H
#define _INSERTSIZE_H

#include "htslib/sam.h"
#include "report.h"
#include "treemap.h"
#include <stdint.h>

/* Insert size metrics structure */
struct insert_size_metrics {
    double mean;
    double std_dev;
    uint64_t median;
    uint64_t mode;
    uint64_t sum_sq;
    uint16_t filter_incl;
    uint16_t filter_excl;
    tree_map_t *insert_size_map;
};
typedef struct insert_size_metrics insert_size_metrics_t;

insert_size_metrics_t *insert_size_metrics_init();
void insert_size_metrics_destroy(insert_size_metrics_t *ism);

/* Insert size metrics calculation */
void insert_size_process_record(bam1_t *rec, insert_size_metrics_t *ism);
void insert_size_finalize(insert_size_metrics_t *ism);
void insert_size_report(report_t *report, insert_size_metrics_t *ism, char *key_buffer, char *value_buffer);

#endif /* _INSERTSIZE_H */
