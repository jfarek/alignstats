# AlignStats

## Overview

AlignStats produces various alignment, whole genome coverage, and capture
coverage metrics for sequence alignment files in SAM, BAM, and CRAM format.
This program is designed to serve reporting and quality control purposes in
sequencing analysis pipelines at the Baylor College of Medicine Human Genome
Sequencing Center (BCM-HGSC). Descriptions of reported metrics can be found in
``doc/``.

## Building

AlignStats requires a recent version of HTSlib (versions 1.3 and higher are
known to work, version 1.4 or higher is suggested). The ``-p`` option for
multithreading requires libpthread support.

To build AlignStats, run:

    autoconf
    ./configure
    make

AlignStats multithreading support is enabled by running ``configure`` as:

    ./configure --enable-multithreading

## Usage

    alignstats [-i INPUT] [-j FORMAT] [-o OUTPUT]
               [-h] [-v] [-n NUMREADS] [-p] [-P INT]
               [-r REGIONS] [-t TARGET] [-m COVMASK] [-T REFFASTA]
               [-q INT] [-f INT] [-F INT]
               [-D] [-U] [-A] [-C] [-W]

    Runtime options:
        -h          Print usage information.
        -v          Print verbose runtime information output to stderr.
        -n INT      Maximum number of records to keep in memory.
        -p          Use separate threads for reading and processing records
                    (requires builtin pthread support).
        -P INT      Number of HTSlib decompression threads to spawn.

    File options:
        -i INPUT    Read INPUT as the input SAM, BAM, or CRAM file (stdin). Input
                    must be coordinate-sorted for accurate results.
        -j FORMAT   Specify file format of input alignment file ("sam", "bam", or
                    "cram" available, default guessed from filename or "sam").
        -o OUTPUT   Write report to OUTPUT (stdout).
        -r REGIONS  File in BED format listing which regions to process. By
                    default, all available records are processed. This option
                    requires the alignment file to be indexed.
        -t TARGET   File in BED format listing capture coverage regions. Required
                    if capture coverage statistics are enabled.
        -m COVMASK  File in BED format listing regions of N bases in reference.
                    Coverage counts will be suppressed for these regions.
        -T REFFASTA Indexed FASTA reference file for CRAM input alignment.

    Processing options:
        -q INT      Only process records with minimum read quality of INT.
        -f INT      Only process records with all bits in INT set in FLAG.
        -F INT      Only process records with none of bits in INT set in FLAG.

    Reporting options:
        -D          Disable excluding duplicate reads from coverage statistics.
        -U          Disable processing unplaced unmapped reads (CHROM "*") when
                    using the -r option.
        -A          Disable reporting alignment statistics.
        -C          Disable reporting capture coverage statistics.
        -W          Disable reporting whole genome coverage statistics.

### Usage Examples

    alignstats -i sample.bam -o report.txt -r regions.bed -t targets.bed -m n_regions.bed

Run AlignStats with ``sample.bam`` as the input alignment file,
``regions.bed`` as the input regions file, ``targets.bed`` as the input target
file, ``n_regions.bed`` as the input coverage mask file, and ``report.txt`` as
the output report file. Both whole genome and capture statistics will be
produced for coverage metrics.

    alignstats -v -p -i sample.bam -r regions.bed -t targets.bed -m n_regions.bed -o report.txt

Run the first command in verbose mode and enable multithreading (if supported).

    alignstats -i sample.bam -r regions.bed -t targets.bed -m n_regions.bed

Same as the first command, but the report will be printed to stdout.

    alignstats -C -i sample.bam -o report.txt

Only produce alignment statistics and whole genome coverage statistics
(suppress capture coverage statistics) for ``sample.bam``.

## Disclaimers

While extensive testing and validation have been performed at the BCM-HGSC to
confirm the accuracy of metrics reported by AlignStats, no guarantees of
accuracy are made or implied.

The set of metrics reported by AlignStats may change in the course of continued
development by the BCM-HGSC.
