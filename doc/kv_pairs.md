# AlignStats Key-Value Pairs

<style>
table, th, td {
    border-color:black;
    border-style:solid;
    border-width:1px;
}
</style>

<table>
<tr><th>Key</th><th>Type</th><th>Description</th></tr>
<tr><td>TotalRecords</td><td>Integer</td><td>Total number of reads from input</td></tr>
<tr><td>UnfilteredRecords</td><td>Integer</td><td>Total number of unfiltered reads</td></tr>
<tr><td>UnfilteredRecordsPct</td><td>Float</td><td>Number of unfiltered reads as a percentage of TotalReads</td></tr>
<tr><td>FilteredRecords</td><td>Integer</td><td>Total number of filtered reads</td></tr>
<tr><td>FilteredRecordsPct</td><td>Float</td><td>Number of filtered reads as a percentage of TotalReads</td></tr>
<tr><td>YieldReads</td><td>Integer</td><td>Total number of unfiltered reads processed</td></tr>
<tr><td>YieldBases</td><td>Integer</td><td>Total number of unfiltered bases processed</td></tr>
<tr><td>UnmappedReads</td><td>Integer</td><td>Number of unmapped reads</td></tr>
<tr><td>UnmappedReadsPct</td><td>Float</td><td>Unmapped reads as percentage of YieldReads</td></tr>
<tr><td>UnmappedBases</td><td>Integer</td><td>Sum of bases in unmapped reads</td></tr>
<tr><td>UnmappedBasesPct</td><td>Float</td><td>Bases in unmapped reads as percentage of YieldBases</td></tr>
<tr><td>DuplicateReads</td><td>Integer</td><td>Number of mapped duplicate reads</td></tr>
<tr><td>DuplicateReadsPct</td><td>Float</td><td>Mapped duplicate reads as percentage of mapped reads passing QC</td></tr>
<tr><td>DuplicateBases</td><td>Integer</td><td>Sum of bases in duplicate reads</td></tr>
<tr><td>DuplicateBasesPct</td><td>Float</td><td>Sum of bases in duplicate reads as percentage of bases in mapped reads passing QC</td></tr>
<tr><td>MappedReads</td><td>Integer</td><td>Number of mapped reads</td></tr>
<tr><td>MappedReadsPct</td><td>Float</td><td>Mapped reads as percentage of YieldReads</td></tr>
<tr><td>MappedBases</td><td>Integer</td><td>Sum of bases in mapped reads</td></tr>
<tr><td>MappedBasesPct</td><td>Float</td><td>Bases in mapped reads as percentage of YieldReads</td></tr>
<tr><td>AlignedBases</td><td>Integer</td><td>MappedBases minus Soft-Clipped Bases</td></tr>
<tr><td>AlignedBasesPct</td><td>Float</td><td>AlignedBases as percentage of YieldReads</td></tr>
<tr><td>MatchedBases</td><td>Integer</td><td>Number of aligned bases matching reference base</td></tr>
<tr><td>MatchedBasesPct</td><td>Float</td><td>Aligned bases matching reference base as percentage of AlignedBases</td></tr>
<tr><td>MismatchedBases</td><td>Integer</td><td>Number of aligned bases not matching reference base</td></tr>
<tr><td>MismatchedBasesPct</td><td>Float</td><td>Aligned bases not matching reference base as percentage of AlignedBases</td></tr>
<tr><td>InsertedBases</td><td>Integer</td><td>Number of insertion error bases</td></tr>
<tr><td>InsertedBasesPct</td><td>Float</td><td>Insertion error bases as percentage of AlignedBases</td></tr>
<tr><td>DeletedBases</td><td>Integer</td><td>Number of deletion error bases</td></tr>
<tr><td>DeletedBasesPct</td><td>Float</td><td>Deletion error bases as percentage of AlignedBases</td></tr>
<tr><td>SoftClippedReads</td><td>Integer</td><td>Number of reads with soft-clipping</td></tr>
<tr><td>SoftClippedReadsPct</td><td>Float</td><td>Reads with soft-clipping as percentage of MappedReads</td></tr>
<tr><td>SoftClippedBases</td><td>Integer</td><td>Number of soft-clipped error bases</td></tr>
<tr><td>SoftClippedBasesPct</td><td>Float</td><td>Soft-clipped error bases as percentage of MappedBases</td></tr>
<tr><td>PerfectReads</td><td>Integer</td><td>Number of aligned reads without any mismatches, insertions, or deletions</td></tr>
<tr><td>PerfectReadsPct</td><td>Float</td><td>Aligned reads without any mismatches, insertions, or deletions as percentage of MappedReads</td></tr>
<tr><td>PerfectBases</td><td>Integer</td><td>Sum of bases in PerfectBases</td></tr>
<tr><td>PerfectBasesPct</td><td>Float</td><td>Bases in PerfectBases as percentage of MappedBases</td></tr>
<tr><td>Q20Bases</td><td>Integer</td><td>Number of bases with quality scores of at least 20</td></tr>
<tr><td>Q20BasesPct</td><td>Float</td><td>Bases with quality scores of at least 20 as percentage of AlignedBases</td></tr>
<tr><td>Q30Bases</td><td>Integer</td><td>Number of bases with quality scores of at least 30</td></tr>
<tr><td>Q30BasesPct</td><td>Float</td><td>Bases with quality scores of at least 30 as percentage of AlignedBases</td></tr>
<tr><td>R1YieldReads</td><td>Integer</td><td>Total number of reads for R1</td></tr>
<tr><td>R1YieldBases</td><td>Integer</td><td>Total number of bases for R1</td></tr>
<tr><td>R1UnmappedReads</td><td>Integer</td><td>Number of unmapped reads for R1</td></tr>
<tr><td>R1UnmappedReadsPct</td><td>Float</td><td>Unmapped reads as percentage of YieldReads for R1</td></tr>
<tr><td>R1UnmappedBases</td><td>Integer</td><td>Sum of bases in unmapped reads for R1</td></tr>
<tr><td>R1UnmappedBasesPct</td><td>Float</td><td>Bases in unmapped reads as percentage of YieldBases for R1</td></tr>
<tr><td>R1MappedReads</td><td>Integer</td><td>Number of mapped reads for R1</td></tr>
<tr><td>R1MappedReadsPct</td><td>Float</td><td>Mapped reads as percentage of YieldReads for R1</td></tr>
<tr><td>R1MappedBases</td><td>Integer</td><td>Sum of bases in mapped reads for R1</td></tr>
<tr><td>R1MappedBasesPct</td><td>Float</td><td>Sum of bases in mapped reads as percentage of YieldReads for R1</td></tr>
<tr><td>R1AlignedBases</td><td>Integer</td><td>MappedBases minus Soft-Clipped Bases for R1</td></tr>
<tr><td>R1AlignedBasesPct</td><td>Float</td><td>AlignedBases as percentage of YieldReads for R1</td></tr>
<tr><td>R1MatchedBases</td><td>Integer</td><td>Number of aligned bases matching reference base for R1</td></tr>
<tr><td>R1MatchedBasesPct</td><td>Float</td><td>Aligned bases matching reference base as percentage of AlignedBases for R1</td></tr>
<tr><td>R1MismatchedBases</td><td>Integer</td><td>Number of aligned bases not matching reference base for R1</td></tr>
<tr><td>R1MismatchedBasesPct</td><td>Float</td><td>Aligned bases not matching reference base as percentage of AlignedBases for R1</td></tr>
<tr><td>R1InsertedBases</td><td>Integer</td><td>Number of insertion error bases for R1</td></tr>
<tr><td>R1InsertedBasesPct</td><td>Float</td><td>Insertion error bases as percentage of AlignedBases for R1</td></tr>
<tr><td>R1DeletedBases</td><td>Integer</td><td>Number of deletion error bases for R1</td></tr>
<tr><td>R1DeletedBasesPct</td><td>Float</td><td>Deletion error bases as percentage of AlignedBases for R1</td></tr>
<tr><td>R1SoftClippedReads</td><td>Integer</td><td>Number of reads with soft-clipping for R1</td></tr>
<tr><td>R1SoftClippedReadsPct</td><td>Float</td><td>Reads with soft-clipping as percentage of MappedReads for R1</td></tr>
<tr><td>R1SoftClippedBases</td><td>Integer</td><td>Number of soft-clipped error bases for R1</td></tr>
<tr><td>R1SoftClippedBasesPct</td><td>Float</td><td>Soft-clipped error bases as percentage of MappedBases for R1</td></tr>
<tr><td>R1PerfectReads</td><td>Integer</td><td>Number of aligned reads without any mismatches, insertions, or deletions for R1</td></tr>
<tr><td>R1PerfectReadsPct</td><td>Float</td><td>Aligned reads without any mismatches, insertions, or deletions as percentage of MappedReads for R1</td></tr>
<tr><td>R1PerfectBases</td><td>Integer</td><td>Sum of bases in PerfectBases for R1</td></tr>
<tr><td>R1PerfectBasesPct</td><td>Float</td><td>Bases in PerfectBases as percentage of MappedBases for R1</td></tr>
<tr><td>R1Q20Bases</td><td>Integer</td><td>Number of bases with quality scores of at least 20 for R1</td></tr>
<tr><td>R1Q20BasesPct</td><td>Float</td><td>Bases with quality scores of at least 20 as percentage of AlignedBases for R1</td></tr>
<tr><td>R1Q30Bases</td><td>Integer</td><td>Number of bases with quality scores of at least 30 for R1</td></tr>
<tr><td>R1Q30BasesPct</td><td>Float</td><td>Bases with quality scores of at least 30 as percentage of AlignedBases for R1</td></tr>
<tr><td>R2YieldReads</td><td>Integer</td><td>Total number of reads for R2</td></tr>
<tr><td>R2YieldBases</td><td>Integer</td><td>Total number of bases for R2</td></tr>
<tr><td>R2UnmappedReads</td><td>Integer</td><td>Number of unmapped reads for R2</td></tr>
<tr><td>R2UnmappedReadsPct</td><td>Float</td><td>Unmapped reads as percentage of YieldReads for R2</td></tr>
<tr><td>R2UnmappedBases</td><td>Integer</td><td>Sum of bases in unmapped reads for R2</td></tr>
<tr><td>R2UnmappedBasesPct</td><td>Float</td><td>Bases in unmapped reads as percentage of YieldBases for R2</td></tr>
<tr><td>R2MappedReads</td><td>Integer</td><td>Number of mapped reads for R2</td></tr>
<tr><td>R2MappedReadsPct</td><td>Float</td><td>Mapped reads as percentage of YieldReads for R2</td></tr>
<tr><td>R2MappedBases</td><td>Integer</td><td>Sum of bases in mapped reads for R2</td></tr>
<tr><td>R2MappedBasesPct</td><td>Float</td><td>Bases in mapped reads as percentage of YieldReads for R2</td></tr>
<tr><td>R2AlignedBases</td><td>Integer</td><td>MappedBases minus Soft-Clipped Bases for R2</td></tr>
<tr><td>R2AlignedBasesPct</td><td>Float</td><td>AlignedBases as percentage of YieldReads for R2</td></tr>
<tr><td>R2MatchedBases</td><td>Integer</td><td>Number of aligned bases matching reference base for R2</td></tr>
<tr><td>R2MatchedBasesPct</td><td>Float</td><td>Aligned bases matching reference base as percentage of AlignedBases for R2</td></tr>
<tr><td>R2MismatchedBases</td><td>Integer</td><td>Number of aligned bases not matching reference base for R2</td></tr>
<tr><td>R2MismatchedBasesPct</td><td>Float</td><td>Aligned bases not matching reference base as percentage of AlignedBases for R2</td></tr>
<tr><td>R2InsertedBases</td><td>Integer</td><td>Number of insertion error bases for R2</td></tr>
<tr><td>R2InsertedBasesPct</td><td>Float</td><td>Insertion error bases as percentage of AlignedBases for R2</td></tr>
<tr><td>R2DeletedBases</td><td>Integer</td><td>Number of deletion error bases for R2</td></tr>
<tr><td>R2DeletedBasesPct</td><td>Float</td><td>Deletion error bases as percentage of AlignedBases for R2</td></tr>
<tr><td>R2SoftClippedReads</td><td>Integer</td><td>Number of reads with soft-clipping for R2</td></tr>
<tr><td>R2SoftClippedReadsPct</td><td>Float</td><td>Reads with soft-clipping as percentage of MappedReads for R2</td></tr>
<tr><td>R2SoftClippedBases</td><td>Integer</td><td>Number of soft-clipped error bases for R2</td></tr>
<tr><td>R2SoftClippedBasesPct</td><td>Float</td><td>Soft-clipped error bases as percentage of MappedBases for R2</td></tr>
<tr><td>R2PerfectReads</td><td>Integer</td><td>Number of aligned reads without any mismatches, insertions, or deletions for R2</td></tr>
<tr><td>R2PerfectReadsPct</td><td>Float</td><td>Aligned reads without any mismatches, insertions, or deletions as percentage of MappedReads for R2</td></tr>
<tr><td>R2PerfectBases</td><td>Integer</td><td>Sum of bases in PerfectBases for R2</td></tr>
<tr><td>R2PerfectBasesPct</td><td>Float</td><td>Bases in PerfectBases as percentage of MappedBases for R2</td></tr>
<tr><td>R2Q20Bases</td><td>Integer</td><td>Number of bases with quality scores of at least 20 for R2</td></tr>
<tr><td>R2Q20BasesPct</td><td>Float</td><td>Bases with quality scores of at least 20 as percentage of AlignedBases for R2</td></tr>
<tr><td>R2Q30Bases</td><td>Integer</td><td>Number of bases with quality scores of at least 30 for R2</td></tr>
<tr><td>R2Q30BasesPct</td><td>Float</td><td>Bases with quality scores of at least 30 as percentage of AlignedBases for R2</td></tr>
<tr><td>AlignedReadLengthMean</td><td>Float</td><td>Mean alignment length</td></tr>
<tr><td>AlignedReadLengthMedian</td><td>Integer</td><td>Median alignment length</td></tr>
<tr><td>AlignedReadLengthMode</td><td>Integer</td><td>Mode alignment length</td></tr>
<tr><td>AlignedReadLengthStandardDeviation</td><td>Float</td><td>Standard deviation of alignment length</td></tr>
<tr><td>R1AlignedReadLengthMean</td><td>Float</td><td>Mean alignment length for R1</td></tr>
<tr><td>R1AlignedReadLengthMedian</td><td>Integer</td><td>Median alignment length for R1</td></tr>
<tr><td>R1AlignedReadLengthMode</td><td>Integer</td><td>Mode alignment length for R1</td></tr>
<tr><td>R1AlignedReadLengthStandardDeviation</td><td>Float</td><td>Standard deviation of alignment length for R1</td></tr>
<tr><td>R2AlignedReadLengthMean</td><td>Float</td><td>Mean alignment length for R2</td></tr>
<tr><td>R2AlignedReadLengthMedian</td><td>Integer</td><td>Median alignment length for R2</td></tr>
<tr><td>R2AlignedReadLengthMode</td><td>Integer</td><td>Mode alignment length for R2</td></tr>
<tr><td>R2AlignedReadLengthStandardDeviation</td><td>Float</td><td>Standard deviation of alignment length for R2</td></tr>
<tr><td>TotalPairs</td><td>Integer</td><td>Number of reads with both pairs mapped</td></tr>
<tr><td>TotalSameChrPairs</td><td>Integer</td><td>Number of reads with both pairs mapped to the same chromosome</td></tr>
<tr><td>TotalSameChrPairsPct</td><td>Float</td><td>Reads with both pairs mapped to the same chromosome as percentage of TotalPairs</td></tr>
<tr><td>UnpairedReads</td><td>Integer</td><td>Number of unpaired reads</td></tr>
<tr><td>UnpairedReadsPct</td><td>Float</td><td>Number of unpaired reads as percentage of TotalPairs</td></tr>
<tr><td>R1UnpairedReads</td><td>Integer</td><td>Number of unpaired reads for R1</td></tr>
<tr><td>R1UnpairedReadsPct</td><td>Float</td><td>Number of unpaired reads as percentage of TotalPairs for R1</td></tr>
<tr><td>R2UnpairedReads</td><td>Integer</td><td>Number of unpaired reads for R2</td></tr>
<tr><td>R2UnpairedReadsPct</td><td>Float</td><td>Number of unpaired reads as percentage of TotalPairs for R2</td></tr>
<tr><td>ChimericReadPairPct</td><td>Float</td><td>Reads in improper pair as percentage of MappedReads (samtools view -F 1806 reads / samtools view -F 1804)</td></tr>
<tr><td>InsertSizeMean</td><td>Float</td><td>Mean observed insert size (read must be in proper pair)</td></tr>
<tr><td>InsertSizeMedian</td><td>Integer</td><td>Median observed insert size (read must be in proper pair)</td></tr>
<tr><td>InsertSizeMode</td><td>Integer</td><td>Mode observed insert size (read must be in proper pair)</td></tr>
<tr><td>InsertSizeStandardDeviation</td><td>Float</td><td>Standard deviation of observed insert size (read must be in proper pair)</td></tr>
<tr><td>WgsTotalReads</td><td>Integer</td><td>Total number of reads</td></tr>
<tr><td>WgsCovDuplicateReads</td><td>Integer</td><td>Number of duplicate reads</td></tr>
<tr><td>WgsCovDuplicateReadsPct</td><td>Float</td><td>Duplicate reads as percentage of total reads</td></tr>
<tr><td>WgsAlignedReads</td><td>Integer</td><td>Number of aligned reads</td></tr>
<tr><td>WgsAlignedReadsPct</td><td>Float</td><td>Aligned reads as percentage of total reads</td></tr>
<tr><td>WgsReadsPaired</td><td>Integer</td><td>Number of paired reads</td></tr>
<tr><td>WgsReadsPairedWithMates</td><td>Integer</td><td>Number of paired reads with mate mapped</td></tr>
<tr><td>WgsCoverageMean</td><td>Float</td><td>Mean coverage</td></tr>
<tr><td>WgsCoverageMedian</td><td>Integer</td><td>Median coverage</td></tr>
<tr><td>WgsCoverageStandardDeviation</td><td>Float</td><td>Standard deviation of coverage</td></tr>
<tr><td>WgsExpectedAlignedReads</td><td>Integer</td><td>Number of aligned reads</td></tr>
<tr><td>WgsCalculatedAlignedReads</td><td>Integer</td><td>Sum of reads in target, reads in buffer, and reads out of target and buffer</td></tr>
<tr><td>WgsCoverageBases1</td><td>Integer</td><td>Number of bases with coverage of at least 1 read</td></tr>
<tr><td>WgsCoverageBases1Pct</td><td>Float</td><td>Bases with coverage of at least 1 read as percentage of total bases</td></tr>
<tr><td>WgsCoverageBases10</td><td>Integer</td><td>Number of bases with coverage of at least 10 reads</td></tr>
<tr><td>WgsCoverageBases10Pct</td><td>Float</td><td>Bases with coverage of at least 10 reads as percentage of total bases</td></tr>
<tr><td>WgsCoverageBases20</td><td>Integer</td><td>Number of bases with coverage of at least 20 reads</td></tr>
<tr><td>WgsCoverageBases20Pct</td><td>Float</td><td>Bases with coverage of at least 20 reads as percentage of total bases</td></tr>
<tr><td>WgsCoverageBases30</td><td>Integer</td><td>Number of bases with coverage of at least 30 reads</td></tr>
<tr><td>WgsCoverageBases30Pct</td><td>Float</td><td>Bases with coverage of at least 30 reads as percentage of total bases</td></tr>
<tr><td>WgsCoverageBases40</td><td>Integer</td><td>Number of bases with coverage of at least 40 reads</td></tr>
<tr><td>WgsCoverageBases40Pct</td><td>Float</td><td>Bases with coverage of at least 40 reads as percentage of total bases</td></tr>
<tr><td>WgsCoverageBases50</td><td>Integer</td><td>Number of bases with coverage of at least 50 reads</td></tr>
<tr><td>WgsCoverageBases50Pct</td><td>Float</td><td>Bases with coverage of at least 50 reads as percentage of total bases</td></tr>
<tr><td>WgsCoverageBases60</td><td>Integer</td><td>Number of bases with coverage of at least 60 reads</td></tr>
<tr><td>WgsCoverageBases60Pct</td><td>Float</td><td>Bases with coverage of at least 60 reads as percentage of total bases</td></tr>
<tr><td>WgsCoverageBases70</td><td>Integer</td><td>Number of bases with coverage of at least 70 reads</td></tr>
<tr><td>WgsCoverageBases70Pct</td><td>Float</td><td>Bases with coverage of at least 70 reads as percentage of total bases</td></tr>
<tr><td>WgsCoverageBases100</td><td>Integer</td><td>Number of bases with coverage of at least 100 reads</td></tr>
<tr><td>WgsCoverageBases100Pct</td><td>Float</td><td>Bases with coverage of at least 100 reads as percentage of total bases</td></tr>
<tr><td>WgsCoverageBases500</td><td>Integer</td><td>Number of bases with coverage of at least 500 reads</td></tr>
<tr><td>WgsCoverageBases500Pct</td><td>Float</td><td>Bases with coverage of at least 100 reads as percentage of total bases</td></tr>
<tr><td>WgsCoverageBases1000</td><td>Integer</td><td>Number of bases with coverage of at least 1000 reads</td></tr>
<tr><td>WgsCoverageBases1000Pct</td><td>Float</td><td>Bases with coverage of at least 1000 reads as percentage of total bases</td></tr>
<tr><td>CapTotalReads</td><td>Integer</td><td>Total number of reads</td></tr>
<tr><td>CapCovDuplicateReads</td><td>Integer</td><td>Number of duplicate reads</td></tr>
<tr><td>CapCovDuplicateReadsPct</td><td>Float</td><td>Duplicate reads as percentage of total reads</td></tr>
<tr><td>CapAlignedReads</td><td>Integer</td><td>Number of aligned reads</td></tr>
<tr><td>CapAlignedReadsPct</td><td>Float</td><td>Aligned reads as percentage of total reads</td></tr>
<tr><td>CapReadsPaired</td><td>Integer</td><td>Number of paired reads</td></tr>
<tr><td>CapReadsPairedWithMates</td><td>Integer</td><td>Number of paired reads with mate mapped</td></tr>
<tr><td>CapCoverageMean</td><td>Float</td><td>Mean coverage within target regions</td></tr>
<tr><td>CapCoverageMedian</td><td>Integer</td><td>Median coverage within target regions</td></tr>
<tr><td>CapCoverageStandardDeviation</td><td>Float</td><td>Standard deviation of coverage within target regions</td></tr>
<tr><td>CapExpectedAlignedReads</td><td>Integer</td><td>Number of aligned reads</td></tr>
<tr><td>CapCalculatedAlignedReads</td><td>Integer</td><td>Sum of reads in target, reads in buffer, and reads out of target and buffer</td></tr>
<tr><td>CapCoverageBases1</td><td>Integer</td><td>Bases with coverage of at least 1 read within target regions</td></tr>
<tr><td>CapCoverageBases1Pct</td><td>Float</td><td>Bases with coverage of at least 1 read as percentage of total bases within target regions</td></tr>
<tr><td>CapCoverageBases10</td><td>Integer</td><td>Bases with coverage of at least 10 reads within target regions</td></tr>
<tr><td>CapCoverageBases10Pct</td><td>Float</td><td>Bases with coverage of at least 10 reads as percentage of total bases within target regions</td></tr>
<tr><td>CapCoverageBases20</td><td>Integer</td><td>Bases with coverage of at least 20 reads within target regions</td></tr>
<tr><td>CapCoverageBases20Pct</td><td>Float</td><td>Bases with coverage of at least 20 reads as percentage of total bases within target regions</td></tr>
<tr><td>CapCoverageBases30</td><td>Integer</td><td>Bases with coverage of at least 30 reads within target regions</td></tr>
<tr><td>CapCoverageBases30Pct</td><td>Float</td><td>Bases with coverage of at least 30 reads as percentage of total bases within target regions</td></tr>
<tr><td>CapCoverageBases40</td><td>Integer</td><td>Bases with coverage of at least 40 reads within target regions</td></tr>
<tr><td>CapCoverageBases40Pct</td><td>Float</td><td>Bases with coverage of at least 40 reads as percentage of total bases within target regions</td></tr>
<tr><td>CapCoverageBases50</td><td>Integer</td><td>Bases with coverage of at least 50 reads within target regions</td></tr>
<tr><td>CapCoverageBases50Pct</td><td>Float</td><td>Bases with coverage of at least 50 reads as percentage of total bases within target regions</td></tr>
<tr><td>CapCoverageBases60</td><td>Integer</td><td>Bases with coverage of at least 60 reads within target regions</td></tr>
<tr><td>CapCoverageBases60Pct</td><td>Float</td><td>Bases with coverage of at least 60 reads as percentage of total bases within target regions</td></tr>
<tr><td>CapCoverageBases70</td><td>Integer</td><td>Bases with coverage of at least 70 reads within target regions</td></tr>
<tr><td>CapCoverageBases70Pct</td><td>Float</td><td>Bases with coverage of at least 70 reads as percentage of total bases within target regions</td></tr>
<tr><td>CapCoverageBases100</td><td>Integer</td><td>Bases with coverage of at least 100 reads within target regions</td></tr>
<tr><td>CapCoverageBases100Pct</td><td>Float</td><td>Bases with coverage of at least 100 reads as percentage of total bases within target regions</td></tr>
<tr><td>CapCoverageBases500</td><td>Integer</td><td>Bases with coverage of at least 500 reads within target regions</td></tr>
<tr><td>CapCoverageBases500Pct</td><td>Float</td><td>Bases with coverage of at least 500 reads as percentage of total bases within target regions</td></tr>
<tr><td>CapCoverageBases1000</td><td>Integer</td><td>Bases with coverage of at least 1000 reads within target regions</td></tr>
<tr><td>CapCoverageBases1000Pct</td><td>Float</td><td>Bases with coverage of at least 1000 reads as percentage of total bases within target regions</td></tr>
<tr><td>CapBufferAlignedReads</td><td>Integer</td><td>Number of aligned reads in buffer regions but not in target regions</td></tr>
<tr><td>CapBufferAlignedReadsPct</td><td>Float</td><td>Number of aligned reads in buffer regions but not in target regions as percentage of aligned reads</td></tr>
<tr><td>CapTargetAlignedReads</td><td>Integer</td><td>Number of aligned reads in target regions</td></tr>
<tr><td>CapTargetAlignedReadsPct</td><td>Float</td><td>Aligned reads in target regions as percentage of aligned reads</td></tr>
<tr><td>CapTargetsHit</td><td>Integer</td><td>Number of targets capturing at least one aligned read</td></tr>
<tr><td>CapTargetsHitPct</td><td>Float</td><td>Targets capturing at least one aligned read as percentage of total targets</td></tr>
<tr><td>CapTargetBuffersHit</td><td>Integer</td><td>Number of targets not capturing any reads whose buffer regions capture at least one aligned read</td></tr>
<tr><td>CapTargetBuffersHitPct</td><td>Float</td><td>Targets not capturing any reads whose buffer regions capture at least one aligned read as percentage of total targets</td></tr>
<tr><td>CapTotalTargets</td><td>Integer</td><td>Total number of targets from capture file</td></tr>
<tr><td>CapHighCoverageNonTargetHits</td><td>Integer</td><td>Number of contiguous regions of coverage outside of target regions or target buffer regions with coverage of at least 20 for at least one base</td></tr>
<tr><td>CapBasesOnTarget</td><td>Integer</td><td>Number of bases within target regions</td></tr>
<tr><td>CapBasesOnBuffer</td><td>Integer</td><td>Number of bases within target buffer regions but outside of target regions</td></tr>
<tr><td>CapReadsOnTargetOrBuffer</td><td>Integer</td><td>Number of reads captured by either target regions or target buffer regions</td></tr>
<tr><td>CapReadsOnTargetOrBufferPct</td><td>Float</td><td>Reads captured by either target regions or target buffer regions as percentage of aligned reads</td></tr>
<tr><td>FilteredOverlapBases</td><td>Integer</td><td>Total per-read bases filtered from coverage statistics as overlapping with bases in the mate pair read</td></tr>
<tr><td>FilteredLowBaseQualityBases</td><td>Integer</td><td>Total per-read bases filtered from coverage statistics as having a base quality below the minimum base quality</td></tr>
</table>
