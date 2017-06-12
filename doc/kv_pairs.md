# AlignStats Key-Value Pairs

<style>
table, th, td {
    border-color:black;
    border-style:solid;
    border-width:1px;
}
</style>

<table>
<tr><th>Key</th><th>Type</th><th>Section</th><th>Description</th></tr>
<tr><td>Yield_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Total number of reads</td></tr>
<tr><td>Yield_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Total number of bases</td></tr>
<tr><td>Unmapped_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of unmapped reads</td></tr>
<tr><td>Unmapped_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Unmapped reads as percentage of Yield_Reads</td></tr>
<tr><td>Unmapped_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Sum of bases in unmapped reads</td></tr>
<tr><td>Unmapped_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Bases in unmapped reads as percentage of Yield_Bases</td></tr>
<tr><td>Duplicate_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of mapped duplicate reads</td></tr>
<tr><td>Duplicate_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Mapped duplicate reads as percentage of aligned reads</td></tr>
<tr><td>Duplicate_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Sum of bases in duplicate reads</td></tr>
<tr><td>Duplicate_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Sum of bases in duplicate reads as percentage of Aligned_Bases</td></tr>
<tr><td>Mapped_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of mapped reads</td></tr>
<tr><td>Mapped_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Mapped reads as percentage of Yield_Reads</td></tr>
<tr><td>Mapped_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Sum of bases in mapped reads</td></tr>
<tr><td>Mapped_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Bases in mapped reads as percentage of Yield_Reads</td></tr>
<tr><td>Aligned_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Mapped_Bases minus Soft-Clipped Bases</td></tr>
<tr><td>Aligned_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Aligned_Bases as percentage of Yield_Reads</td></tr>
<tr><td>Matched_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of aligned bases matching reference base</td></tr>
<tr><td>Matched_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Aligned bases matching reference base as percentage of Aligned_Bases</td></tr>
<tr><td>Mismatched_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of aligned bases not matching reference base</td></tr>
<tr><td>Mismatched_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Aligned bases not matching reference base as percentage of Aligned_Bases</td></tr>
<tr><td>Inserted_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of insertion error bases</td></tr>
<tr><td>Inserted_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Insertion error bases as percentage of Aligned_Bases</td></tr>
<tr><td>Deleted_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of deletion error bases</td></tr>
<tr><td>Deleted_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Deletion error bases as percentage of Aligned_Bases</td></tr>
<tr><td>SoftClipped_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of reads with soft-clipping</td></tr>
<tr><td>SoftClipped_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Reads with soft-clipping as percentage of Mapped_Reads</td></tr>
<tr><td>SoftClipped_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of soft-clipped error bases</td></tr>
<tr><td>SoftClipped_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Soft-clipped error bases as percentage of Mapped_Bases</td></tr>
<tr><td>Perfect_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of aligned reads without any mismatches, insertions, or deletions</td></tr>
<tr><td>Perfect_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Aligned reads without any mismatches, insertions, or deletions as percentage of Mapped_Reads</td></tr>
<tr><td>Perfect_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Sum of bases in Perfect_Bases</td></tr>
<tr><td>Perfect_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Bases in Perfect_Bases as percentage of Mapped_Bases</td></tr>
<tr><td>Q20_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of bases with quality scores of at least 20</td></tr>
<tr><td>Q20_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Bases with quality scores of at least 20 as percentage of Aligned_Bases</td></tr>
<tr><td>R1_Yield_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Total number of reads for R1</td></tr>
<tr><td>R1_Yield_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Total number of bases for R1</td></tr>
<tr><td>R1_Unmapped_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of unmapped reads for R1</td></tr>
<tr><td>R1_Unmapped_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Unmapped reads as percentage of Yield_Reads for R1</td></tr>
<tr><td>R1_Unmapped_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Sum of bases in unmapped reads for R1</td></tr>
<tr><td>R1_Unmapped_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Bases in unmapped reads as percentage of Yield_Bases for R1</td></tr>
<tr><td>R1_Mapped_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of mapped reads for R1</td></tr>
<tr><td>R1_Mapped_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Mapped reads as percentage of Yield_Reads for R1</td></tr>
<tr><td>R1_Mapped_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Sum of bases in mapped reads for R1</td></tr>
<tr><td>R1_Mapped_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Sum of bases in mapped reads as percentage of Yield_Reads for R1</td></tr>
<tr><td>R1_Aligned_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Mapped_Bases minus Soft-Clipped Bases for R1</td></tr>
<tr><td>R1_Aligned_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Aligned_Bases as percentage of Yield_Reads for R1</td></tr>
<tr><td>R1_Matched_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of aligned bases matching reference base for R1</td></tr>
<tr><td>R1_Matched_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Aligned bases matching reference base as percentage of Aligned_Bases for R1</td></tr>
<tr><td>R1_Mismatched_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of aligned bases not matching reference base for R1</td></tr>
<tr><td>R1_Mismatched_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Aligned bases not matching reference base as percentage of Aligned_Bases for R1</td></tr>
<tr><td>R1_Inserted_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of insertion error bases for R1</td></tr>
<tr><td>R1_Inserted_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Insertion error bases as percentage of Aligned_Bases for R1</td></tr>
<tr><td>R1_Deleted_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of deletion error bases for R1</td></tr>
<tr><td>R1_Deleted_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Deletion error bases as percentage of Aligned_Bases for R1</td></tr>
<tr><td>R1_SoftClipped_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of reads with soft-clipping for R1</td></tr>
<tr><td>R1_SoftClipped_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Reads with soft-clipping as percentage of Mapped_Reads for R1</td></tr>
<tr><td>R1_SoftClipped_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of soft-clipped error bases for R1</td></tr>
<tr><td>R1_SoftClipped_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Soft-clipped error bases as percentage of Mapped_Bases for R1</td></tr>
<tr><td>R1_Perfect_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of aligned reads without any mismatches, insertions, or deletions for R1</td></tr>
<tr><td>R1_Perfect_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Aligned reads without any mismatches, insertions, or deletions as percentage of Mapped_Reads for R1</td></tr>
<tr><td>R1_Perfect_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Sum of bases in Perfect_Bases for R1</td></tr>
<tr><td>R1_Perfect_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Bases in Perfect_Bases as percentage of Mapped_Bases for R1</td></tr>
<tr><td>R1_Q20_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of bases with quality scores of at least 20 for R1</td></tr>
<tr><td>R1_Q20_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Bases with quality scores of at least 20 as percentage of Aligned_Bases for R1</td></tr>
<tr><td>R2_Yield_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Total number of reads for R2</td></tr>
<tr><td>R2_Yield_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Total number of bases for R2</td></tr>
<tr><td>R2_Unmapped_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of unmapped reads for R2</td></tr>
<tr><td>R2_Unmapped_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Unmapped reads as percentage of Yield_Reads for R2</td></tr>
<tr><td>R2_Unmapped_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Sum of bases in unmapped reads for R2</td></tr>
<tr><td>R2_Unmapped_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Bases in unmapped reads as percentage of Yield_Bases for R2</td></tr>
<tr><td>R2_Mapped_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of mapped reads for R2</td></tr>
<tr><td>R2_Mapped_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Mapped reads as percentage of Yield_Reads for R2</td></tr>
<tr><td>R2_Mapped_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Sum of bases in mapped reads for R2</td></tr>
<tr><td>R2_Mapped_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Bases in mapped reads as percentage of Yield_Reads for R2</td></tr>
<tr><td>R2_Aligned_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Mapped_Bases minus Soft-Clipped Bases for R2</td></tr>
<tr><td>R2_Aligned_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Aligned_Bases as percentage of Yield_Reads for R2</td></tr>
<tr><td>R2_Matched_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of aligned bases matching reference base for R2</td></tr>
<tr><td>R2_Matched_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Aligned bases matching reference base as percentage of Aligned_Bases for R2</td></tr>
<tr><td>R2_Mismatched_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of aligned bases not matching reference base for R2</td></tr>
<tr><td>R2_Mismatched_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Aligned bases not matching reference base as percentage of Aligned_Bases for R2</td></tr>
<tr><td>R2_Inserted_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of insertion error bases for R2</td></tr>
<tr><td>R2_Inserted_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Insertion error bases as percentage of Aligned_Bases for R2</td></tr>
<tr><td>R2_Deleted_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of deletion error bases for R2</td></tr>
<tr><td>R2_Deleted_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Deletion error bases as percentage of Aligned_Bases for R2</td></tr>
<tr><td>R2_SoftClipped_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of reads with soft-clipping for R2</td></tr>
<tr><td>R2_SoftClipped_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Reads with soft-clipping as percentage of Mapped_Reads for R2</td></tr>
<tr><td>R2_SoftClipped_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of soft-clipped error bases for R2</td></tr>
<tr><td>R2_SoftClipped_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Soft-clipped error bases as percentage of Mapped_Bases for R2</td></tr>
<tr><td>R2_Perfect_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of aligned reads without any mismatches, insertions, or deletions for R2</td></tr>
<tr><td>R2_Perfect_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Aligned reads without any mismatches, insertions, or deletions as percentage of Mapped_Reads for R2</td></tr>
<tr><td>R2_Perfect_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Sum of bases in Perfect_Bases for R2</td></tr>
<tr><td>R2_Perfect_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Bases in Perfect_Bases as percentage of Mapped_Bases for R2</td></tr>
<tr><td>R2_Q20_Bases</td><td>Integer</td><td>Aligment Statistics</td><td>Number of bases with quality scores of at least 20 for R2</td></tr>
<tr><td>R2_Q20_Bases_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Bases with quality scores of at least 20 as percentage of Aligned_Bases for R2</td></tr>
<tr><td>Mean_Aligned_Read_Length</td><td>Float</td><td>Aligment Statistics</td><td>Mean alignment length</td></tr>
<tr><td>Median_Aligned_Read_Length</td><td>Integer</td><td>Aligment Statistics</td><td>Median alignment length</td></tr>
<tr><td>Mode_Aligned_Read_Length</td><td>Integer</td><td>Aligment Statistics</td><td>Mode alignment length</td></tr>
<tr><td>R1_Mean_Aligned_Read_Length</td><td>Float</td><td>Aligment Statistics</td><td>Mean alignment length for R1</td></tr>
<tr><td>R1_Median_Aligned_Read_Length</td><td>Integer</td><td>Aligment Statistics</td><td>Median alignment length for R1</td></tr>
<tr><td>R1_Mode_Aligned_Read_Length</td><td>Integer</td><td>Aligment Statistics</td><td>Mode alignment length for R1</td></tr>
<tr><td>R2_Mean_Aligned_Read_Length</td><td>Float</td><td>Aligment Statistics</td><td>Mean alignment length for R2</td></tr>
<tr><td>R2_Median_Aligned_Read_Length</td><td>Integer</td><td>Aligment Statistics</td><td>Median alignment length for R2</td></tr>
<tr><td>R2_Mode_Aligned_Read_Length</td><td>Integer</td><td>Aligment Statistics</td><td>Mode alignment length for R2</td></tr>
<tr><td>Total_Pairs</td><td>Integer</td><td>Aligment Statistics</td><td>Number of reads with both pairs mapped</td></tr>
<tr><td>Total_Same_Chr_Pairs</td><td>Integer</td><td>Aligment Statistics</td><td>Number of reads with both pairs mapped to the same chromosome</td></tr>
<tr><td>Total_Same_Chr_Pairs_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Reads with both pairs mapped to the same chromosome as percentage of Total_Pairs</td></tr>
<tr><td>Unpaired_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of unpaired reads</td></tr>
<tr><td>Unpaired_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Number of unpaired reads as percentage of Total_Pairs</td></tr>
<tr><td>R1_Unpaired_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of unpaired reads for R1</td></tr>
<tr><td>R1_Unpaired_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Number of unpaired reads as percentage of Total_Pairs for R1</td></tr>
<tr><td>R2_Unpaired_Reads</td><td>Integer</td><td>Aligment Statistics</td><td>Number of unpaired reads for R2</td></tr>
<tr><td>R2_Unpaired_Reads_Pct</td><td>Float</td><td>Aligment Statistics</td><td>Number of unpaired reads as percentage of Total_Pairs for R2</td></tr>
<tr><td>Chimeric_Rate</td><td>Float</td><td>Aligment Statistics</td><td>Reads in improper pair as percentage of Mapped_Reads (samtools view -F 1806 reads / samtools view -F 1804)</td></tr>
<tr><td>Mean_Insert_Size</td><td>Float</td><td>Aligment Statistics</td><td>Mean observed insert size (read must be in proper pair)</td></tr>
<tr><td>Median_Insert_Size</td><td>Integer</td><td>Aligment Statistics</td><td>Median observed insert size (read must be in proper pair)</td></tr>
<tr><td>Mode_Insert_Size</td><td>Integer</td><td>Aligment Statistics</td><td>Mode observed insert size (read must be in proper pair)</td></tr>
<tr><td>Wgs_Total_Reads</td><td>Integer</td><td>Whole Genome Metrics</td><td>Total number of reads</td></tr>
<tr><td>Wgs_Coverage_Duplicate_Reads</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of duplicate reads</td></tr>
<tr><td>Wgs_Coverage_Duplicate_Reads_Pct</td><td>Float</td><td>Whole Genome Metrics</td><td>Duplicate reads as percentage of total reads</td></tr>
<tr><td>Wgs_Aligned_Reads</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of aligned reads</td></tr>
<tr><td>Wgs_Aligned_Reads_Pct</td><td>Float</td><td>Whole Genome Metrics</td><td>Aligned reads as percentage of total reads</td></tr>
<tr><td>Wgs_Reads_Paired</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of paired reads</td></tr>
<tr><td>Wgs_Reads_Paired_With_Mates</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of paired reads with mate mapped</td></tr>
<tr><td>Wgs_Average_Coverage</td><td>Float</td><td>Whole Genome Metrics</td><td>Average coverge</td></tr>
<tr><td>Wgs_Median_Coverage</td><td>Integer</td><td>Whole Genome Metrics</td><td>Median coverage</td></tr>
<tr><td>Wgs_Expected_Aligned_Reads</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of aligned reads</td></tr>
<tr><td>Wgs_Calculated_Aligned_Reads</td><td>Integer</td><td>Whole Genome Metrics</td><td>Sum of reads in target, reads in buffer, and reads out of target and buffer</td></tr>
<tr><td>Wgs_Coverage_Bases_1</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of bases with coverage of at least 1 read</td></tr>
<tr><td>Wgs_Coverage_Bases_1_Pct</td><td>Float</td><td>Whole Genome Metrics</td><td>Bases with coverage of at least 1 read as percentage of total bases</td></tr>
<tr><td>Wgs_Coverage_Bases_10</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of bases with coverage of at least 10 reads</td></tr>
<tr><td>Wgs_Coverage_Bases_10_Pct</td><td>Float</td><td>Whole Genome Metrics</td><td>Bases with coverage of at least 10 reads as percentage of total bases</td></tr>
<tr><td>Wgs_Coverage_Bases_20</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of bases with coverage of at least 20 reads</td></tr>
<tr><td>Wgs_Coverage_Bases_20_Pct</td><td>Float</td><td>Whole Genome Metrics</td><td>Bases with coverage of at least 20 reads as percentage of total bases</td></tr>
<tr><td>Wgs_Coverage_Bases_30</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of bases with coverage of at least 30 reads</td></tr>
<tr><td>Wgs_Coverage_Bases_30_Pct</td><td>Float</td><td>Whole Genome Metrics</td><td>Bases with coverage of at least 30 reads as percentage of total bases</td></tr>
<tr><td>Wgs_Coverage_Bases_40</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of bases with coverage of at least 40 reads</td></tr>
<tr><td>Wgs_Coverage_Bases_40_Pct</td><td>Float</td><td>Whole Genome Metrics</td><td>Bases with coverage of at least 40 reads as percentage of total bases</td></tr>
<tr><td>Wgs_Coverage_Bases_50</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of bases with coverage of at least 50 reads</td></tr>
<tr><td>Wgs_Coverage_Bases_50_Pct</td><td>Float</td><td>Whole Genome Metrics</td><td>Bases with coverage of at least 50 reads as percentage of total bases</td></tr>
<tr><td>Wgs_Coverage_Bases_100</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of bases with coverage of at least 100 reads</td></tr>
<tr><td>Wgs_Coverage_Bases_100_Pct</td><td>Float</td><td>Whole Genome Metrics</td><td>Bases with coverage of at least 100 reads as percentage of total bases</td></tr>
<tr><td>Wgs_Coverage_Bases_1000</td><td>Integer</td><td>Whole Genome Metrics</td><td>Number of bases with coverage of at least 1000 reads</td></tr>
<tr><td>Wgs_Coverage_Bases_1000_Pct</td><td>Float</td><td>Whole Genome Metrics</td><td>Bases with coverage of at least 1000 reads as percentage of total bases</td></tr>
<tr><td>Cap_Total_Reads</td><td>Integer</td><td>Capture Metrics</td><td>Total number of reads</td></tr>
<tr><td>Cap_Coverage_Duplicate_Reads</td><td>Integer</td><td>Capture Metrics</td><td>Number of duplicate reads</td></tr>
<tr><td>Cap_Coverage_Duplicate_Reads_Pct</td><td>Float</td><td>Capture Metrics</td><td>Duplicate reads as percentage of total reads</td></tr>
<tr><td>Cap_Aligned_Reads</td><td>Integer</td><td>Capture Metrics</td><td>Number of aligned reads</td></tr>
<tr><td>Cap_Aligned_Reads_Pct</td><td>Float</td><td>Capture Metrics</td><td>Aligned reads as percentage of total reads</td></tr>
<tr><td>Cap_Reads_Paired</td><td>Integer</td><td>Capture Metrics</td><td>Number of paired reads</td></tr>
<tr><td>Cap_Reads_Paired_With_Mates</td><td>Integer</td><td>Capture Metrics</td><td>Number of paired reads with mate mapped</td></tr>
<tr><td>Cap_Average_Coverage</td><td>Float</td><td>Capture Metrics</td><td>Average coverge within target regions</td></tr>
<tr><td>Cap_Median_Coverage</td><td>Integer</td><td>Capture Metrics</td><td>Median coverage within target regions</td></tr>
<tr><td>Cap_Expected_Aligned_Reads</td><td>Integer</td><td>Capture Metrics</td><td>Number of aligned reads</td></tr>
<tr><td>Cap_Calculated_Aligned_Reads</td><td>Integer</td><td>Capture Metrics</td><td>Sum of reads in target, reads in buffer, and reads out of target and buffer</td></tr>
<tr><td>Cap_Coverage_Bases_1</td><td>Integer</td><td>Capture Metrics</td><td>Bases with coverage of at least 1 read within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_1_Pct</td><td>Float</td><td>Capture Metrics</td><td>Bases with coverage of at least 1 read as percentage of total bases within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_10</td><td>Integer</td><td>Capture Metrics</td><td>Bases with coverage of at least 10 reads within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_10_Pct</td><td>Float</td><td>Capture Metrics</td><td>Bases with coverage of at least 10 reads as percentage of total bases within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_20</td><td>Integer</td><td>Capture Metrics</td><td>Bases with coverage of at least 20 reads within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_20_Pct</td><td>Float</td><td>Capture Metrics</td><td>Bases with coverage of at least 20 reads as percentage of total bases within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_30</td><td>Integer</td><td>Capture Metrics</td><td>Bases with coverage of at least 30 reads within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_30_Pct</td><td>Float</td><td>Capture Metrics</td><td>Bases with coverage of at least 30 reads as percentage of total bases within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_40</td><td>Integer</td><td>Capture Metrics</td><td>Bases with coverage of at least 40 reads within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_40_Pct</td><td>Float</td><td>Capture Metrics</td><td>Bases with coverage of at least 40 reads as percentage of total bases within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_50</td><td>Integer</td><td>Capture Metrics</td><td>Bases with coverage of at least 50 reads within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_50_Pct</td><td>Float</td><td>Capture Metrics</td><td>Bases with coverage of at least 50 reads as percentage of total bases within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_100</td><td>Integer</td><td>Capture Metrics</td><td>Bases with coverage of at least 100 reads within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_100_Pct</td><td>Float</td><td>Capture Metrics</td><td>Bases with coverage of at least 100 reads as percentage of total bases within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_1000</td><td>Integer</td><td>Capture Metrics</td><td>Bases with coverage of at least 1000 reads within target regions</td></tr>
<tr><td>Cap_Coverage_Bases_1000_Pct</td><td>Float</td><td>Capture Metrics</td><td>Bases with coverage of at least 1000 reads as percentage of total bases within target regions</td></tr>
<tr><td>Cap_Buffer_Aligned_Reads</td><td>Integer</td><td>Capture Metrics</td><td>Number of aligned reads in buffer regions but not in target regions</td></tr>
<tr><td>Cap_Buffer_Aligned_Reads_Pct</td><td>Float</td><td>Capture Metrics</td><td>Number of aligned reads in buffer regions but not in target regions as percentage of aligned reads</td></tr>
<tr><td>Cap_Target_Aligned_Reads</td><td>Integer</td><td>Capture Metrics</td><td>Number of aligned reads in target regions</td></tr>
<tr><td>Cap_Target_Aligned_Reads_Pct</td><td>Float</td><td>Capture Metrics</td><td>Aligned reads in target regions as percentage of aligned reads</td></tr>
<tr><td>Cap_Targets_Hit</td><td>Integer</td><td>Capture Metrics</td><td>Number of targets capturing at least one aligned read</td></tr>
<tr><td>Cap_Targets_Hit_Pct</td><td>Float</td><td>Capture Metrics</td><td>Targets capturing at least one aligned read as percentage of total targets</td></tr>
<tr><td>Cap_Target_Buffers_Hit</td><td>Integer</td><td>Capture Metrics</td><td>Number of targets not capturing any reads whose buffer regions capture at least one alinged read</td></tr>
<tr><td>Cap_Target_Buffers_Hit_Pct</td><td>Float</td><td>Capture Metrics</td><td>Targets not capturing any reads whose buffer regions capture at least one alinged read as percentage of total targets</td></tr>
<tr><td>Cap_Total_Targets</td><td>Integer</td><td>Capture Metrics</td><td>Total number of targets from capture file</td></tr>
<tr><td>Cap_High_Coverage_Non_Target_Hits</td><td>Integer</td><td>Capture Metrics</td><td>Number of contiguous regions of coverage outside of target regions or target buffer regions with coverage of at least 20 for at least one base</td></tr>
<tr><td>Cap_Bases_On_Target</td><td>Integer</td><td>Capture Metrics</td><td>Number of bases within target regions</td></tr>
<tr><td>Cap_Bases_On_Buffer</td><td>Integer</td><td>Capture Metrics</td><td>Number of bases within target buffer regions but outside of target regions</td></tr>
<tr><td>Cap_Reads_On_Target_Or_Buffer</td><td>Integer</td><td>Capture Metrics</td><td>Number of reads captured by either target regions or target buffer regions</td></tr>
<tr><td>Cap_Reads_On_Target_Or_Buffer_Pct</td><td>Float</td><td>Capture Metrics</td><td>Reads captured by either target regions or target buffer regions as percentage of aligned reads</td></tr>
</table>
