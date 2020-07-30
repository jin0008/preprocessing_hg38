


REF=/media/hanjinu/SS200/db/refs/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
INTERVAL=/media/hanjinu/SS200/db/refs/interval_list/whole.exome.hg38.gene.interval_list


dbSNP_vcf=/media/hanjinu/SS200/db/refs/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
known_indels_sites_VCF1=/media/hanjinu/SS200/db/refs/hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
known_indels_sites_VCF2=/media/hanjinu/SS200/db/refs/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz
recalibration_report_filename=recalibration_report

#define Sample name
Sample=10075733

#define thread number
thread=20

F1=${Sample}_R1.fastq.gz
F2=${Sample}_R2.fastq.gz
sequencing_center=Severance
compression_level=5
java_opt=-Xmx32g
metrix_filename=mark_duplicate_metrix


#unmapped bam
gatk --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" FastqToSam \
-F1 $F1 \
-F2 $F2 \
-O ${Sample}.unmapped.bam \
-SM ${Sample} \
-PL illumina \
-LB ${Sample} \
-RG ${Sample} \
-CN ${sequencing_center} \
--TMP_DIR $PWD \
--VERBOSITY ERROR

gatk --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" MarkIlluminaAdapters \
-I ${Sample}.unmapped.bam \
-O ${Sample}.markilluminaadapters.bam \
-M ${Sample}.markilluminaadapters_metrics.txt \
--TMP_DIR $PWD \
--VERBOSITY ERROR

gatk --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" SamToFastq \
-I ${Sample}.markilluminaadapters.bam \
-FASTQ /dev/stdout \
-CLIPPING_ATTRIBUTE XT \
-CLIPPING_ACTION 2 \
-INTERLEAVE true \
-NON_PF true \
--VERBOSITY ERROR \
--TMP_DIR $PWD | \
bwa mem -t ${thread} -M -p ${REF} /dev/stdin | \
gatk --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" MergeBamAlignment \
--VALIDATION_STRINGENCY SILENT \
--EXPECTED_ORIENTATIONS FR \
--ATTRIBUTES_TO_RETAIN X0 \
--ALIGNED_BAM /dev/stdin \
--UNMAPPED_BAM ${Sample}.unmapped.bam \
--OUTPUT ${Sample}.bam \
--REFERENCE_SEQUENCE ${REF} \
--PAIRED_RUN true \
--SORT_ORDER "unsorted" \
--IS_BISULFITE_SEQUENCE false \
--ALIGNED_READS_ONLY false \
--CLIP_ADAPTERS false \
--MAX_RECORDS_IN_RAM 2000000 \
--ADD_MATE_CIGAR true \
--MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
--UNMAPPED_READ_STRATEGY COPY_TO_TAG \
--ALIGNER_PROPER_PAIR_FLAGS true \
--UNMAP_CONTAMINANT_READS true \
--TMP_DIR $PWD \
--VERBOSITY ERROR 

#sortsam
gatk --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
SortSam \
--INPUT ${Sample}.bam \
--OUTPUT /dev/stdout \
--SORT_ORDER "coordinate" \
--CREATE_INDEX false \
--CREATE_MD5_FILE false \
--TMP_DIR $PWD \
--VERBOSITY ERROR \
| \
gatk --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
SetNmMdAndUqTags \
--INPUT /dev/stdin \
--OUTPUT ${Sample}.sorted.bam \
--CREATE_INDEX true \
--CREATE_MD5_FILE true \
--REFERENCE_SEQUENCE ${REF} \
--TMP_DIR $PWD \
--VERBOSITY ERROR \

rm -rf ${Sample}.bam

#markduplicatesSpark (#--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \ for patterned flow cell model)
gatk --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" MarkDuplicates \
-I ${Sample}.sorted.bam \
-O ${Sample}.markduplicates.bam \
-M ${Sample}.mark_duplicate_metrix.txt \
-ASO queryname \
--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
--CREATE_INDEX true \
--TMP_DIR $PWD \
--VERBOSITY ERROR 

rm -rf ${Sample}.sorted.bam
rm -rf ${Sample}.sorted.bai
rm -rf ${Sample}.sorted.bam.md5

samtools index -@ {thread} ${Sample}.markduplicates.bam

#BaseRecalibrator 
gatk --java-options "${java_opt}" BaseRecalibrator \
-R ${REF} \
-I ${Sample}.markduplicates.bam \
--use-original-qualities \
-O ${Sample}.recalibration_table \
--known-sites ${dbSNP_vcf} \
--known-sites ${known_indels_sites_VCF1} \
--known-sites ${known_indels_sites_VCF2} \
-L ${INTERVAL} \
--tmp-dir $PWD \
--VERBOSITY ERROR

#ApplyBQSR
gatk --java-options "${java_opt}" ApplyBQSR \
-R ${REF} \
-I ${Sample}.markduplicates.bam \
-O ${Sample}.analysisready.bam \
-L ${INTERVAL} \
-bqsr ${Sample}.recalibration_table \
--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
--add-output-sam-program-record \
--create-output-bam-md5 \
--use-original-qualities \
--tmp-dir $PWD \
--VERBOSITY ERROR

samtools index -@ {thread} ${Sample}.analysisready.bam

#Haplotypecaller
gatk --java-options "${java_opt}" HaplotypeCaller \
-R ${REF} \
-I ${Sample}.analysisready.bam \
-L ${INTERVAL} \
-O ${Sample}.g.vcf.gz \
-bamout ${Sample}.bamout.bam \
-ERC GVCF \
--tmp-dir $PWD \
--VERBOSITY ERROR

rm -rf ${Sample}.markduplicates.bam
rm -rf ${Sample}.markduplicates.bam.bai
rm -rf ${Sample}.realibration_table
rm -rf ${Sample}.mark_duplicate_metrix.txt
rm -rf ${Sample}.markilluminaadapters_metrics.txt
rm -rf ${Sample}.markilluminaadapters.bam
rm -rf ${Sample}.mapped.bam






