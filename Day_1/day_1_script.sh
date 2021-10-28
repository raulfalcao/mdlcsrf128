
# Passo 1: Ao baixar os arquivos fornecidos a primeira coisa a ser feita é a:
# 1) descompactação da sequência referência
# 2) criação de um dicionário bem como o índice usado pelo alinhador bwa e pelo samtools


bgzip -d grch38.chr22.fasta.gz

samtools faidx grch38.chr22.fasta

java -jar gatk-package-4.2.0.0-local.jar CreateSequenceDictionary --TMP_DIR . -R grch38.chr22.fasta -O grch38.chr22.dict

bwa index grch38.chr22.fasta

# Passo 2: Controle de Qualidade com o fastqc e identificar possíveis artefatos

fastqc amostra-lbb_R* -o .


# Passo 3: Na etapa de Controle de Qualidade foi identificado sequências de adaptadores.
# Por essa razão será feito a eliminação dos adaptadores usando a ferramenta Trim Galore.
# Em seguida será feita o Controle de Qualidade novamente para termos certeza da eliminição
# dos adaptadores.


trim_galore --path_to_cutadapt ~/Library/Python/3.9/bin/cutadapt --paired --gzip amostra-lbb_R1.fq.gz amostra-lbb_R2.fq.gz

fastqc amostra-lbb_R*val* -o .

# Passo 4: Alinhar as sequências resultantes da etapa anterior com a ferramenta BWA.
# A etapa de alinhamento foi combinada com um filtro do samtools pela flag 2828. Ou seja,
# foram filtrados os reads que são:
# 1) read unmapped
# 2) mate unmapped
# 3) not primary alignment
# 4) read fails platform/vendor quality checks
# 5) supplementary alignment  

bwa mem -t 40 -M -R '@RG\tID:amostra-lbb\tSM:amostra-lbb\tPL:ILLUMINA' grch38.chr22.fasta amostra-lbb_R1_val_1.fq.gz amostra-lbb_R2_val_2.fq.gz | samtools view -bS -h -F 2828 - > amostra-lbb.filtered.bam

# Passo 5: Os comandos a seguir fazem parte das melhores práticas do pipeline do GATK 4.2.
# 1) Ordenar o arquivo bam;
# 2) Indexar o arquivo bam ordenado;
# 3) Marcar reads duplicados usando Picard;
# 4) Recalibrar as chamadas de bases usando a base de dados de SNPs de alta confiança do 1000 Genomas bem como a base de dados de INDELs Mills e 1000 Genomas;
# 5) Aplicar o método BQSR para a recalibração
# 6) Chamar as variantes com HaplotypeCaller passando como argumento a base de dados do dbSNP152. O output nessa etapa será um gvcf
# 7) Gerar o arquivo vcf a partir do gvcf

java -XX:ParallelGCThreads=16 -Xmx8G -jar gatk-package-4.2.0.0-local.jar SortSam I=amostra-lbb.filtered.bam O=amostra-lbb.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=.

java -XX:ParallelGCThreads=16 -Xmx8G -jar gatk-package-4.2.0.0-local.jar BuildBamIndex I=amostra-lbb.sorted.bam O=amostra-lbb.sorted.bam.bai TMP_DIR=.

java -XX:ParallelGCThreads=16 -Xmx8G -jar picard.jar MarkDuplicates TMP_DIR=. I=amostra-lbb.sorted.bam O=amostra-lbb.duplicates_marked.bam M=amostra-lbb.duplicates_marked.metrics.txt

java -XX:ParallelGCThreads=16 -Xmx8G -jar gatk-package-4.2.0.0-local.jar BaseRecalibrator --tmp-dir . -R  grch38.chr22.fasta -I amostra-lbb.duplicates_marked.bam  -O amostra-lbb.aligned.BaseRecalReport.grp --known-sites ~/database/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites ~/database/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

java -Xmx8G -jar gatk-package-4.2.0.0-local.jar ApplyBQSRSpark --tmp-dir . -R grch38.chr22.fasta -I amostra-lbb.duplicates_marked.bam  --bqsr-recal-file amostra-lbb.aligned.BaseRecalReport.grp -O amostra-lbb.processed.bam

java -XX:ParallelGCThreads=16 -Xmx8G -jar gatk-package-4.2.0.0-local.jar HaplotypeCaller --tmp-dir . -R grch38.chr22.fasta -I amostra-lbb.processed.bam --dbsnp ~/database/dbsnp152_GCF_000001405.38.vcf.gz -ERC GVCF -O amostra-lbb.gvcf.gz

java -XX:ParallelGCThreads=16 -Xmx8G -jar gatk-package-4.2.0.0-local.jar GenotypeGVCFs  -V amostra-lbb.gvcf.gz --tmp-dir . -R grch38.chr22.fasta -O amostra-lbb.vcf.gz

# Passo 6: As boas práticas do pipeline GATK sugere em analisar separadamente os SNPs e os INDELs. Dessa forma o arquivo vcf será dividido em 2 arquivos vcfs.

java -jar gatk-package-4.2.0.0-local.jar SelectVariants -V amostra-lbb.vcf.gz --tmp-dir . -select-type SNP -O amostra-lbb_SNPs.vcf.gz

java -jar gatk-package-4.2.0.0-local.jar SelectVariants -V amostra-lbb.vcf.gz --tmp-dir . -select-type INDEL -O amostra-lbb_INDELs.vcf.gz


# Passo 7: Usar filtros recomendados pelo GATK tanto para SNPs quanto para INDELs

java -jar gatk-package-4.2.0.0-local.jar VariantFiltration -V amostra-lbb_SNPs.vcf.gz --tmp-dir . -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O amostra-lbb_SNPs_filtered.vcf.gz

java -jar gatk-package-4.2.0.0-local.jar VariantFiltration -V amostra-lbb_INDELs.vcf.gz --tmp-dir . -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O amostra-lbb_INDELs_filtered.vcf.gz

# Passo 8: Concatenar os arquivos de SNPs e INDELs. Em seguida ordenar o arquivo vcf.

vcf-concat amostra-lbb_SNPs_filtered.vcf.gz amostra-lbb_INDELs_filtered.vcf.gz > amostra-lbb_filtered.vcf

vcf-sort amostra-lbb_filtered.vcf > amostra-lbb_filtered_sorted.vcf

# Passo 9: Quebrar em múltiplas linhas as possíveis variantes multi-alélicas

vt decompose -s -o amostra-lbb_decomposed.vcf amostra-lbb_filtered_sorted.vcf

vt normalize -r grch38.chr22.fasta -o amostra-lbb_normalized.vcf amostra-lbb_decomposed.vcf

# Passo 10: Compactar e indexar o resultado final. Então, finalmente o arquivo vcf.gz junto com o script será submetido para análise da Tarefa 1.

bgzip -c amostra-lbb_normalized.vcf > amostra-lbb_normalized.vcf.gz

tabix -p vcf amostra-lbb_normalized.vcf.vcf


