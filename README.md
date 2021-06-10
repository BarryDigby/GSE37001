# M6A robustness test

<p markdown="1" align="center">
  <img src="assets/images/M6a_paper.png" alt="paper_header">
</p>

This repository contains the results of differential gene, isofrom, exon and intron usage in siRNA directed `METTL3` knockdown & control HepG2 cells (n=2) in an attempt to faithfully reproduce the RNA-Seq analysis performed in the Nature paper [Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq](https://www.nature.com/articles/nature11112).

Methods and results are detailed below for posterity, with comments describing missing information that would have aided the analysis from a reviewers perpective. The results from this repository are to be used in conjunction with a M6A peak calling analysis of HepG2 cells, to interrogate correlations between M6A peaks and differential gene/isoform/exon/intron usage, a topic that has 13 [contrasting citations](https://scite.ai/reports/topology-of-the-human-and-WDmMRO?contradicting=true&mentioning=false&page=1&supporting=false) to date.

## Analysis overview

* [Download data](#download-data)
* [Quantification](#quantification)
* [Differentially expressed genes](#differentially-expressed-genes)
* [Differentially expressed isoforms](#differentially-expressed-isoforms)
* [Differentially expressed exons](#differentially-expressed-exons)
* [Differentially expressed introns](#differentially-expressed-introns)

> Please note: The direction of fold change is with respect to control. I.e a log2FC value of 2 refers to up-regulation in the METTL3 knockdown samples.

# Download data

Raw sequencing data was downloaded using `SRAtools` `fastq-dump` via a singularity container. The nextflow script to download the reads is provided in [`scripts/`](https://github.com/BarryDigby/GSE37001/tree/main/scripts) and the set of commands used is given below:

<details markdown="1">
<summary>Download raw reads</summary>

```bash
singularity pull sratoolkit.img docker://pegi3s:sratoolkit
nextflow -bg run dl_sra.nf --sra_id 'SRP012096' -with-singularity 'sratoolkit.img'
```

</details>

SRA ID's can be mapped to the corresponding experiment metadata using the table below:


|       SRA ID       |   Experimental design  |
|:------------------:|:----------------------:|
| SRR456526.fastq.gz |   METTL3_KD1.fastq.gz  |
| SRR456527.fastq.gz |   METTL3_KD2.fastq.gz  |
| SRR456528.fastq.gz | Mock_control1.fastq.gz |
| SRR456529.fastq.gz | Mock_control2.fastq.gz |

***

Reference genome and GTF files were prepared as per the paper, using `H. sapiens ENSEMBL release 54 (NCBI36/hg18)`.

<details markdown="1">
<summary>Download reference files</summary>

```bash
wget http://ftp.ensembl.org/pub/release-54/fasta/homo_sapiens/dna/Homo_sapiens.NCBI36.54.dna.toplevel.fa.gz
guzip Homo_sapiens.NCBI36.54.dna.toplevel.fa.gz

wget http://ftp.ensembl.org/pub/release-54/gtf/homo_sapiens/Homo_sapiens.NCBI36.54.gtf.gz
gunzip Homo_sapiens.NCBI36.54.gtf.gz
```

</details>

# Quantification

 RNA-Seq analysis was performed using `nf-core/rnaseq v3.1` with default parameters (except for the reference files provided). The `samples.csv` file provided to `nf-core/rnaseq` is given below, of note the dataset is single-end and unstranded:

 | sample        | fastq_1                      | fastq_2 | strandedness |
 |---------------|------------------------------|---------|--------------|
 | METTL3_KD1    | fastq/METTL3_KD1.fastq.gz    |         | unstranded   |
 | METTL3_KD2    | fastq/METTL3_KD2.fastq.gz    |         | unstranded   |
 | Mock_control1 | fastq/Mock_control1.fastq.gz |         | unstranded   |
 | Mock_control2 | fastq/Mock_control2.fastq.gz |         | unstranded   |

***

 <details markdown="1">
 <summary>Nextflow command</summary>

 ```bash
nextflow pull nf-core/rnaseq
nextflow -bg run nf-core/rnaseq -profile singularity --input 'samples.csv' --fasta 'assets/Homo_sapiens.NCBI36.54.dna.toplevel.fa' --gtf 'Homo_sapiens.NCBI36.54.gtf' --max_memory '62.GB' --max_cpus 16
 ```

 </details>


# Differentially expressed genes

Pairwise comparisons of `knockdown` vs `control` populations were generated using an additive linear model with replicates as the blocking factor in `DESeq2`.

The `R` code used to generate results are given below:

<details markdown="1">
<summary>Differential gene analysis</summary>

```R
# stage files
meta <- data.frame(row.names=c("METTL3_KD1", "METTL3_KD2", "Mock_control1", "Mock_control2"),
                   "condition"=c("KD", "KD", "CTRL", "CTRL"),
                   "replicates"=as.factor(c("1", "2", "1", "2")))
dir <- "/data/projects/leipzig/results/star_salmon"
files <- file.path(dir, rownames(meta), "quant.sf")
# Use release 54 here *crucial for accurate coordinates across DE analyses
mart <- useMart("ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl",
                host="may2009.archive.ensembl.org",
                path="/biomart/martservice",
                archive=FALSE)
tx2gene <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"), mart = mart)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = FALSE)
# DDS object
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ replicates + condition )
dds$condition <- relevel(dds$condition, ref="CTRL")
dds <- DESeq(dds)
# DGE
res <- results(dds, filterFun=ihw, alpha=0.05, c("condition", "KD", "CTRL"))
LFC <- lfcShrink(dds = dds, res= res, coef = 3, type = "apeglm")
LFC_df <- as.data.frame(LFC)
# functions
# use paper cutoff here (FC > 2, FDR 5%)
get_upregulated <- function(df){

    key <- intersect(rownames(df)[which(df$log2FoldChange>=2)], rownames(df)[which(df$padj<=0.05)])

    results <- as.data.frame((df)[which(rownames(df) %in% key),])
    results <- results[order(-results$log2FoldChange),]
    return(results)
}
get_downregulated <- function(df){

    key <- intersect(rownames(df)[which(df$log2FoldChange<=-2)], rownames(df)[which(df$padj<=0.05)])

    results <- as.data.frame((df)[which(rownames(df) %in% key),])
    results <- results[order(results$log2FoldChange),]
    return(results)
}
annotate_de_genes <- function(df){

    df$hgnc_symbol <- rownames(df)
    mart <- useMart("ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl",
                host="may2009.archive.ensembl.org",
                path="/biomart/martservice",
                archive=FALSE)
    info <- getBM(attributes=c("hgnc_symbol",
                               "ensembl_gene_id",
                               "chromosome_name",
                               "start_position",
                               "end_position",
                               "strand"),
                  filters = c("hgnc_symbol"),
                  values = df$hgnc_symbol,
                  mart = mart,
                  useCache=FALSE)

    tmp <- merge(df, info, by="hgnc_symbol")
    tmp$strand <- gsub("-1", "-", tmp$strand)
    tmp$strand <- gsub("1", "+", tmp$strand)
    tmp$hgnc_symbol <- make.names(tmp$hgnc_symbol, unique = T)
    #tmp <- tmp[!grepl("CHR", tmp$chromosome_name),]

    output_col <- c("hgnc", "ensembl_gene_id", "chromosome", "start", "end", "strand", "log2FC", "pvalue", "padj")
    tmp <- subset(tmp, select=c(hgnc_symbol, ensembl_gene_id, chromosome_name, start_position, end_position, strand, log2FoldChange, pvalue, padj))
    colnames(tmp) <- output_col

    if(min(tmp$Log2FC) > 0){
        tmp <- tmp[order(-tmp$log2FC),]
    }else{
        tmp <- tmp[order(tmp$log2FC),]
    }

    return(tmp)

}
# get up regulated
up <- get_upregulated(LFC)
up$hgnc_symbol <- rownames(up)
up <- annotate_de_genes(up)
DT::datatable(up, rownames = FALSE)
# get down regulated
down <- get_downregulated(LFC)
down$hgnc_symbol <- rownames(down)
down <- annotate_de_genes(down)
DT::datatable(down, rownames = FALSE)
# write to file
write.table(up, "/data/github/GSE37001/gene/DESeq2_gene_upregulated.txt", sep="\t", quote = FALSE, row.names = FALSE)
write.table(down, "/data/github/GSE37001/gene/DESeq2_gene_downregulated.txt", sep="\t", quote = FALSE, row.names = FALSE)
```

</details>

### Comments

The number of differentially expressed genes returned by the analysis (223) were significantly lower than those reported by the study (1977), despite using the same cut-off values (LFC > 2 & FDR 5%). In my opinion, this is due to using `filterFun=ihw` when extracting the `DESeq2` results, and by using `apeglm` `LFCShrink` penalised regression to reduce low confidence DEGs.

Despite the discrepancy in results, this was a relatively simple analysis to perform. The paper stated which reference genome files were used (`ENSEMBL release 54`) which is crucial in returning the same genomic coordinates for M6a peak overlap analysis.

# Differentially expressed isoforms

Initially attempted the analysis using `Stringtie` output files in `Ballgown`, but was not satisfied with the results generated. Specifying `adjustVars="replicates"` in an attempt to construct an additive linear model with replicates as the blocking factor (as in the DEG analysis) produced `NAN` pvalue and adjusted pvalues. This is most likely due to a lack of variance amongst transcripts after correcting for `replicates` and due to the fact the study was underpowered.

To overcome this, the analysis was performed using the `DESeq2` workflow above with one change: `txi <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = TRUE)` to use transcript counts in the analysis instead of gene counts. The workflow is given below.

<details markdown="1">
<summary>Differential isoform analysis</summary>

```R
# same as above, but use TX counts
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = TRUE)
# DDS object
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ replicates + condition )
dds$condition <- relevel(dds$condition, ref="CTRL")
dds <- DESeq(dds)
# DGE
res <- results(dds, filterFun=ihw, alpha=0.05, c("condition", "KD", "CTRL"))
LFC <- lfcShrink(dds = dds, res= res, coef = 3, type = "apeglm")
LFC_df <- as.data.frame(LFC)
# functions
# use paper cutoff here (FC > 2, FDR 5%)
# use paper cutoff here (FC > 2, FDR 5%)
get_upregulated <- function(df){

    key <- intersect(rownames(df)[which(df$log2FoldChange>=2)], rownames(df)[which(df$padj<=0.05)])

    results <- as.data.frame((df)[which(rownames(df) %in% key),])
    results <- results[order(-results$log2FoldChange),]
    return(results)
}
get_downregulated <- function(df){

    key <- intersect(rownames(df)[which(df$log2FoldChange<=-2)], rownames(df)[which(df$padj<=0.05)])

    results <- as.data.frame((df)[which(rownames(df) %in% key),])
    results <- results[order(results$log2FoldChange),]
    return(results)
}
annotate_de_genes <- function(df){

    df$hgnc_symbol <- rownames(df)
    mart <- useMart("ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl",
                host="may2009.archive.ensembl.org",
                path="/biomart/martservice",
                archive=FALSE)
    info <- getBM(attributes=c("ensembl_transcript_id",
                               "ensembl_gene_id",
                               "chromosome_name",
                               "start_position",
                               "end_position",
                               "strand"),
                  filters = c("ensembl_transcript_id"),
                  values = df$ensembl_transcript_id,
                  mart = mart,
                  useCache=FALSE)

    tmp <- merge(df, info, by="ensembl_transcript_id")
    tmp$strand <- gsub("-1", "-", tmp$strand)
    tmp$strand <- gsub("1", "+", tmp$strand)
    tmp$ensembl_transcript_id <- make.names(tmp$ensembl_transcript_id, unique = T)
    #tmp <- tmp[!grepl("CHR", tmp$chromosome_name),]

    output_col <- c("ensembl_transcript_id", "ensembl_gene_id", "chromosome", "start", "end", "strand", "log2FC", "pvalue", "padj")
    tmp <- subset(tmp, select=c(ensembl_transcript_id, ensembl_gene_id, chromosome_name, start_position, end_position, strand, log2FoldChange, pvalue, padj))
    colnames(tmp) <- output_col

    if(min(tmp$Log2FC) > 0){
        tmp <- tmp[order(-tmp$log2FC),]
    }else{
        tmp <- tmp[order(tmp$log2FC),]
    }

    return(tmp)

}
# get up regulated
up <- get_upregulated(LFC)
up$ensembl_transcript_id <- rownames(up)
up <- annotate_de_genes(up)
DT::datatable(up, rownames = FALSE)
# get down regulated
down <- get_downregulated(LFC)
down$ensembl_transcript_id <- rownames(down)
down <- annotate_de_genes(down)
DT::datatable(down, rownames = FALSE)
# write to file
write.table(up, "/data/github/GSE37001/isoform/DESeq2_isoform_upregulated.txt", sep="\t", quote = FALSE, row.names = FALSE)
write.table(down, "/data/github/GSE37001/isoform/DESeq2_isoform_downregulated.txt", sep="\t", quote = FALSE, row.names = FALSE)
```

</details>

### Comments

Once again the number of differentially expressed isoforms returned by the analysis (408) is substantially lower than those reported by the paper (7521!). This is likely due to the independent filtering method, apeglm methods employed by the analysis.

Recreating this analysis was difficult, in the end I compromised by using `DESeq2`. Per the paper:

> Differentially expressed isoforms: Ensembl gtf file of all human genes (hg18 release 54) was re-processed using Cuffcompare v1.0.3 in order to add the missing tss_id and p_id attributes according to the user guide. The resulting gtf annotation file created by Cuffcompare was used as input to Cuffdiff v1.0.3 tool together with the fragment alignment files. Both Cuffcompare and Cuffdiff are part of the Cufflinks package.

The `Stringtie` and `Ballgown` tuxedo suite superseded `TopHat` and `Cufflinks` originally used in the study. Perhaps advanced users comfortable with the `Tuxedo` suite could faithfully reproduce the differentially expressed isoforms produced by the paper (e.g interpret and implement the `tss_id` and `p_id` values generated by the authors).

# Differentially expressed exons

The standard `DEXSeq` analysis workflow was followed to produce results for differentially expressed exons. Prior to analysis in R, reference GFF files and exon counts were generated. The code is given below.

<details markdown="1">
<summary>Prepare annotation</summary>

```bash
python /home/barry/R/x86_64-pc-linux-gnu-library/4.1/DEXSeq/python_scripts/dexseq_prepare_annotation.py Homo_sapiens.NCBI36.54.gtf Homo_sapiens.NCBI36.54.gff -r no
```

</details>

<details markdown="1">
<summary>Counting Reads</summary>

```bash
python /home/barry/R/x86_64-pc-linux-gnu-library/4.1/DEXSeq/python_scripts/dexseq_count.py Homo_sapiens.NCBI36.54.gff /data/projects/leipzig/results/star_salmon/METTL3_KD1.markdup.sorted.bam METTL3_KD1.txt -r pos -s no -f bam -a 0

python /home/barry/R/x86_64-pc-linux-gnu-library/4.1/DEXSeq/python_scripts/dexseq_count.py Homo_sapiens.NCBI36.54.gff /data/projects/leipzig/results/star_salmon/METTL3_KD2.markdup.sorted.bam METTL3_KD2.txt -r pos -s no -f bam -a 0

python /home/barry/R/x86_64-pc-linux-gnu-library/4.1/DEXSeq/python_scripts/dexseq_count.py Homo_sapiens.NCBI36.54.gff /data/projects/leipzig/results/star_salmon/Mock_control1.markdup.sorted.bam Mock_control1.txt -r pos -s no -f bam -a 0

python /home/barry/R/x86_64-pc-linux-gnu-library/4.1/DEXSeq/python_scripts/dexseq_count.py Homo_sapiens.NCBI36.54.gff /data/projects/leipzig/results/star_salmon/Mock_control2.markdup.sorted.bam Mock_control2.txt -r pos -s no -f bam -a 0
```

</details>

<details markdown="1">
<summary>DEXSeq</summary>

```R
# stage files
inDir = "/data/projects/leipzig/dexseq/"
countFiles = list.files(inDir, pattern=".txt$", full.names=TRUE)
# stage GFF
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
# define full, reduced models
# "To detect differences in exon usage that affect both replicates in the same manner due to condition"
formulaFullModel    =  ~ sample + exon + replicates:exon + condition:exon
formulaReducedModel =  ~ sample + exon + replicates:exon
# dxd object
dxd = DEXSeqDataSetFromHTSeq(
      countFiles,
      sampleData=meta,
      design= ~ sample + exon + replicates:exon + condition:exon,
      flattenedfile=flattenedFile )
# relevel, normalize
dxd$condition <- relevel(dxd$condition, ref="CTRL")
dxd <- estimateSizeFactors(dxd)
# gimme power!
BPPARAM = BiocParallel::MulticoreParam(4)
dxd <- estimateDispersions(dxd, formula = formulaFullModel, BPPARAM=BPPARAM)
dxd <- testForDEU(dxd, fullModel = formulaFullModel, reducedModel = formulaReducedModel, BPPARAM = BPPARAM)
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM = BPPARAM, independentFiltering = TRUE)
dex_res <- DEXSeqResults(dxd)
dex_df <- as.data.frame(dex_res)
# subset results
get_upregulated_dex <- function(df){

    key <- intersect(rownames(df)[which(df$log2fold_KD_CTRL>=2)], rownames(df)[which(df$padj<=0.05)])

    results <- as.data.frame((df)[which(rownames(df) %in% key),])
    results <- results[order(-results$log2fold_KD_CTRL),]
    return(results)
}
get_downregulated_dex <- function(df){

    key <- intersect(rownames(df)[which(df$log2fold_KD_CTRL<=-2)], rownames(df)[which(df$padj<=0.05)])

    results <- as.data.frame((df)[which(rownames(df) %in% key),])
    results <- results[order(results$log2fold_KD_CTRL),]
    return(results)
}
# run functions
up_dex <- get_upregulated_dex(dex_df)
down_dex <- get_downregulated_dex(dex_df)
# tidy
up_dex <- subset(up_dex, select=c(groupID, transcripts, featureID, genomicData.seqnames, genomicData.start, genomicData.end, genomicData.strand, log2fold_KD_CTRL, pvalue, padj))
colnames(up_dex) <- c("ensembl_gene_id", "ensembl_transcript_id", "exon_id", "chromosome", "start", "end", "strand", "log2FC", "pvalue", "padj" )
down_dex <- subset(down_dex, select=c(groupID, transcripts, featureID, genomicData.seqnames, genomicData.start, genomicData.end, genomicData.strand, log2fold_KD_CTRL, pvalue, padj))
colnames(down_dex) <- c("ensembl_gene_id", "ensembl_transcript_id", "exon_id", "chromosome", "start", "end", "strand", "log2FC", "pvalue", "padj" )
# write to rda file (transcripts column is a list of transcripts overlapping the exon, very hard to read into R by TSV/CSV)
saveRDS(up_dex, file="/data/github/GSE37001/exon/DEXSeq_exons_upregulated.rda")
saveRDS(down_dex, file="/data/github/GSE37001/exon/DEXSeq_exons_downregulated.rda")
```

### Comments

The analysis yielded less differentially expressed exons (126) compared to the paper (474). It is less clear why this is the case, as `IHW` nor `apeglm` filtering were applied to the workflow.

The analysis was straight forward, the information provided by the paper was sufficient to perfrom the analysis.

</details>

# Differentially expressed introns

The methods section of the paper states they used a customised in-house script to prepare introns for DESeq. I have decided to reformat the reference GTF file to contain introns (not exons) and use this as the starting point for a `DEXSeq` analysis. In theory, this should work the exact same, producing differentially expressed introns as outputs of the analysis. Any criticisms of the approach taken with the analysis are welcome!

<details markdown="1">
<summary>Extract introns from GTF</summary>

```R
library(GenomicFeatures)
library(rtracklayer)
txdb <- makeTxDbFromGFF('/data/projects/leipzig/introns/Homo_sapiens.NCBI36.54.gtf')
introns <- intronicParts(txdb)
rtracklayer::export(introns, "/data/projects/leipzig/introns/introns.gtf")
```

</details>

Next, using a customised version of `DEXSeq` `prepare_annotations.py` (available in [`scripts/`](https://github.com/BarryDigby/GSE37001/tree/main/scripts)), produce a GFF file containing non-overlapping introns for `DEXSeq` analysis. (The results will follow the naming convention 'exon', thus they have been edited within the R script to output 'intron').

<details markdown="1">
<summary>Prepare (intron) annotations</summary>

```bash
sed 's/tx_name/transcript_id/g' introns.gtf > introns_rename.gtf
python prepare_annotation_introns.py introns_rename.gtf introns.gff -r no
```

</details>

Finally, use the sequencing BAM files in conjunction with `dexseq_count.py` to produce counts for each intron.

<details markdown="1">
<summary>Counting Reads</summary>

```bash
python /home/barry/R/x86_64-pc-linux-gnu-library/4.1/DEXSeq/python_scripts/dexseq_count.py introns.gff /data/projects/leipzig/results/star_salmon/METTL3_KD1.markdup.sorted.bam METTL3_KD1.txt -r pos -s no -f bam -a 0

python /home/barry/R/x86_64-pc-linux-gnu-library/4.1/DEXSeq/python_scripts/dexseq_count.py introns.gff /data/projects/leipzig/results/star_salmon/METTL3_KD2.markdup.sorted.bam METTL3_KD2.txt -r pos -s no -f bam -a 0

python /home/barry/R/x86_64-pc-linux-gnu-library/4.1/DEXSeq/python_scripts/dexseq_count.py introns.gff /data/projects/leipzig/results/star_salmon/Mock_control1.markdup.sorted.bam Mock_control1.txt -r pos -s no -f bam -a 0

python /home/barry/R/x86_64-pc-linux-gnu-library/4.1/DEXSeq/python_scripts/dexseq_count.py introns.gff /data/projects/leipzig/results/star_salmon/Mock_control2.markdup.sorted.bam Mock_control2.txt -r pos -s no -f bam -a 0
```

</details>


<details markdown="1">
<summary>DEXSeq (introns)</summary>

```R
# stage files
inDir = "/data/projects/leipzig/introns/"
countFiles = list.files(inDir, pattern=".txt$", full.names=TRUE)
# stage GFF
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
# define full, reduced models
# "To detect differences in exon usage that affect both replicates in the same manner due to condition"
formulaFullModel    =  ~ sample + exon + replicates:exon + condition:exon
formulaReducedModel =  ~ sample + exon + replicates:exon
# dxd object
dxd = DEXSeqDataSetFromHTSeq(
      countFiles,
      sampleData=meta,
      design= ~ sample + exon + replicates:exon + condition:exon,
      flattenedfile=flattenedFile)
# relevel, normalize
dxd$condition <- relevel(dxd$condition, ref="CTRL")
dxd <- estimateSizeFactors(dxd)
# gimme power!
BPPARAM = BiocParallel::MulticoreParam(4)
dxd <- estimateDispersions(dxd, formula = formulaFullModel, BPPARAM=BPPARAM)
plotDispEsts(dxd)
dxd <- testForDEU(dxd, fullModel = formulaFullModel, reducedModel = formulaReducedModel, BPPARAM = BPPARAM)
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM = BPPARAM, independentFiltering = TRUE)
dex_res <- DEXSeqResults(dxd)
dex_df <- as.data.frame(dex_res)
# subset results
get_upregulated_dex <- function(df){

    key <- intersect(rownames(df)[which(df$log2fold_KD_CTRL>=2)], rownames(df)[which(df$padj<=0.05)])

    results <- as.data.frame((df)[which(rownames(df) %in% key),])
    results <- results[order(-results$log2fold_KD_CTRL),]
    return(results)
}
get_downregulated_dex <- function(df){

    key <- intersect(rownames(df)[which(df$log2fold_KD_CTRL<=-2)], rownames(df)[which(df$padj<=0.05)])

    results <- as.data.frame((df)[which(rownames(df) %in% key),])
    results <- results[order(results$log2fold_KD_CTRL),]
    return(results)
}
# run functions
up_dex <- get_upregulated_dex(dex_df)
down_dex <- get_downregulated_dex(dex_df)
# tidy
up_dex <- subset(up_dex, select=c(groupID, transcripts, featureID, genomicData.seqnames, genomicData.start, genomicData.end, genomicData.strand, log2fold_KD_CTRL, pvalue, padj))
colnames(up_dex) <- c("ensembl_gene_id", "ensembl_transcript_id", "intron_id", "chromosome", "start", "end", "strand", "log2FC", "pvalue", "padj" )
down_dex <- subset(down_dex, select=c(groupID, transcripts, featureID, genomicData.seqnames, genomicData.start, genomicData.end, genomicData.strand, log2fold_KD_CTRL, pvalue, padj))
colnames(down_dex) <- c("ensembl_gene_id", "ensembl_transcript_id", "intron_id", "chromosome", "start", "end", "strand", "log2FC", "pvalue", "padj" )
# write to rda file (transcripts column is a list of transcripts overlapping the exon, very hard to read into R by TSV/CSV)
saveRDS(up_dex, file="/data/github/GSE37001/intron/DEXSeq_introns_upregulated.rda")
saveRDS(down_dex, file="/data/github/GSE37001/intron/DEXSeq_introns_downregulated.rda")
```

</details>
