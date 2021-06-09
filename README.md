# M6A robustness test

<p markdown="1" align="center">
  <img src="assets/images/M6a_paper.png" alt="paper_header">
</p>

This repository contains the results of Differential Gene, Transcript, Exon and Intron usage in siRNA directed `METTL3` knockdown & control HepG2 cells (n=2) in an attempt to faithfully reproduce the RNA-Seq analysis performed in the Nature paper [Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq](https://www.nature.com/articles/nature11112).

Methods and results are detailed below, with comments describing missing information that would have aided the analysis from a reviewers perpective. The results from this repository are to be used in conjunction with a M6A analysis of HepG2 cells, to interrogate correlations between M6A peaks and differential gene/transcript/exon/intron usage.

## Analysis overview

* [Quantification](#quantification)
* [Differentially expressed genes](#differentially-expressed-genes)
* [Differentially expressed isoforms](#differentially-expressed-isoforms)
* [Differentially expressed exons](#differentially-expressed-exons)
* [Differentially expressed introns](#differentially-expressed-introns)

# Quantification

# Differentially expressed genes

# Differentially expressed isoforms

# Differentially expressed exons

<details markdown="1">
<summary>Prepare annotation</summary>
```bash
python /home/barry/R/x86_64-pc-linux-gnu-library/4.1/DEXSeq/python_scripts/dexseq_prepare_annotation.py Homo_sapiens.NCBI36.54.gtf Homo_sapiens.NCBI36.54.gff -r no
```
</details>

# Differentially expressed introns
