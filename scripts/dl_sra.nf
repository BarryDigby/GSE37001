#!/usr/bin/env nextflow

/*
 * Data is single end
 */


Channel
     .fromSRA(params.sra_id)
     .map{ it -> it[0] }
     .view()
     .set{ SRA_ID }


process dl{

  publishDir "./fastq", mode:'copy'

  input:
  val(base) from SRA_ID

  output:
  file("*") into fastqs

  script:
  """
  fastq-dump --gzip ${base}
  """
}
