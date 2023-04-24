#! /usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// Log info
log.info "\n\n          S I P O R E  ~  version 0.1          "
log.info "==============================================="
log.info "raw reads folder                : ${params.raw_reads_folder}"
log.info "demultiplexed reads folder      : ${params.demulti_folder}"
log.info "trimmed reads folder            : ${params.trim_folder}"
log.info "\n"


/* Help message (from https://bioinformaticsworkbook.org/dataAnalysis/nextflow/02_creatingAworkflow.html#gsc.tab=0)
def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --query QUERY.fasta --dbDir "blastDatabaseDirectory" --dbName "blastPrefixName"

        Mandatory arguments:
         --query                        Query fasta file of sequences you wish to BLAST
         --dbDir                        BLAST database directory (full path required)
         --dbName                       Prefix name of the BLAST database

       Optional arguments:
        --outdir                       Output directory to place final BLAST output
        --outfmt                       Output format ['6']
        --options                      Additional options for BLAST command [-evalue 1e-3]
        --outFileName                  Prefix name for BLAST output [input.blastout]
        --threads                      Number of CPUs to use during blast job [16]
        --chunkSize                    Number of fasta records to use when splitting the query fasta file
        --app                          BLAST program to use [blastn;blastp,tblastn,blastx]
        --help                         This usage statement.
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}
*/


// Porechop
process PORECHOP {
    publishDir params.demulti_folder, mode: "copy"
    conda 'bioconda::porechop=0.2.4'

    input:
        path fastq

    output:
        path "demultiplexed_reads/*.fastq"
        path "demultiplexed_reads/none.fastq"

    script:
    """
    porechop -i $fastq -b $params.demulti_folder
    """
}

// Filtlong
process FILTLONG {
    publishDir params.trim_folder, mode: "copy"
    conda 'bioconda::filtlong=0.2.1'
    
    input:
        path demulti_reads
        path trim_folder

    output:
        path "*"

    script:
    """
    filtlong --min_length 1000 $demulti_reads > $params.trim_folder
    """

}


workflow {
    fastq_files_ch = Channel.fromPath(params.raw_reads_folder)
    demulti_reads_ch = PORECHOP(fastq_files_ch)
    trimmed_reads_ch = FILTLONG(demulti_reads_ch) // .flatten.collect() gives error
}


/*
workflow {
    PORECHOP | flatten | FILTLONG   // if processes were just bash commands
}
*/