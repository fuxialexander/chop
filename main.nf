nextflow.preview.dsl = 2

process preprocess {
    conda '/home/xfu/miniconda3/envs/chop'
    publishDir "${params.outdir}/expression/"
    input:
    file microarray

    output:
    file "exprs.ann.txt"
    script:
    """
    array.R $microarray
    """
}

process deg {
    conda '/home/xfu/miniconda3/envs/chop'
    publishDir "${params.outdir}/deg", saveAs: { 
        filename -> filename.endsWith(".pdf") 
        ? "plots/$filename" 
        : "lists/$filename"
    }
    input:
    file expression
    output:
    file "*"
    script:
    """
    deg.py
    """
}

process get_promoter {
    input:
    file deg
    output:
    file "*.bed"
    script:
    """
    cut -d, -f1 $deg \
        | tail -n +2 \
        | awk '{print \$1"\t"\$1"_1"}'  \
        | sort -k2,2 \
        | join -1 2 -2 4 - <(sort -k4,4 ${params.genomes.ucsc_mm10.tss}) -t\$'\t' \
        | awk '{OFS="\t";print \$3,\$4-500,\$5+500,\$2,\$6,\$7}' \
        > ${deg.baseName}.bed
    """
}

process get_fasta {
    input:
    file promoter
    output:
    file "*.fa"
    script:
    """
    bedtools getfasta -name -fi ${params.genomes.ucsc_mm10.fasta} -bed $promoter -fo ${promoter.baseName}.fa
    """
}

process fimo {
    conda '/home/xfu/miniconda3/envs/meme'
    publishDir "${params.outdir}/fimo", saveAs: { filename -> "${fasta.baseName}.fimo.tsv" }
    input:
    file fasta
    output:
    file("fimo_out/fimo.tsv")
    script:
    """
    fimo ${params.genomes.ucsc_mm10.motif} $fasta
    """
}

process ame {
    conda '/home/xfu/miniconda3/envs/meme'
    publishDir "${params.outdir}/ame", saveAs: { filename -> "${fasta.baseName}.ame.tsv" }
    input:
    file fasta
    output:
    file("ame_out/ame.tsv")
    script:
    """
    ame $fasta ${params.genomes.ucsc_mm10.motif}
    """
}

workflow {
    microarray = Channel
        .fromPath( "raw/microarray/" )
        .ifEmpty { exit 1, "Microarray file not found: ${params.microarray}" }

    microarray \
        | preprocess \
        | deg \
        | flatten \
        | filter ( ~/^.*.csv/ ) \
        | get_promoter \
        | get_fasta \
        | ( fimo & ame )
}
