nextflow.preview.dsl = 2

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}
params.motif = params.genome ? params.genomes[ params.genome ].motif ?: false : false
params.motif_annot = params.genome ? params.genomes[ params.genome ].motif_annot ?: false : false
params.tss = params.genome ? params.genomes[ params.genome ].tss ?: false : false
params.size = params.genome ? params.genomes[ params.genome ].size ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false

process preprocess {
    conda '/home/xfu/miniconda3/envs/chop'
    publishDir "$params.outdir/expression"
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
    publishDir "$params.outdir/deg", saveAs: { 
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

process get_deg_promoter {
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
        | join -1 2 -2 4 - <(sort -k4,4 $params.tss) -t\$'\t' \
        | awk '{OFS="\t";print \$3,\$4-500,\$5+500,\$2,\$6,\$7}' \
        > ${deg.baseName}.bed
    """
}

process get_all_promoter {
    output:
    file "*.bed"
    script:
    """
    cat $params.tss \
        | awk -F"[\t_]" '{OFS="\t";print \$1,\$2-500,\$3+500,\$4,\$6,\$7}' \
        > ${params.genome}.promoter.bed
    """
}

include fimo from './nf_modules/meme/fimo.nf' params (
    conda: '/home/xfu/miniconda3/envs/meme',
    motif: "$params.motif"
)

include ame from './nf_modules/meme/ame.nf' params (
    conda: '/home/xfu/miniconda3/envs/meme',
    motif: "$params.motif"
)

process find_coregulators {
  publishDir "${params.outdir}/coregulators"
  input:
  tuple val(name), file(fimo), val(target)
  output:
  path "*.csv"
  script:
  """
  find_coregulators.py -f $fimo -t $target -a $params.motif_annot -i 10 -o 40 -O ${name}.coregulators_of_${target}.csv
  """
}

include {get_fasta; get_fasta as get_fasta_1} \
    from './nf_modules/bedtools/get_fasta.nf' \
    params (fasta: "$params.fasta")

workflow {
    microarray = Channel
        .fromPath( "raw/microarray/" )
        .ifEmpty { exit 1, "Microarray file not found: ${params.microarray}" }

    promoter_fa = get_all_promoter | get_fasta

    microarray \
        | preprocess \
        | deg \
        | flatten \
        | filter ( ~/^.*.csv/ ) \
        | get_deg_promoter \
        | get_fasta_1 \
        | fimo
    
    get_fasta_1.out \
        | combine (promoter_fa) \
        | ame

    // invariant: co-regulators of Atf4, Sox9
    interested_tfs = Channel.fromPath('./interested_tfs.txt') \
        | splitCsv \
        | map { x -> x[0] } \

    fimo.out \
        | filter ( ~/^.*invariant.*/ ) \
        | combine (interested_tfs) \
        | find_coregulators
}
