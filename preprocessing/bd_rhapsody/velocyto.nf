#!/usr/bin/env nextflow

/**
 * Workflow to run velocyto.py on BD Rhapsody BAM files.
 *
 * The BAM files are maninpulated such that they are compatible with velocyto.
 * The respective code snipped has been provided by the BD Rhapsody support team.
 *
 * USAGE:
 *   - adjust all required parameters in `nf_velocyto.config`
 *   - execute the workflow with `nextflow run velocyto.py -c nf_velocyto.config`
 */

nextflow.enable.dsl = 2

process prepare_bam {

    input:
    path(bam)

    output:
    path("*.for_velocyto.bam"), emit: bam

    script:
    """
    cat \\
        <(samtools view -HS $bam) \\
        <(samtools view $bam | grep "MA:Z:*"  | sed  "s/MA:Z:/UB:Z:/" ) | \\
    samtools view -Sb -@6 > ${bam.baseName}.for_velocyto.bam
    """
}

process velocyto {
    // Velocyto otherwise resolves the symbolic link and writes into the original directory - need to make
    // a copy to stay clean.
    stageInMode 'copy'

    input:
    path bam
    path gtf
    path repeat_mask

    output:
    path "velocyto/*.loom", emit: loom

    script:
    """
    # uncompress gtf (works also if it is already uncompressed)
    gzip -cdf $gtf > GTFFILE.gtf

    # run velocyto
    velocyto run \\
        -m ${repeat_mask} \\
        --samtools-threads ${task.cpus} \\
        --samtools-memory ${task.memory.toMega()} \\
        $bam \\
        GTFFILE.gtf
    """
}

workflow {
    prepare_bam(Channel.fromPath(params.bam_files))

    velocyto(
        prepare_bam.out.bam,
        file(params.gtf, checkIfExists: true),
        file(params.repeat_mask, checkIfExists: true)
    )
}

