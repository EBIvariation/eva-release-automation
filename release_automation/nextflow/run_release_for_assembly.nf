#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.active_chunk_size = 10000000
workflow {

    initiate_release_status_for_assembly()
    run_dump_active_rs_for_assembly(initiate_release_status_for_assembly.out.flag)
    split_release_active_for_assembly(run_dump_active_rs_for_assembly.out.release_active_rs, params.active_chunk_size)
    release_active_rs_for_assembly(split_release_active_for_assembly.out.release_active_chunks)
    sort_and_index_chunk(release_active_rs_for_assembly.out.release_active_chunk)
    merge_active_chunks(sort_and_index_chunk.out.sorted_release_active_chunk.collect(), sort_and_index_chunk.out.index_release_active_chunk.collect())
}

process initiate_release_status_for_assembly {

    label 'short_time', 'med_mem'

    output:
    val true, emit: flag

    script:
    """
    $params.executable.python_interpreter -m release_automation.initiate_release_status_for_assembly --taxonomy-id $params.taxonomy --assembly-accession $params.assembly --release-version $params.release_version 1>> $params.log_file 2>&1
    """
}

process run_dump_active_rs_for_assembly {

    label 'long_time', 'med_mem'

    input:
    val flag

    output:
    path "release_active_dump", emit: release_active_rs

    script:
    def pipeline_parameters = " --spring.batch.job.names=DUMP_ACTIVE_ACCESSIONS_JOB"
    pipeline_parameters += " --parameters.outputFolder=\$PWD"
    pipeline_parameters += " --parameters.rsAccDumpFile=" + "release_active_dump"
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.release_pipeline --spring.config.location=file:$params.release_job_props $pipeline_parameters
    """
}

process split_release_active_for_assembly {

    label 'long_time', 'med_mem'

    input:
    path "release_active_dump"
    val chunk_size

    output:
    path "rs-chunk-*", emit: release_active_chunks

    script:
    """
    split -a 5 -d -l ${chunk_size} release_active_dump rs-chunk-
    """
}

process release_active_rs_for_assembly {

    label 'med_time', 'med_mem'

    input:
    each path(rs_chunk)

    output:
    path "current_ids_${rs_chunk}.vcf", emit: release_active_chunk

    script:
    def pipeline_parameters = " --spring.batch.job.names=ACTIVE_ACCESSIONS_RELEASE_FROM_DB_JOB"
    pipeline_parameters += " --parameters.rsAccFile=" + rs_chunk
    pipeline_parameters += " --parameters.outputFolder=\$PWD"
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.release_pipeline --spring.config.location=file:$params.release_job_props $pipeline_parameters
    mv *_current_ids.vcf current_ids_${rs_chunk}.vcf
    """
}


process sort_and_index_chunk {

    label 'med_time', 'med_mem'

    input:
    path(release_active_chunk)

    output:
    path "*_sorted.vcf.gz", emit: sorted_release_active_chunk
    path "*_sorted.vcf.gz.csi", emit: index_release_active_chunk

    script:
    """
    awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' $release_active_chunk | $params.executable.bcftools view --no-version  -O z -o ${release_active_chunk}_sorted.vcf.gz -
    $params.executable.bcftools index ${release_active_chunk}_sorted.vcf.gz
    """
}


process merge_active_chunks {

    label 'med_time', 'med_mem'

    input:
    path(release_active_chunks)
    path(index_release_active_chunk)

    output:
    path "active.vcf.gz", emit: release_active_merged

    script:
    """
    $params.executable.bcftools concat -a -o active.vcf.gz -O z $release_active_chunks
    """
}