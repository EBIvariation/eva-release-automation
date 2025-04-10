#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.active_chunk_size = 10000000
params.merged_chunk_size = 100000
params.deprecated_chunk_size = 100000
workflow {

    initiate_release_status_for_assembly()
    run_dump_active_rs_for_assembly(initiate_release_status_for_assembly.out.flag)
    split_release_active_for_assembly(run_dump_active_rs_for_assembly.out.release_active_rs, params.active_chunk_size)
    release_active_rs_for_assembly(split_release_active_for_assembly.out.release_active_chunks)
    sort_and_index_chunk_active(release_active_rs_for_assembly.out.release_active_chunk)
    merge_active_chunks(sort_and_index_chunk_active.out.sorted_release_active_chunk.collect(), sort_and_index_chunk_active.out.index_release_active_chunk.collect())

    run_dump_merged_rs_for_assembly(initiate_release_status_for_assembly.out.flag)
    split_release_merged_for_assembly(run_dump_merged_rs_for_assembly.out.release_merged_rs, params.merged_chunk_size)
    release_merged_rs_for_assembly(split_release_active_for_assembly.out.release_merged_chunks)
    sort_and_index_chunk_merged(release_merged_rs_for_assembly.out.release_merged_chunk)
    merge_merged_chunks(sort_and_index_chunk_merged.out.sorted_release_merged_chunk.collect(), sort_and_index_chunk_merged.out.index_release_merged_chunk.collect())


    run_dump_deprecated_rs_for_assembly(initiate_release_status_for_assembly.out.flag)
    split_release_deprecated_for_assembly(run_dump_deprecated_rs_for_assembly.out.release_deprecated_rs, params.deprecated_chunk_size)
    release_deprecated_rs_for_assembly(split_release_deprecated_for_assembly.out.release_deprecated_chunks)
    merge_deprecated_chunks(release_deprecated_rs_for_assembly.out.release_deprecated_chunk.collect(), release_merged_rs_for_assembly.out.release_merged_deprecated_chunk.collect())

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
    path "active-rs-chunk-*", emit: release_active_chunks

    script:
    """
    split -a 5 -d -l ${chunk_size} release_active_dump active-rs-chunk-
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


process sort_and_index_chunk_active {

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


process run_dump_merged_rs_for_assembly {

    label 'long_time', 'med_mem'

    input:
    val flag

    output:
    path "release_merge_dump", emit: release_merged_rs

    script:
    def pipeline_parameters = " --spring.batch.job.names=DUMP_MERGED_ACCESSIONS_JOB"
    pipeline_parameters += " --parameters.outputFolder=\$PWD"
    pipeline_parameters += " --parameters.rsAccDumpFile=" + "release_merged_dump"
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.release_pipeline --spring.config.location=file:$params.release_job_props $pipeline_parameters
    """
}

process split_release_merged_for_assembly {

    label 'long_time', 'med_mem'

    input:
    path "release_merged_dump"
    val chunk_size

    output:
    path "merged-rs-chunk-*", emit: release_merged_chunks

    script:
    """
    split -a 5 -d -l ${chunk_size} release_merged_dump merged-rs-chunk-
    """
}

process release_merged_rs_for_assembly {

    label 'med_time', 'med_mem'

    input:
    each path(rs_chunk)

    output:
    path "merged_ids_${rs_chunk}.vcf", emit: release_merged_chunk
    path "merged_deprecated_ids_${rs_chunk}.txt",  emit: release_merged_deprecated_chunk

    script:
    def pipeline_parameters = " --spring.batch.job.names=MERGED_ACCESSIONS_RELEASE_FROM_DB_JOB"
    pipeline_parameters += " --parameters.rsAccFile=" + rs_chunk
    pipeline_parameters += " --parameters.outputFolder=\$PWD"
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.release_pipeline --spring.config.location=file:$params.release_job_props $pipeline_parameters
    mv *_merged_ids.vcf merged_ids_${rs_chunk}.vcf
    mv *_deprecated_ids.unsorted.txt merged_deprecated_ids_${rs_chunk}.txt
    """
}


process sort_and_index_chunk_merged {

    label 'med_time', 'med_mem'

    input:
    path(release_merged_chunk)

    output:
    path "*_sorted.vcf.gz", emit: sorted_release_merged_chunk
    path "*_sorted.vcf.gz.csi", emit: index_release_merged_chunk

    script:
    """
    awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' $rrelease_merged_chunk | $params.executable.bcftools view --no-version  -O z -o ${release_merged_chunk}_sorted.vcf.gz -
    $params.executable.bcftools index ${release_active_chunk}_sorted.vcf.gz
    """
}


process merge_merged_chunks {

    label 'med_time', 'med_mem'

    input:
    path(release_merged_chunks)
    path(index_release_merged_chunk)

    output:
    path "merged.vcf.gz", emit: release_merged_merged

    script:
    """
    $params.executable.bcftools concat -a -o merged.vcf.gz -O z $release_merged_chunks
    """
}


process run_dump_deprecated_rs_for_assembly {

    label 'long_time', 'med_mem'

    input:
    val flag

    output:
    path "release_deprecated_dump", emit: release_deprecated_rs

    script:
    def pipeline_parameters = " --spring.batch.job.names=DUMP_DEPRECATED_ACCESSIONS_JOB"
    pipeline_parameters += " --parameters.outputFolder=\$PWD"
    pipeline_parameters += " --parameters.rsAccDumpFile=" + "release_deprecated_dump"
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.release_pipeline --spring.config.location=file:$params.release_job_props $pipeline_parameters
    """
}

process split_release_deprecated_for_assembly {

    label 'long_time', 'med_mem'

    input:
    path "release_deprecated_dump"
    val chunk_size

    output:
    path "deprecated-rs-chunk-*", emit: release_deprecated_chunks

    script:
    """
    split -a 5 -d -l ${chunk_size} release_deprecated_dump deprecated-rs-chunk-
    """
}



process release_deprecated_rs_for_assembly {

    label 'med_time', 'med_mem'

    input:
    each path(rs_chunk)

    output:
    path "deprecated_ids_${rs_chunk}.vcf", emit: release_deprecated_chunk

    script:
    def pipeline_parameters = " --spring.batch.job.names=DEPRECATED_ACCESSIONS_RELEASE_FROM_DB_JOB"
    pipeline_parameters += " --parameters.rsAccFile=" + rs_chunk
    pipeline_parameters += " --parameters.outputFolder=\$PWD"
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.release_pipeline --spring.config.location=file:$params.release_job_props $pipeline_parameters
    mv *_deprecated_ids.unsorted.txt deprecated_ids_${rs_chunk}.txt
    """
}


process merge_deprecated_chunks {

    label 'med_time', 'med_mem'

    input:
    path(release_deprecated_chunks)
    path(release_merged_deprecated_chunks)


    output:
    path "deprecated.txt", emit: release_merged_deprecated

    script:
    """
    cat $release_deprecated_chunks $release_merged_deprecated_chunks | sort -u  > deprecated_rs.txt
    """
}
