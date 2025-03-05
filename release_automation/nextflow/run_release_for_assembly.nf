#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
    initiate_release_status_for_assembly('initiate') |
    run_dump_active_rs_for_assembly(initiate_release_status_for_assembly.out.flag)
    split_release_active_for_assembly(run_dump_active_rs_for_assembly.release_active_rs)
    run_release_active_rs_for_assembly(split_release_active_for_assembly.out.release_active_chunks)

}

process initiate_release_status_for_assembly {

    label 'short_time', 'med_mem'


    output:
    val true, emit: flag

    script:
    """
    export PYTHONPATH=$params.python_path
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
    pipeline_parameters += " --rsAccDumpFile=" + "release_active_dump"
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

process run_release_active_rs_for_assembly {

    label 'long_time', 'med_mem'

    input:
    each path(rslist)

    output:
    path "release_active_dump", emit: release_active_dump

    script:
    def pipeline_parameters = " --spring.batch.job.names=ACTIVE_ACCESSIONS_RELEASE_FROM_DB_JOB"
    pipeline_parameters += " --rsAccFile=" + rslist
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.release_pipeline --spring.config.location=file:$params.release_job_props $pipeline_parameters
    """
}



