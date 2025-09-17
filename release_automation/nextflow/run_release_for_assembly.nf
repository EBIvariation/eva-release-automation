#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.active_chunk_size = 10000000
params.merged_and_deprecated_chunk_size = 100000
output_log = params.output_dir + '/log_files'

workflow {

    initiate_release_status_for_assembly()

    // Release active variants
    run_dump_active_rs_for_assembly(initiate_release_status_for_assembly.out.flag)
    split_release_active_for_assembly(run_dump_active_rs_for_assembly.out.release_active_rs, params.active_chunk_size)
    release_active_rs_for_assembly(split_release_active_for_assembly.out.release_active_chunks)
    sort_and_index_chunk_active(release_active_rs_for_assembly.out.release_active_chunk)
    merge_active_chunks(sort_and_index_chunk_active.out.sorted_release_active_chunk.collect(), sort_and_index_chunk_active.out.index_release_active_chunk.collect())

    // Release merged and deprecated variants
    run_dump_merged_and_deprecated_rs_for_assembly(initiate_release_status_for_assembly.out.flag)
    split_release_merged_and_deprecated_for_assembly(run_dump_merged_and_deprecated_rs_for_assembly.out.release_merged_and_deprecated_rs, params.merged_and_deprecated_chunk_size)
    release_merged_and_deprecated_rs_for_assembly(split_release_merged_and_deprecated_for_assembly.out.release_merged_and_deprecated_chunks)
    sort_and_index_chunk_merged(release_merged_and_deprecated_rs_for_assembly.out.release_merged_chunk)
    merge_merged_chunks(sort_and_index_chunk_merged.out.sorted_release_merged_chunk.collect(), sort_and_index_chunk_merged.out.index_release_merged_chunk.collect())
    merge_deprecated_chunks(release_merged_and_deprecated_rs_for_assembly.out.release_deprecated_chunk.collect())

    vcf_channel = channel.of().concat(merge_active_chunks.out.release_active_merged, merge_merged_chunks.out.release_merged_merged)
    // copy the files with the updated sequence to the output directory
    update_sequence_names_to_ena(vcf_channel)

    // Validate VCF files
    vcf_validator_release_vcf_files(vcf_channel)
    assembly_check_release_vcf_files(vcf_channel)
    analyze_vcf_validator_results(vcf_validator_release_vcf_files.out.vcf_validator_results)
    analyze_assembly_checker_results(assembly_check_release_vcf_files.out.assembly_check_report)

    // count
    count_rs_ids_in_release_vcf(update_sequence_names_to_ena.out.release_vcf_output_file)
    count_rs_ids_in_release_txt(merge_deprecated_chunks.out.release_deprecated)
    merge_count_files(count_rs_ids_in_release_vcf.out.count_vcf.collect(), count_rs_ids_in_release_txt.out.count_txt)

    update_release_status_for_assembly(merge_count_files.out.readme_count)
}


process initiate_release_status_for_assembly {

    label 'short_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    output:
    val true, emit: flag

    script:
    log_file = "initiate_release_${params.taxonomy}_${params.assembly}_${task.index}.log"
    """
    $params.executable.python_interpreter -m release_automation.initiate_release_status_for_assembly --taxonomy_id $params.taxonomy --assembly_accession $params.assembly --release_version $params.release_version 1>> $log_file 2>&1
    """
}

process run_dump_active_rs_for_assembly {

    label 'long_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    val flag

    output:
    path "release_active_dump", emit: release_active_rs
    path "dump_active_rs_${params.taxonomy}_${params.assembly}_${task.index}.log", emit: log_file

    script:
    log_file = "dump_active_rs_${params.taxonomy}_${params.assembly}_${task.index}.log"
    def pipeline_parameters = " --spring.batch.job.names=DUMP_ACTIVE_ACCESSIONS_JOB"
    pipeline_parameters += " --parameters.outputFolder=\$PWD"
    pipeline_parameters += " --parameters.rsAccDumpFile=" + "release_active_dump"
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.release_pipeline --spring.config.location=file:$params.release_job_props $pipeline_parameters 1>> $log_file 2>&1
    """
}

process split_release_active_for_assembly {

    label 'long_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    path "release_active_dump"
    val chunk_size

    output:
    path "active-rs-chunk-*", emit: release_active_chunks
    path "split_active_rs_${params.taxonomy}_${params.assembly}_${task.index}.log", emit: log_file

    script:
    log_file = "split_active_rs_${params.taxonomy}_${params.assembly}_${task.index}.log"
    """
    if [ -s release_active_dump ]
    then
      split -a 5 -d -l ${chunk_size} release_active_dump active-rs-chunk- 1>> $log_file 2>&1
    else
      touch active-rs-chunk-0 > $log_file
    fi
    """
}

process release_active_rs_for_assembly {
    maxForks 50

    label 'med_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    each path(rs_chunk)

    output:
    path "current_ids_${rs_chunk}.vcf", emit: release_active_chunk
    path "release_active_rs_${params.taxonomy}_${params.assembly}_${task.index}.log", emit: log_file

    script:
    log_file = "release_active_rs_${params.taxonomy}_${params.assembly}_${task.index}.log"

    def pipeline_parameters = " --spring.batch.job.names=ACTIVE_ACCESSIONS_RELEASE_FROM_DB_JOB"
    pipeline_parameters += " --parameters.rsAccFile=" + rs_chunk
    pipeline_parameters += " --parameters.outputFolder=\$PWD"
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.release_pipeline --spring.config.location=file:$params.release_job_props $pipeline_parameters 1>> $log_file 2>&1
    mv *_current_ids.vcf current_ids_${rs_chunk}.vcf
    """
}


process sort_and_index_chunk_active {

    label 'med_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    path(release_active_chunk)

    output:
    path "*_sorted.vcf.gz", emit: sorted_release_active_chunk
    path "*_sorted.vcf.gz.csi", emit: index_release_active_chunk
    path "sort_and_index_active_rs_${params.taxonomy}_${params.assembly}_${task.index}.log", emit: log_file

    script:
    log_file = "sort_and_index_active_rs_${params.taxonomy}_${params.assembly}_${task.index}.log"
    """
    awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' $release_active_chunk | $params.executable.bcftools view --no-version  -O z -o ${release_active_chunk}_sorted.vcf.gz -  1>> $log_file 2>&1
    $params.executable.bcftools index ${release_active_chunk}_sorted.vcf.gz 1>> $log_file 2>&1
    """
}


process merge_active_chunks {

    label 'med_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    path(release_active_chunks)
    path(index_release_active_chunks)

    output:
    path "$active_vcf", emit: release_active_merged
    path "merge_active_rs_${params.taxonomy}_${params.assembly}_${task.index}.log", emit: log_file

    script:
    active_vcf = "${params.taxonomy}_${params.assembly}_current_ids.before_rename.vcf.gz"
    log_file = "merge_active_rs_${params.taxonomy}_${params.assembly}_${task.index}.log"
    """
    $params.executable.bcftools concat --no-version -a -o $active_vcf -O z $release_active_chunks 1>> $log_file 2>&1
    """
}

process run_dump_merged_and_deprecated_rs_for_assembly {

    label 'long_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    val flag

    output:
    path "release_merged_and_deprecated_dump", emit: release_merged_and_deprecated_rs
    path "dump_merged_and_deprecated_rs_${params.taxonomy}_${params.assembly}_${task.index}.log", emit: log_file

    script:
    log_file = "dump_merged_and_deprecated_rs_${params.taxonomy}_${params.assembly}_${task.index}.log"
    def pipeline_parameters = " --spring.batch.job.names=DUMP_MERGED_AND_DEPRECATED_ACCESSIONS_JOB"
    pipeline_parameters += " --parameters.outputFolder=\$PWD"
    pipeline_parameters += " --parameters.rsAccDumpFile=" + "release_merged_and_deprecated_dump"
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.release_pipeline --spring.config.location=file:$params.release_job_props $pipeline_parameters 1>> $log_file 2>&1
    """
}

process split_release_merged_and_deprecated_for_assembly {

    label 'long_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    path "release_merged_and_deprecated_dump"
    val chunk_size

    output:
    path "merged-and-deprecated-rs-chunk-*", emit: release_merged_and_deprecated_chunks
    path "split_merged_and_deprecated_rs_${params.taxonomy}_${params.assembly}_${task.index}.log", emit: log_file

    script:
    log_file = "split_merged_and_deprecated_rs_${params.taxonomy}_${params.assembly}_${task.index}.log"
    """
    set -o pipefail

    if [ -s release_merged_and_deprecated_dump ]
    then
      sort -u release_merged_and_deprecated_dump | split -a 5 -d -l ${chunk_size} - merged-and-deprecated-rs-chunk- 1>> $log_file 2>&1
    else
      touch merged-and-deprecated-rs-chunk-0 > $log_file
    fi
    """
}

process release_merged_and_deprecated_rs_for_assembly {
    maxForks 50

    label 'med_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    each path(rs_chunk)

    output:
    path "merged_ids_${rs_chunk}.vcf", emit: release_merged_chunk
    path "deprecated_ids_${rs_chunk}.txt",  emit: release_deprecated_chunk
    path "release_merged_and_deprecated_rs_${params.taxonomy}_${params.assembly}_${task.index}.log", emit: log_file

    script:
    log_file = "release_merged_and_deprecated_rs_${params.taxonomy}_${params.assembly}_${task.index}.log"
    def pipeline_parameters = " --spring.batch.job.names=MERGED_AND_DEPRECATED_ACCESSIONS_RELEASE_FROM_DB_JOB"
    pipeline_parameters += " --parameters.rsAccFile=" + rs_chunk
    pipeline_parameters += " --parameters.outputFolder=\$PWD"
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.release_pipeline --spring.config.location=file:$params.release_job_props $pipeline_parameters 1>> $log_file 2>&1
    mv *_merged_ids.vcf merged_ids_${rs_chunk}.vcf
    mv *_deprecated_ids.unsorted.txt deprecated_ids_${rs_chunk}.txt
    """
}


process sort_and_index_chunk_merged {

    label 'med_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    path(release_merged_chunk)

    output:
    path "*_sorted.vcf.gz", emit: sorted_release_merged_chunk
    path "*_sorted.vcf.gz.csi", emit: index_release_merged_chunk
    path "sort_index_merged_rs_${params.taxonomy}_${params.assembly}_${task.index}.log", emit: log_file

    script:
    log_file = "sort_index_merged_rs_${params.taxonomy}_${params.assembly}_${task.index}.log"
    """
    awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' $release_merged_chunk | $params.executable.bcftools view --no-version  -O z -o ${release_merged_chunk}_sorted.vcf.gz - 1>> $log_file 2>&1
    $params.executable.bcftools index ${release_merged_chunk}_sorted.vcf.gz 1>> $log_file 2>&1
    """
}


process merge_merged_chunks {

    label 'med_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    path(release_merged_chunks)
    path(index_release_merged_chunk)

    output:
    path "$merged_vcf", emit: release_merged_merged
    path "merge_merged_rs_${params.taxonomy}_${params.assembly}_${task.index}.log", emit: log_file

    script:
    log_file = "merge_merged_rs_${params.taxonomy}_${params.assembly}_${task.index}.log"
    merged_vcf = "${params.taxonomy}_${params.assembly}_merged_ids.before_rename.vcf.gz"

    """
    $params.executable.bcftools concat --no-version -a -o $merged_vcf -O z $release_merged_chunks 1>> $log_file 2>&1
    """
}


process merge_deprecated_chunks {

    label 'med_time', 'med_mem'

    publishDir path: params.output_dir, pattern: '*.txt.gz', mode: 'copy', overwrite: true

    input:
    path(release_deprecated_chunks)


    output:
    path "$final_file_deprecated", emit: release_deprecated

    script:
    final_file_deprecated = "${params.taxonomy}_${params.assembly}_deprecated_ids.txt.gz"
    """
    cat $release_deprecated_chunks | sort -u | gzip -c  > $final_file_deprecated
    """
}



process vcf_validator_release_vcf_files {

    label 'med_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    path vcf_file

    output:
    path "vcf_format/*.errors.*.txt", emit: vcf_validator_results
    path "vcf_format/*.vcf_format.log", emit: vcf_validator_log

    script:
    """
    trap 'if [[ \$? == 1 ]]; then exit 0; fi' EXIT

    mkdir -p vcf_format
    $params.executable.vcf_validator -i ${vcf_file}  -r text -o vcf_format > vcf_format/${vcf_file}.vcf_format.log 2>&1
    """
}

process assembly_check_release_vcf_files {

    label 'med_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    path(vcf_file)

    output:
    path "assembly_check/*text_assembly_report*", emit: assembly_check_report
    path "assembly_check/*.assembly_check.log", emit: assembly_check_log

    script:
    """
    trap 'if [[ \$? == 1 || \$? == 139 ]]; then exit 0; fi' EXIT

    mkdir -p assembly_check
    $params.executable.vcf_assembly_checker -i $vcf_file -f $params.fasta_file -a $params.assembly_report -r summary,text  -o assembly_check > assembly_check/${vcf_file}.assembly_check.log 2>&1
    """
}

process analyze_vcf_validator_results {

    label 'med_time', 'default_mem'

    input:
    path vcf_validator_log

    output:
    val true, emit: flag

    script:
    """
    echo "Duplicated variant" > allowed_errors.txt
    echo "Warning: Reference and alternate alleles " >> allowed_errors.txt
    echo "do not share the first nucleotide" >> allowed_errors.txt
    echo "the input file is not valid" >> allowed_errors.txt
    echo "the input file is valid" >> allowed_errors.txt
    echo "not listed in a valid meta-data ALT entry" >> allowed_errors.txt
    echo "Chromosome is not a string without colons or whitespaces" >> allowed_errors.txt

    NB_ERROR=\$(cat $vcf_validator_log | grep -vFf allowed_errors.txt | wc -l)
    if [ ! "\$NB_ERROR" -eq "0" ]; then
        echo "\$NB_ERROR Unusual error(s) found in VCF validation log: $vcf_validator_log"
        exit 1
    fi
    """
}


process analyze_assembly_checker_results {

    label 'med_time', 'med_mem'

    input:
    path assembly_checker_log

    output:
    val true, emit: flag

    script:
    """
    echo "not present in FASTA file"> allowed_errors.txt
    echo "does not match the reference sequence" >> allowed_errors.txt
    echo "Multiple synonyms  found for contig" >> allowed_errors.txt

    NB_ERROR=\$(cat $assembly_checker_log | grep -vFf allowed_errors.txt | wc -l)
    if [ ! "\$NB_ERROR" -eq "0" ]; then
        echo "\$NB_ERROR Unusual error(s) found in assembly report log: $assembly_checker_log"
        exit 1
    fi
    """
}

process count_rs_ids_in_release_vcf {

    label 'med_time', 'med_mem'

    input:
    path vcf_file

    output:
    path "${vcf_file}.count", emit: count_vcf

    script:
    """
    COUNT=\$(gunzip -c $vcf_file | grep -v "^#" | cut -f 3 | sed s/";"/\\n/g | sort -T . -u --parallel=4 | wc -l)
    echo -e "$vcf_file\\t\${COUNT}" > ${vcf_file}.count
    """
}

process count_rs_ids_in_release_txt {

    label 'med_time', 'med_mem'

    input:
    path txt_file

    output:
    path "${txt_file}.count", emit: count_txt

    script:
    """
    COUNT=\$(gunzip -c $txt_file | cut -f1 | uniq | wc -l)
    echo -e "$txt_file\\t\${COUNT}" > ${txt_file}.count
    """
}


process merge_count_files {

    label 'long_time', 'med_mem'

    publishDir path: params.output_dir, mode: 'copy', overwrite: true

    input:
    path vcf_based_ids_counts
    path txt_based_ids_count

    output:
    path "README_rs_ids_counts.txt", emit: readme_count

    script:
    log_file = "validate_rs_release_files_${params.taxonomy}_${params.assembly}_${task.index}.log"
    """
    cat  $vcf_based_ids_counts $txt_based_ids_count | sort > README_rs_ids_counts.txt
    """


}


process update_sequence_names_to_ena {

    label 'med_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true
    publishDir path: params.output_dir, pattern: '*.vcf.gz', mode: 'copy', overwrite: true


    input:
    path release_vcf_file

    output:
    path release_vcf_output_file, emit: release_vcf_output_file
    path release_vcf_output_index, emit: release_vcf_output_index

    script:
    log_file = "update_sequence_names_to_ena_${release_vcf_file.getSimpleName()}.log"
    // Remove three extensions ".before_rename.vcf.gz". Cannot use getSimpleName as there is a "." in the assembly accession
    release_vcf_output_file=release_vcf_file.getName().replace('before_rename.','')
    release_vcf_output_index="${release_vcf_output_file}.csi"
    """
    $params.executable.convert_vcf_file -i $release_vcf_file -o $release_vcf_output_file -c enaSequenceName 1>> $log_file 2>&1
    $params.executable.bcftools index $release_vcf_output_file 1>> $log_file 2>&1
    """
}

process update_release_status_for_assembly {

    label 'short_time', 'med_mem'

    publishDir path: output_log, pattern: '*.log', mode: 'copy', overwrite: true

    input:
    path readme_count

    output:
    val true, emit: flag

    script:
    log_file = "finalise_release_${params.taxonomy}_${params.assembly}_${task.index}.log"
    """
    $params.executable.python_interpreter -m release_automation.update_release_status_for_assembly --taxonomy_id $params.taxonomy --assembly_accession $params.assembly --release_version $params.release_version 1>> $log_file 2>&1

    """
}