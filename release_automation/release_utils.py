# Copyright 2020 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import datetime
import os
from functools import lru_cache

from ebi_eva_common_pyutils.taxonomy import taxonomy
from ebi_eva_internal_pyutils.pg_utils import get_all_results_for_query

collections_assembly_attribute_map = {
    "dbsnpSubmittedVariantEntity": "seq",
    "dbsnpSubmittedVariantOperationEntity": "inactiveObjects.seq",
    "submittedVariantEntity": "seq",
    "submittedVariantOperationEntity": "inactiveObjects.seq",
    "dbsnpClusteredVariantEntity": "asm",
    "dbsnpClusteredVariantOperationEntity": "inactiveObjects.asm",
    "clusteredVariantEntity": "asm",
    "clusteredVariantOperationEntity": "inactiveObjects.asm"
}

submitted_collections_taxonomy_attribute_map = {
    "dbsnpSubmittedVariantEntity": "tax",
    "dbsnpSubmittedVariantOperationEntity": "inactiveObjects.tax",
    "submittedVariantEntity": "tax",
    "submittedVariantOperationEntity": "inactiveObjects.tax"
}

def update_release_progress_status(metadata_connection_handle, release_species_inventory_table, taxonomy,
                                   assembly_accession, release_version, release_status):
    if release_status == 'Started':
        date_to_change = 'release_start'
    else:
        date_to_change = 'release_end'
    now = datetime.datetime.now().isoformat()
    update_status_query = (
        f"update {release_species_inventory_table} "
        f"set release_status = '{release_status}', {date_to_change} = '{now}' "
        f"where taxonomy = {taxonomy} "
        f"and assembly_accession = '{assembly_accession}' "
        f"and release_version = {release_version};"
    )
    with metadata_connection_handle.cursor() as cursor:
        cursor.execute(update_status_query)
    metadata_connection_handle.commit()


def get_release_assemblies_for_taxonomy(taxonomy_id, release_species_inventory_table,
                                        release_version, metadata_connection_handle):
    results = get_all_results_for_query(metadata_connection_handle, "select assembly_accession from {0} "
                                                                    "where taxonomy = '{1}' "
                                                                    "and release_version = {2} and should_be_released "
                                                                    "and num_rs_to_release > 0"
                                        .format(release_species_inventory_table, taxonomy_id, release_version))
    if len(results) == 0:
        raise Exception("Could not find assemblies pertaining to taxonomy ID: " + taxonomy_id)
    return [result[0] for result in results]


def get_release_inventory_info_for_assembly(taxonomy_id, assembly_accession, release_species_inventory_table,
                                            release_version, metadata_connection_handle):
    results = get_all_results_for_query(metadata_connection_handle, "select row_to_json(row) from "
                                                                    "(select * from {0} where "
                                                                    "taxonomy = '{1}' "
                                                                    "and assembly_accession = '{2}' "
                                                                    "and release_version = {3} "
                                                                    "and should_be_released "
                                                                    "and num_rs_to_release > 0) row"
                                        .format(release_species_inventory_table, taxonomy_id, assembly_accession,
                                                release_version))
    if len(results) == 0:
        raise Exception("Could not find release inventory pertaining to taxonomy ID: {0} and assembly: {1} "
                        .format(taxonomy_id, assembly_accession))
    return results[0][0]


def get_release_for_status_and_version(release_species_inventory_table, metadata_connection_handle, status=None,
                                       release_version=None, taxonomy_id=None, assembly_accessions=None):
    def format_list(list_to_format):
        return f"({str(list_to_format).strip('[]')})"

    if status:
        if 'Pending' in status:
            status.remove('Pending')
            if status:
                status_statement = f"and (release_status in {format_list(status)} or release_status is null) "
            else:
                status_statement = f"and release_status is null "
        else:
            status_statement = f"and release_status in {format_list(status)} "
    else:
        status_statement = None

    query = ''.join((
        f"select taxonomy, assembly_accession, release_version, release_status from {release_species_inventory_table} ",
        f"where should_be_released ",
        f"and num_rs_to_release > 0 ",
        f"{status_statement}" if status_statement else '',
        f"and release_version={release_version} " if release_version else '',
        f"and taxonomy={taxonomy_id} " if taxonomy_id else '',
        f"and assembly_accessions in {format_list(assembly_accessions)} " if assembly_accessions else '',
        "ORDER BY release_version, taxonomy, assembly_accession"
    ))
    return get_all_results_for_query(metadata_connection_handle, query)


@lru_cache
def get_release_folder_name(taxonomy_id):
    return taxonomy.get_normalized_scientific_name_from_ensembl(taxonomy_id)
