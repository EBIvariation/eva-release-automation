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

import os
import signal
import traceback
from functools import lru_cache

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.network_utils import get_available_local_port, forward_remote_port_to_local_port
from ebi_eva_common_pyutils.taxonomy import taxonomy
from ebi_eva_internal_pyutils.metadata_utils import get_metadata_connection_handle

from run_release_in_embassy.release_metadata import get_target_mongo_instance_for_assembly

logger = logging_config.get_logger(__name__)


def get_bgzip_bcftools_index_commands_for_file(bgzip_path, bcftools_path, file):
    commands = ["rm -f {0}.gz".format(file), "({0} -cf {1} > {1}.gz)".format(bgzip_path, file),
                "({0} index --csi {1}.gz)".format(bcftools_path, file)]
    return commands


def get_release_vcf_file_name(assembly_release_folder, taxonomy_id, assembly_accession, vcf_file_category):
    return os.path.join(assembly_release_folder, "{0}_{1}_{2}.vcf".format(taxonomy_id, assembly_accession,
                                                                          vcf_file_category))


def get_release_vcf_file_name_genbank(assembly_release_folder, taxonomy_id, assembly_accession, vcf_file_category):
    return os.path.join(
        assembly_release_folder,
        "{0}_{1}_{2}_with_genbank.vcf".format(taxonomy_id, assembly_accession, vcf_file_category)
    )


def get_unsorted_release_vcf_file_name(assembly_release_folder, taxonomy_id, assembly_accession, vcf_file_category):
    vcf_file_path = get_release_vcf_file_name(assembly_release_folder, taxonomy_id, assembly_accession,
                                              vcf_file_category)
    filename = os.path.basename(vcf_file_path)
    return vcf_file_path.replace(filename, filename.replace(".vcf", "_unsorted.vcf"))


def get_release_text_file_name(assembly_release_folder, taxonomy_id, assembly_accession, release_text_file_category):
    return os.path.join(assembly_release_folder,
                        "{0}_{1}_{2}.txt".format(taxonomy_id, assembly_accession, release_text_file_category))


def get_unsorted_release_text_file_name(assembly_release_folder, taxonomy_id, assembly_accession,
                                        release_text_file_category):
    release_text_file_path = get_release_text_file_name(assembly_release_folder, taxonomy_id, assembly_accession,
                                                        release_text_file_category)
    filename = os.path.basename(release_text_file_path)
    return release_text_file_path.replace(filename, filename.replace(".txt", ".unsorted.txt"))


def get_release_db_name_in_tempmongo_instance(taxonomy_id, assembly_accession):
    return "acc_" + str(taxonomy_id) + "_" + assembly_accession.replace('.', '_')


@lru_cache
def get_release_folder_name(taxonomy_id):
    return taxonomy.get_normalized_scientific_name_from_ensembl(taxonomy_id)
