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
import sys
import traceback
from argparse import ArgumentParser

import click
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_internal_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_internal_pyutils.spring_properties import SpringPropertiesGenerator

from release_automation.release_config import load_config
from release_automation.release_utils import get_release_inventory_info_for_assembly

logger = logging_config.get_logger(__name__)


def create_release_properties_file_for_assembly(taxonomy_id, assembly_accession, release_version,
                                                assembly_release_folder, job_name='ACCESSION_RELEASE_JOB',
                                                file_name='release'):
    os.makedirs(assembly_release_folder, exist_ok=True)
    output_file = os.path.join(assembly_release_folder, f"{assembly_accession}_{file_name}.properties")
    with get_metadata_connection_handle(cfg['maven']['environment'], cfg['maven']['settings_file']) as metadata_connection_handle:
        release_species_inventory_table = cfg['release']['release_species_inventory_table'],

        release_inventory = get_release_inventory_info_for_assembly(
            taxonomy_id, assembly_accession, release_species_inventory_table, release_version,
            metadata_connection_handle
        )

    properties_string = SpringPropertiesGenerator(cfg['maven']['environment'], cfg['maven']['settings_file']).get_release_properties(
        temp_mongo_db='', # No need for temp mongo
        assembly_accession=assembly_accession,
        taxonomy_accession=taxonomy_id,
        fasta=release_inventory['fasta_path'],
        assembly_report=release_inventory['report_path'],
        contig_naming='SEQUENCE_NAME',
        output_folder=assembly_release_folder
    )
    with open(output_file, "w") as open_file:
        open_file.write(properties_string)

def run_release_for_assembly(taxonomy_id, assembly_accession, release_version):
    exit_code = -1
    try:
        release_properties_file = create_release_properties_file_for_assembly(private_config_xml_file, profile,
                                                                              taxonomy_id, assembly_accession,
                                                                              release_species_inventory_table,
                                                                              release_version, assembly_release_folder)
        release_command = 'java -Xmx{0}g -jar {1} --spring.config.location=file:{2} --spring.data.mongodb.port={3}' \
            .format(memory, release_jar_path, release_properties_file, mongo_port)
        run_command_with_output("Running release pipeline for assembly: " + assembly_accession, release_command)
        exit_code = 0
    except Exception as ex:
        logger.error("Encountered an error while running release for assembly: " + assembly_accession + "\n"
                     + traceback.format_exc())
        exit_code = -1
    finally:
        close_mongo_port_to_tempmongo(port_forwarding_process_id)
        logger.info("Java release pipeline run completed with exit_code: " + str(exit_code))
        sys.exit(exit_code)


def main():
    parser = ArgumentParser()
    parser.add_argument("--taxonomy-id", help="ex: 9913", required=True)
    parser.add_argument("--assembly-accession", help="ex: GCA_000003055.6", required=True)
    parser.add_argument("--release-version", help="ex: 2", type=int, required=True)
    parser.add_argument("--assembly-release-folder", required=True)
    args = parser.parse_args()
    logging_config.add_stdout_handler()

    load_config()

    run_release_for_assembly(args.taxonomy_id, args.assembly_accession, args.release_version)


if __name__ == "__main__":
    main()
