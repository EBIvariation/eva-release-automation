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
from argparse import ArgumentParser
from functools import lru_cache
from random import choice
from string import ascii_lowercase

import yaml
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_common_pyutils.common_utils import pretty_print
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_internal_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_internal_pyutils.spring_properties import SpringPropertiesGenerator

from release_automation.release_config import load_config
from release_automation.release_utils import get_release_inventory_info_for_assembly, \
    get_release_assemblies_for_taxonomy, get_release_for_status_and_version, get_release_folder_name

logger = logging_config.get_logger(__name__)


def get_nextflow_params(taxonomy_id, assembly_accession, release_version):
    release_dir = get_assembly_release_folder(release_version, taxonomy_id, assembly_accession)
    config_param = os.path.join(release_dir, f'nextflow_params_{taxonomy_id}_{assembly_accession}.yaml')
    # Add the same python interpreter as the one we're using to use with the python step scripts
    cfg['executable']['python_interpreter'] = sys.executable
    release_job_props = get_release_java_props(release_dir, assembly_accession, taxonomy_id, release_version, 'release')
    yaml_data = {
        'assembly': assembly_accession,
        'executable': cfg['executable'],
        'jar': cfg['jar'],
        'log_file': get_release_log_file_name(release_version, taxonomy_id, assembly_accession),
        'maven': cfg['maven'],
        'release_version': release_version,
        'assembly_folder': release_dir,
        'taxonomy': taxonomy_id,
        'release_job_props': release_job_props,
        'output_dir': release_dir
    }
    with open(config_param, 'w') as open_file:
        yaml.safe_dump(yaml_data, open_file)
    return config_param

def get_release_java_props(assembly_release_folder, assembly_accession, taxonomy_id, release_version, file_name):
    os.makedirs(assembly_release_folder, exist_ok=True)
    output_file = f"{assembly_release_folder}/{assembly_accession}_{file_name}.properties"
    private_config_xml_file = cfg.query("maven", "settings_file")
    profile = cfg.query("maven", "environment")
    release_species_inventory_table = release_species_inventory_table = cfg.query('release', 'inventory_table')
    with get_metadata_connection_handle(profile, private_config_xml_file) as metadata_connection_handle:
        release_properties = get_release_inventory_info_for_assembly(taxonomy_id, assembly_accession, release_species_inventory_table,
                                                release_version, metadata_connection_handle)
    properties_string = SpringPropertiesGenerator(profile, private_config_xml_file).get_release_properties(
        temp_mongo_db=None,
        job_name=None,
        assembly_accession=assembly_accession,
        taxonomy_accession=taxonomy_id,
        fasta=release_properties['fasta_path'],
        assembly_report=release_properties['report_path'],
        contig_naming='INSDC',
        output_folder=None
    )
    open(output_file, "w").write(properties_string)
    return output_file


def get_nextflow_config():
    if 'RELEASE_NEXTFLOW_CONFIG' in os.environ and os.path.isfile(os.environ['RELEASE_NEXTFLOW_CONFIG']):
        return os.path.abspath(os.environ['RELEASE_NEXTFLOW_CONFIG'])


def get_run_release_for_assembly_nextflow():
    curr_dir = os.path.dirname(__file__)
    return os.path.join(curr_dir, 'nextflow', 'run_release_for_assembly.nf')


def get_release_log_file_name(release_version, taxonomy_id, assembly_accession):
    return f"{get_assembly_release_folder(release_version, taxonomy_id, assembly_accession)}/release_{taxonomy_id}_{assembly_accession}.log"


@lru_cache
def get_release_folder(release_version):
    folder = os.path.join(cfg.query('release', 'release_output'), f'release_{release_version}')
    os.makedirs(folder, exist_ok=True)
    return folder


@lru_cache
def get_species_release_folder(release_version, taxonomy_id):
    folder = os.path.join(get_release_folder(release_version), get_release_folder_name(taxonomy_id))
    os.makedirs(folder, exist_ok=True)
    return folder


@lru_cache
def get_assembly_release_folder(release_version, taxonomy_id, assembly_accession):
    folder = os.path.join(get_species_release_folder(release_version, taxonomy_id), assembly_accession)
    os.makedirs(folder, exist_ok=True)
    return folder


def run_release_for_species(taxonomy_id, release_assemblies, release_version, resume=False):
    private_config_xml_file = cfg.query("maven", "settings_file")
    profile = cfg.query("maven", "environment")
    release_species_inventory_table = cfg.query('release', 'inventory_table')
    with get_metadata_connection_handle(profile, private_config_xml_file) as metadata_connection_handle:
        if not release_assemblies:
            release_assemblies = get_release_assemblies_for_taxonomy(
                taxonomy_id, release_species_inventory_table, release_version, metadata_connection_handle
            )

        for assembly_accession in release_assemblies:
            nextflow_params = get_nextflow_params(taxonomy_id, assembly_accession, release_version)
            workflow_file_path = get_run_release_for_assembly_nextflow()
            release_dir = get_assembly_release_folder(release_version, taxonomy_id, assembly_accession)
            nextflow_config = get_nextflow_config()
            random_string = ''.join(choice(ascii_lowercase) for i in range(4))
            run_name = f'release_{release_version}_{taxonomy_id}_{assembly_accession}_{random_string}'.replace('.', '_')
            workflow_command = ' '.join((
                f"cd {release_dir} &&",
                f"{cfg.query('executable', 'nextflow')} run {workflow_file_path}",
                f"-name {run_name}",
                f"-params-file {nextflow_params}",
                f'-c {nextflow_config}' if nextflow_config else '',
                '-resume' if resume else '',
            ))
            logger.info(f"Running workflow file {workflow_file_path} with the following command: "
                        f"\n {workflow_command} \n")
            run_command_with_output(f"Running workflow file {workflow_file_path}", workflow_command)


def list_release_per_status(status, release_version, taxonomy_id, assembly_accessions):
    private_config_xml_file = cfg.query("maven", "settings_file")
    profile = cfg.query("maven", "environment")
    release_species_inventory_table = cfg.query('release', 'inventory_table')
    with get_metadata_connection_handle(profile, private_config_xml_file) as metadata_connection_handle:
        header = ['taxonomy', 'assembly_accession', 'release_version', 'release_status']
        table = []
        for taxonomy, assembly_accession, release_version, release_status in get_release_for_status_and_version(
                release_species_inventory_table, metadata_connection_handle, status, release_version=release_version,
                taxonomy_id=taxonomy_id, assembly_accessions=assembly_accessions
        ):
            if release_status is None:
                release_status = 'Pending'
            table.append((taxonomy, assembly_accession, release_version, release_status))
        pretty_print(header, table)


def main():
    argparse = ArgumentParser()
    argparse.add_argument("--list_status", nargs='+',
                          help="Generate the list of species and assembly that needs to be released",
                          choices=['Pending', 'Started', 'Completed'])
    argparse.add_argument("--taxonomy_id", help="ex: 9913")
    argparse.add_argument("--assembly_accessions", nargs='+', help="ex: GCA_000003055.3")
    argparse.add_argument("--release_version")
    argparse.add_argument("--resume", default=False, action='store_true',
                          help="Resume the nextflow pipeline for the specified taxonomy and assembly")
    args = argparse.parse_args()
    load_config()
    logging_config.add_stdout_handler()
    if args.list_status:
        list_release_per_status(args.list_status, args.release_version, args.taxonomy_id, args.assembly_accessions)
    else:
        if not args.taxonomy_id:
            logger.error('--taxonomy_id is required when running the release')
            sys.exit(1)
        if not args.release_version:
            logger.error('--release_version is required when running the release')
            sys.exit(1)
        run_release_for_species(args.taxonomy_id, args.assembly_accessions, args.release_version, args.resume)


if __name__ == "__main__":
    main()
