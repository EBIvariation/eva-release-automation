import os
import time
from unittest import TestCase

import yaml
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_internal_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_internal_pyutils.pg_utils import execute_query

from release_automation.run_release_for_species import run_release_for_species

resource_folder = os.path.join(os.path.dirname(__file__), 'resources')
release_inventory = [
    {'taxonomy': 4530, 'scientific_name':'Oryza sativa', 'assembly_accession': 'GCA_000005425.2', 'release_version': 1,
     'sources':'EVA', 'fasta_path': 'path_to_fasta', 'report_path': 'path_to_report', 'tempmongo_instance':1,
     'should_be_released':True, 'num_rs_to_release':1, 'total_num_variants': 1, 'release_folder_name': 'oryza_sativa',
     'release_status': None, 'release_start': None, 'release_end': None}
]

class TestRunReleaseForSpecies(TestCase):

    compose_dir = os.path.join(os.path.dirname(__file__), 'docker')

    @classmethod
    def setUpClass(cls):
        os.chdir(cls.compose_dir)
        command = 'docker compose up -d --wait'
        run_command_with_output('Start docker-compose in the background', command)
        # Wait for the postgres to be ready to accept command
        # time.sleep(10)



    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.compose_dir)
        command = 'docker compose down'
        run_command_with_output('Stop docker-compose running in the background', command)


    def test_run_release_for_species(self):
        config_xml = os.path.join(resource_folder, 'config_xml_file.xml')
        config_yaml = os.path.join(resource_folder, 'release_config.yml')
        nextflow_config = os.path.join(resource_folder, 'nextflow.config')
        output_dir = os.path.join(resource_folder, 'release_output_dir')
        nextflow_path = 'nextflow'
        bcftools_path = 'bcftools'
        profile = 'localhost'
        inventory_table =  'eva_progress_tracker.clustering_release_tracker'
        # cfg.content = {
        #     'maven': {'environment': profile, 'settings_file': config_xml},
        #     'release':{'release_output': output_dir, 'inventory_table': inventory_table},
        #     'executable': {'nextflow': nextflow_path, 'bcftools': bcftools_path},
        #     'jar':{'release_pipeline': '/Users/tcezard/PycharmProjects/eva-release-automation/tests/release_automation/resources/release_output_dir/release_1/oryza_sativa/GCA_000005425.2/eva-accession-release-0.6.47-20250422.150441-1-exec.jar'}
        # }
        fasta_file = os.path.join(resource_folder,'GCA_000005425.2.fa')
        report_file = os.path.join(resource_folder,'GCA_000005425.2_assembly_report.txt')

        with get_metadata_connection_handle(profile, config_xml) as metadata_connection_handle:
            insert_query = (f"insert into {inventory_table} "
                            f"(taxonomy, scientific_name, assembly_accession, release_version, sources, fasta_path, "
                            f"report_path, should_be_released, num_rs_to_release, total_num_variants, release_folder_name) "
                            f"values (4530, 'Oryza sativa', 'GCA_000005425.2', 1, 'EVA', '{fasta_file}', "
                            f"'{report_file}', True, 1, 1, 'oryza_sativa')")
            execute_query(metadata_connection_handle, insert_query)
        # os.environ['RELEASE_CONFIG'] = config_yaml
        # os.environ['RELEASE_NEXTFLOW_CONFIG'] = nextflow_config
        command = 'docker exec executor python3 -m  release_automation.run_release_for_species --list_status Pending'
        output = run_command_with_output('list all pending', command, True)
        expected_output='''| taxonomy | assembly_accession | release_version | release_status |
|     4530 |    GCA_000005425.2 |               1 |        Pending |
'''
        assert expected_output == output
        command = 'docker exec executor python3 -m  release_automation.run_release_for_species --taxonomy_id 4530 --assembly_accessions GCA_000005425.2 --release_version 1'
        output = run_command_with_output('list all pending', command, True)
        print(output)
