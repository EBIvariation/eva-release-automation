import os
import shutil
import time
from unittest import TestCase

import yaml
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_internal_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_internal_pyutils.pg_utils import execute_query

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
        command = 'docker compose up --build -d --wait'
        run_command_with_output('Start docker-compose in the background', command)
        # Wait for the postgres to be ready to accept command
        time.sleep(10)


    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.compose_dir)
        command = 'docker compose down'
        run_command_with_output('Stop docker-compose running in the background', command)

    def setUp(self):
        self.output_directory =  os.path.join(resource_folder, 'release_output_dir', 'release_1')

    def tearDown(self):
        if os.path.isdir(self.output_directory):
            command = 'docker exec executor rm -rf /usr/local/test_eva_release/*'
            run_command_with_output('Remove the files in the volume', command)


    def test_run_release_for_species(self):
        add_data_to_docker()
        expected_output = ('| taxonomy | assembly_accession | release_version | release_status |\n'
                           '|     4530 |    GCA_000005425.2 |               1 |        Pending |\n')
        expected_log = os.path.join(self.output_directory, 'oryza_sativa', 'GCA_000005425.2',
                                    'release_active_rs_4530_GCA_000005425.2_1.log')
        assert not os.path.isfile(expected_log)
        command = 'docker exec executor python3 -m  release_automation.run_release_for_species --list_status Pending'
        output = run_command_with_output('list all pending', command, True)
        assert expected_output == output
        command = 'docker exec executor python3 -m release_automation.run_release_for_species --taxonomy_id 4530 --assembly_accessions GCA_000005425.2 --release_version 1'
        output = run_command_with_output('list all pending', command, True)
        # Check that the run has completed and that log files are present
        assert os.path.isfile(expected_log)


def add_data_to_docker():
    config_xml = os.path.join(resource_folder, 'config_xml_file.xml')
    profile = 'localhost'
    inventory_table = 'eva_progress_tracker.clustering_release_tracker'
    # Path inside the executor Docker image
    fasta_file = '/opt/tests/release_automation/resources/GCA_000005425.2.fa'
    report_file = '/opt/tests/release_automation/resources/GCA_000005425.2_assembly_report.txt'
    with get_metadata_connection_handle(profile, config_xml) as metadata_connection_handle:
        insert_query = (f"insert into {inventory_table} "
                        f"(taxonomy, scientific_name, assembly_accession, release_version, sources, fasta_path, "
                        f"report_path, should_be_released, num_rs_to_release, total_num_variants, release_folder_name) "
                        f"values (4530, 'Oryza sativa', 'GCA_000005425.2', 1, 'EVA', '{fasta_file}', "
                        f"'{report_file}', True, 1, 1, 'oryza_sativa')")
        execute_query(metadata_connection_handle, insert_query)


if __name__ == '__main__':
    add_data_to_docker()

