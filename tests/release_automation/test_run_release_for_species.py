import json
import os
import time
from datetime import datetime
from unittest import TestCase

from bson import  Int64
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_internal_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_internal_pyutils.mongo_utils import get_mongo_connection_handle
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
            # Remove from within Docker to ensure we have permission
            command = 'docker exec executor rm -rf /usr/local/test_eva_release/*'
            run_command_with_output('Remove the files in the volume', command)

    def test_run_release_for_species(self):
        add_data_to_docker()
        expected_output = ('| taxonomy | assembly_accession | release_version | release_status |\n'
                           '|     4530 |    GCA_000005425.2 |               1 |        Pending |\n')
        expected_files = [
            os.path.join(self.output_directory, 'oryza_sativa', 'GCA_000005425.2', f)
            for f in [
                'log_files/release_active_rs_4530_GCA_000005425.2_1.log',
                '4530_GCA_000005425.2_current_ids.vcf.gz',
                '4530_GCA_000005425.2_deprecated_ids.txt.gz',
                '4530_GCA_000005425.2_merged_ids.vcf.gz',
                'README_rs_ids_counts.txt',
            ]
        ]
        command = 'docker exec executor python3 -m  release_automation.run_release_for_species --list_status Pending'
        output = run_command_with_output('list all pending', command, True)
        assert expected_output == output
        command = 'docker exec executor python3 -m release_automation.run_release_for_species --taxonomy_id 4530 --assembly_accessions GCA_000005425.2 --release_version 1'
        output = run_command_with_output('Run release for species', command, True)
        # Check that the run has completed and that log files are present
        for expected_file in expected_files:
            assert os.path.isfile(expected_file)

        # Check that the counts are correct
        count_file = expected_files[-1]
        counts = {}
        with open(count_file) as open_file:
            for line in open_file:
                sp_line = line.strip().split()
                if sp_line:
                    counts[sp_line[0]] = sp_line[1]
        assert counts == {
            '4530_GCA_000005425.2_current_ids.vcf.gz':'6',
            '4530_GCA_000005425.2_deprecated_ids.txt.gz':'1',
            '4530_GCA_000005425.2_merged_ids.vcf.gz':'1'
        }

        command = 'docker exec executor python3 -m  release_automation.run_release_for_species --list_status Completed'
        output = run_command_with_output('list all pending', command, True)
        expected_output = ('| taxonomy | assembly_accession | release_version | release_status |\n'
                           '|     4530 |    GCA_000005425.2 |               1 |      Completed |\n')
        assert expected_output == output



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
    sve_file = os.path.join(resource_folder, 'submittedVariantEntities.json')
    cve_file = os.path.join(resource_folder, 'clusteredVariantEntities.json')
    svoe_file = os.path.join(resource_folder, 'submittedVariantOperationEntities.json')
    cvoe_file = os.path.join(resource_folder, 'clusteredVariantOperationEntities.json')

    sves = read_mongo_data(sve_file)
    cves = read_mongo_data(cve_file)
    svoe = read_mongo_data(svoe_file)
    cvoe = read_mongo_data(cvoe_file)
    with get_mongo_connection_handle(profile, config_xml) as mongo_handle:
        db = mongo_handle['eva_accession_sharded']
        db['dbsnpSubmittedVariantEntity'].insert_many(sves)
        db['dbsnpClusteredVariantEntity'].insert_many(cves)
        db['dbsnpSubmittedVariantOperationEntity'].insert_many(svoe)
        db['dbsnpClusteredVariantOperationEntity'].insert_many(cvoe)

def read_mongo_data(json_file):
    with open(json_file) as open_file:
        json_data = json.load(open_file)
    for document in json_data:
        add_mongo_types(document)
        if 'inactiveObjects' in document:
            for sub_doc in document['inactiveObjects']:
                add_mongo_types(sub_doc)
    return json_data

def add_mongo_types(json_dict):
    for key in json_dict:
        if key in ['start', 'rs', 'accession', 'mergeInto']:
            json_dict[key] = Int64(json_dict[key])
        if key in ['createdDate']:
            json_dict[key] = datetime.strptime(json_dict[key], '%Y-%m-%dT%H:%M:%S.%fZ')

def clean_up_docker_data():
    config_xml = os.path.join(resource_folder, 'config_xml_file.xml')
    profile = 'localhost'
    inventory_table = 'eva_progress_tracker.clustering_release_tracker'

    with get_metadata_connection_handle(profile, config_xml) as metadata_connection_handle:
        delete_query = (f"delete from {inventory_table};")
        execute_query(metadata_connection_handle, delete_query)
    with get_mongo_connection_handle(profile, config_xml) as mongo_handle:
        db = mongo_handle['eva_accession_sharded']
        db['dbsnpSubmittedVariantEntity'].drop()
        db['dbsnpClusteredVariantEntity'].drop()
        db['dbsnpSubmittedVariantOperationEntity'].drop()
        db['dbsnpClusteredVariantOperationEntity'].drop()

if __name__ == '__main__':
    clean_up_docker_data()
    add_data_to_docker()

