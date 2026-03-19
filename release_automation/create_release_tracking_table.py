#!/usr/bin/env python
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
import argparse
from collections import defaultdict
from functools import cached_property
from itertools import cycle

from ebi_eva_common_pyutils.assembly import NCBIAssembly
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.logger import logging_config, AppLogger
from ebi_eva_common_pyutils.taxonomy.taxonomy import normalise_taxon_scientific_name, get_scientific_name_from_ensembl
from ebi_eva_internal_pyutils.config_utils import get_mongo_uri_for_eva_profile
from ebi_eva_internal_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_internal_pyutils.mongodb import MongoDatabase
from ebi_eva_internal_pyutils.pg_utils import get_all_results_for_query, execute_query

from release_automation.release_config import load_config


class ReleaseTracker(AppLogger):

    def __init__(self, release_version):
        self.private_config_xml_file = cfg.query("maven", "settings_file")
        self.maven_profile = cfg.query("maven", "environment")
        self.release_version = release_version
        self.ref_dir = cfg.query('genome_downloader', 'output_directory')

    @cached_property
    def metadata_conn(self):
        return get_metadata_connection_handle(self.maven_profile, self.private_config_xml_file)

    @cached_property
    def mongo_conn(self):
        mongo_uri = get_mongo_uri_for_eva_profile(self.maven_profile, self.private_config_xml_file)
        return MongoDatabase(uri=mongo_uri, db_name="eva_accession_sharded")

    def create_table_if_not_exists(self):
        query_create_table = (
            'create table if not exists eva_progress_tracker.clustering_release_tracker('
            'taxonomy int4 not null, '
            'scientific_name text not null, '
            'assembly_accession text not null, '
            'release_version int8 not null, '
            'sources text not null,'
            'clustering_status text null, '  # unused
            'clustering_start timestamp null, '  # unused
            'clustering_end timestamp null, '  # unused
            'should_be_clustered boolean null, '  # unused
            'fasta_path text null, '
            'report_path text null, '
            'tempmongo_instance text null, '
            'should_be_released boolean null, '
            'num_rs_to_release int8 null, '  # not computed but still used by release automation
            'total_num_variants int8 null, '  # not computed and unused
            'release_folder_name text null, '
            'release_status text null, '
            'primary key (taxonomy, assembly_accession, release_version))'
        )
        execute_query(self.metadata_conn, query_create_table)

    def fill_release_entries(self):
        """Fill in release table based on previous release data, EVA metadata, and supported assembly tracker.
        Also fills in should_be_released values."""
        self._fill_from_previous_release()
        self._fill_from_eva_metadata()
        self._fill_from_supported_assembly_tracker()
        self.fill_should_be_released_from_clustered_variant_update()

    def _fill_from_previous_release(self):
        query = f"""select taxonomy, scientific_name, assembly_accession, sources, fasta_path, report_path, 
                    release_folder_name from eva_progress_tracker.clustering_release_tracker 
                    where release_version = {self.release_version - 1}"""
        for tax, sc_name, asm_acc, src, fs_path, rpt_path, rls_folder_name in get_all_results_for_query(
                self.metadata_conn, query):
            self._insert_entry_for_taxonomy_and_assembly(tax, asm_acc, src, sc_name, fs_path, rpt_path,
                                                         rls_folder_name)

    def _fill_from_eva_metadata(self):
        query = f"""select distinct  pt.taxonomy_id as taxonomy, asm.assembly_accession as assembly_accession
                    from evapro.project_taxonomy pt
                    join evapro.project_analysis pa on pt.project_accession = pa.project_accession 
                    join evapro.analysis a on a.analysis_accession = pa.analysis_accession 
                    join evapro.assembly asm on asm.assembly_set_id = a.assembly_set_id 
                    and asm.assembly_accession is not null and assembly_accession like 'GCA%'"""
        sources = 'EVA'
        for tax, asm_acc in get_all_results_for_query(self.metadata_conn, query):
            self._insert_entry_for_taxonomy_and_assembly(tax, asm_acc, sources)

    def _fill_from_supported_assembly_tracker(self):
        query = f"""select distinct taxonomy_id as taxonomy, assembly_id as assembly_accession
                    from evapro.supported_assembly_tracker"""
        sources = 'DBSNP, EVA'
        for tax, asm_acc in get_all_results_for_query(self.metadata_conn, query):
            self._insert_entry_for_taxonomy_and_assembly(tax, asm_acc, sources)


    def _insert_entry_for_taxonomy_and_assembly(self, tax, asm_acc, sources, sc_name=None, fasta_path=None,
                                                report_path=None, release_folder_name=None):
        sc_name = sc_name if sc_name else get_scientific_name_from_ensembl(tax)
        sc_name = sc_name.replace("'", "\''")
        if asm_acc != 'Unmapped':
            ncbi_assembly = NCBIAssembly(asm_acc, sc_name, self.ref_dir)
            fasta_path = fasta_path if fasta_path else ncbi_assembly.assembly_fasta_path
            report_path = report_path if report_path else ncbi_assembly.assembly_report_path
        release_folder_name = release_folder_name if release_folder_name else normalise_taxon_scientific_name(sc_name)

        tempmongo_instance = 'dummy'
        src_in_db = self.get_sources_for_taxonomy_assembly(tax, asm_acc)

        if not src_in_db:
            # entry does not exist for tax and asm
            insert_query = f"""INSERT INTO eva_progress_tracker.clustering_release_tracker(
                            taxonomy, scientific_name, assembly_accession, release_version, sources,
                            fasta_path, report_path, tempmongo_instance, release_folder_name) 
                            VALUES ({tax}, '{sc_name}', '{asm_acc}', {self.release_version}, '{sources}', 
                            '{fasta_path}', '{report_path}', '{tempmongo_instance}', '{release_folder_name}') 
                            ON CONFLICT DO NOTHING"""
            execute_query(self.metadata_conn, insert_query)
        else:
            # if DB source is equal to what we are trying to insert or if the DB source already contains
            # both EVA and DBSNP, no need to insert again
            if src_in_db == sources or ('EVA' in src_in_db and 'DBSNP' in src_in_db):
                self.info(f"Entry already present for taxonomy {tax} and assembly {asm_acc} with sources {sources}")
            else:
                # We have different sources which means we need to update entry to have both DBNSP and EVA in sources
                update_query = f"""update eva_progress_tracker.clustering_release_tracker set sources='DBSNP, EVA'
                                where taxonomy={tax} and assembly_accession='{asm_acc}' and  
                                release_version={self.release_version}"""
                execute_query(self.metadata_conn, update_query)

    def fill_should_be_released_from_clustered_variant_update(self):
        # Old date to include everything unless we find a previous release
        last_release_date = '2000-01-01'

        query = f'''select release_date from eva_stats.release_rs where release_version={self.release_version -1};'''
        results = list(get_all_results_for_query(self.metadata_conn, query))
        if results:
            last_release_date, last_release_version = results[0]
            assert int(last_release_version) == self.release_version - 1, f'Attempting to prepare release {self.release_version} but previous version is {last_release_version}'

        query = f"""select distinct taxonomy_id, assembly_accession from clustered_variant_update 
                where ingestion_time>'{last_release_date}'; """
        for tax, asm_acc in get_all_results_for_query(self.metadata_conn, query):
            self.fill_should_be_released_for_taxid_assembly(tax, asm_acc)

    def fill_should_be_released_for_taxid_assembly(self, tax, asm_acc):
        should_be_released = True
        num_rs_to_release = 1
        self.info(f"For taxonomy {tax} and assembly {asm_acc}, should_be_released is {should_be_released}")
        update_should_be_released_query = f"""update eva_progress_tracker.clustering_release_tracker 
                                set should_be_released={should_be_released}, num_rs_to_release={num_rs_to_release}
                                where taxonomy={tax} and assembly_accession='{asm_acc}' and release_version={self.release_version}"""
        execute_query(self.metadata_conn, update_should_be_released_query)


    def get_taxonomy_list_for_release(self):
        """Get all taxonomies with assemblies and sources for the current release version."""
        tax_asm = defaultdict(defaultdict)
        query = f"""select distinct taxonomy, assembly_accession, sources 
                    from eva_progress_tracker.clustering_release_tracker
                    where release_version={self.release_version}"""
        for tax, asm_acc, sources in get_all_results_for_query(self.metadata_conn, query):
            tax_asm[tax][asm_acc] = sources
        return tax_asm

    def get_assemblies_and_sources_for_taxonomy(self, taxonomy):
        assembly_source = {}
        query = f"""SELECT distinct assembly_accession, sources from eva_progress_tracker.clustering_release_tracker 
                    where taxonomy = {taxonomy} and release_version = {self.release_version}"""
        for assembly, sources in get_all_results_for_query(self.metadata_conn, query):
            assembly_source[assembly] = sources
        return assembly_source

    def get_sources_for_taxonomy_assembly(self, taxonomy, assembly):
        query = f"""SELECT sources from eva_progress_tracker.clustering_release_tracker 
                    where taxonomy = {taxonomy} and assembly_accession='{assembly}' 
                    and release_version = {self.release_version}"""
        result = get_all_results_for_query(self.metadata_conn, query)
        if not result:
            return None
        else:
            return result[0][0]


def main():
    parser = argparse.ArgumentParser(description='Create and load the clustering and release tracking table',
                                     add_help=False)
    parser.add_argument("--release-version", required=True, type=int, help="version of the release")

    args = parser.parse_args()
    load_config()
    logging_config.add_stdout_handler()

    release_tracker = ReleaseTracker(release_version=args.release_version)

    release_tracker.create_table_if_not_exists()
    release_tracker.fill_release_entries()




if __name__ == '__main__':
    main()
