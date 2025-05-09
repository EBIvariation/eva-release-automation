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

import logging
from argparse import ArgumentParser

from ebi_eva_common_pyutils.config import cfg

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_internal_pyutils.metadata_utils import get_metadata_connection_handle

from release_automation.release_utils import update_release_progress_status
from release_automation.release_config import load_config

logger = logging.getLogger(__name__)


def update_release_status_for_assembly(taxonomy_id, assembly_accession, release_version):
    with get_metadata_connection_handle(cfg['maven']['environment'],
                                        cfg['maven']['settings_file']) as metadata_connection_handle:
        release_species_inventory_table = cfg['release']['inventory_table']
        update_release_progress_status(metadata_connection_handle, release_species_inventory_table,
                                       taxonomy_id, assembly_accession, release_version,
                                       release_status='Completed')
        logger.info("Successfully marked release status as 'Completed' in {0} for taxonomy {1} and assembly {2}"
                    .format(release_species_inventory_table, taxonomy_id, assembly_accession))


def main():
    argparse = ArgumentParser()
    argparse.add_argument("--taxonomy_id", help="ex: 9913", required=True)
    argparse.add_argument("--assembly_accession", help="ex: GCA_000003055.6", required=True)
    argparse.add_argument("--release_version", help="ex: 2", type=int, required=True)
    args = argparse.parse_args()
    load_config()

    logging_config.add_stdout_handler()
    update_release_status_for_assembly(args.taxonomy_id, args.assembly_accession, args.release_version)
    return 0

if __name__ == "__main__":
    main()
