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
from functools import lru_cache

from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.taxonomy import taxonomy

logger = logging_config.get_logger(__name__)


def get_release_text_file_name(assembly_release_folder, taxonomy_id, assembly_accession, release_text_file_category):
    return os.path.join(assembly_release_folder,
                        "{0}_{1}_{2}.txt".format(taxonomy_id, assembly_accession, release_text_file_category))

@lru_cache
def get_release_folder_name(taxonomy_id):
    return taxonomy.get_normalized_scientific_name_from_ensembl(taxonomy_id)
