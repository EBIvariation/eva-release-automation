import glob
import os
import xml.etree.ElementTree

import requests
from ebi_eva_common_pyutils.logger import logging_config

logger = logging_config.get_logger(__name__)


def get_assembly_name_for_accession(assembly_accession):
    """
    Get assembly name for a given assembly accession using ENA XML API.

    Override value for GenScope12.X since it has a nonsensical name 12X from ENA:
    https://github.com/EBIvariation/eva-accession/blob/56a8d6a5c42a7bf477c4623fcba21085c313bade/eva-accession-release/src/main/java/uk/ac/ebi/eva/accession/release/io/AssemblyNameRetriever.java#L60

    :param assembly_accession: GCA assembly accession
    :return: Assembly name from ENA
    """
    if assembly_accession == "GCA_000003745.2":
        return "Genoscope.12X"
    ENA_XML_API_URL = "https://www.ebi.ac.uk/ena/browser/api/xml/" + assembly_accession
    try:
        response = requests.get(ENA_XML_API_URL)
        if response.ok:
            return xml.etree.ElementTree.fromstring(response.content.decode(encoding="utf-8")) \
                .find("ASSEMBLY").find("NAME").text
        else:
            logger.error("API call {0} failed with response code {1}".format(ENA_XML_API_URL, response.status_code))
    except Exception as ex:
        logger.error(ex)


def make_valid_filename_without_spaces(file_name):
    file_name_without_special_chars = "".join(
        [char for char in file_name if char.isalnum() or char in (' ', '.', '_')]).rstrip()
    return file_name_without_special_chars.replace(" ", "_")


def add_assembly_name_symlink_for_dir(assembly_accession_release_dir):
    parent_dir = os.path.dirname(assembly_accession_release_dir)
    assembly_accession = os.path.basename(assembly_accession_release_dir)
    target_dir_to_create = os.path.join(
        parent_dir, make_valid_filename_without_spaces(get_assembly_name_for_accession(assembly_accession)))
    logger.info(f'Creating symlink: "{target_dir_to_create}" to {assembly_accession_release_dir}')
    if os.path.islink(target_dir_to_create) or os.path.exists(target_dir_to_create):
        os.remove(target_dir_to_create)
    os.symlink(assembly_accession, target_dir_to_create)


def create_assembly_name_symlinks(species_dir):
    for assembly_dir in glob.glob(species_dir + os.path.sep + "GCA_*"):
        add_assembly_name_symlink_for_dir(assembly_dir)


