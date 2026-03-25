import gzip
import os
import tempfile
from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import patch

from publish_release_to_ftp.publish_release_files_to_ftp import (
    copy_current_assembly_data_to_ftp,
    copy_current_unmapped_files,
    hardlink_previous_unmapped_files,
    get_folder_path_for_assembly,
    get_folder_path_for_species,
    get_folder_path_for_species_assembly,
    get_release_file_list_for_assembly,
    hardlink_to_previous_release_assembly_files_in_ftp,
    publish_assembly_release_files_to_ftp,
    publish_release_top_level_files_to_ftp,
    release_top_level_files_to_copy,
)

TAXONOMY = 9606
ASSEMBLY_ACCESSION = 'GCA_000001405.15'
SPECIES_FOLDER = 'homo_sapiens'
RELEASE_VERSION = 5


def make_assembly_info(should_be_released=True, num_rs_to_release=10):
    return {
        'taxonomy': TAXONOMY,
        'assembly_accession': ASSEMBLY_ACCESSION,
        'release_version': RELEASE_VERSION,
        'release_folder_name': SPECIES_FOLDER,
        'should_be_released': should_be_released,
        'num_rs_to_release': num_rs_to_release,
    }


def make_release_properties(tmp_dir, release_version=RELEASE_VERSION):
    ftp_base = os.path.join(tmp_dir, 'ftp')
    return SimpleNamespace(
        release_version=release_version,
        staging_release_folder=os.path.join(tmp_dir, 'staging'),
        public_ftp_release_base_folder=ftp_base,
        public_ftp_current_release_folder=os.path.join(ftp_base, f'release_{release_version}'),
        public_ftp_previous_release_folder=os.path.join(ftp_base, f'release_{release_version - 1}'),
    )


class TestPublishReleaseFilesToFTP(TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.release_properties = make_release_properties(self.tmp)
        self.assembly_info = make_assembly_info()

    def tearDown(self):
        self._tmp.cleanup()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _create_previous_assembly_folder(self):
        """Create previous-release assembly folder with all expected files."""
        prev_folder = get_folder_path_for_assembly(
            self.release_properties.public_ftp_previous_release_folder, ASSEMBLY_ACCESSION)
        os.makedirs(prev_folder, exist_ok=True)
        for filename in get_release_file_list_for_assembly(self.assembly_info) + ['md5checksums.txt']:
            with open(os.path.join(prev_folder, filename), 'w') as f:
                f.write('Previous release')
        return prev_folder

    def _create_staging_assembly_folder(self):
        """Create staging folder with all expected files for an assembly."""
        staging = os.path.join(
            self.release_properties.staging_release_folder, SPECIES_FOLDER, ASSEMBLY_ACCESSION)
        os.makedirs(staging, exist_ok=True)
        for filename in get_release_file_list_for_assembly(self.assembly_info):
            with open(os.path.join(staging, filename), 'w') as f:
                f.write('New release')
        return staging

    # ------------------------------------------------------------------
    # publish_release_top_level_files_to_ftp
    # ------------------------------------------------------------------

    def test_publish_release_top_level_files_to_ftp(self):
        dest_dir = self.release_properties.public_ftp_current_release_folder
        os.makedirs(dest_dir)

        publish_release_top_level_files_to_ftp(self.release_properties)

        for filename in release_top_level_files_to_copy:
            self.assertTrue(os.path.exists(os.path.join(dest_dir, filename)),
                            f"Expected {filename} in FTP release folder")

    # ------------------------------------------------------------------
    # hardlink_to_previous_release_assembly_files_in_ftp
    # ------------------------------------------------------------------

    def test_hardlink_to_previous_release_assembly_files_in_ftp(self):
        prev_folder = self._create_previous_assembly_folder()

        hardlink_to_previous_release_assembly_files_in_ftp(self.assembly_info, self.release_properties)

        curr_folder = get_folder_path_for_assembly(
            self.release_properties.public_ftp_current_release_folder, ASSEMBLY_ACCESSION)
        self.assertTrue(os.path.isdir(curr_folder))
        for filename in get_release_file_list_for_assembly(self.assembly_info):
            curr_file = os.path.join(curr_folder, filename)
            prev_file = os.path.join(prev_folder, filename)
            self.assertTrue(os.path.exists(curr_file), f"Expected {filename} to be hard-linked")
            self.assertEqual(os.stat(curr_file).st_ino, os.stat(prev_file).st_ino,
                             f"{filename} should share inode with previous release copy")

    def test_hardlink_raises_when_previous_release_missing(self):
        with self.assertRaises(Exception) as ctx:
            hardlink_to_previous_release_assembly_files_in_ftp(self.assembly_info, self.release_properties)
        self.assertIn('does not exist', str(ctx.exception))

    # ------------------------------------------------------------------
    # copy_current_assembly_data_to_ftp
    # ------------------------------------------------------------------

    @patch('publish_release_to_ftp.publish_release_files_to_ftp.get_release_folder_name',
           return_value=SPECIES_FOLDER)
    def test_copy_current_assembly_data_to_ftp(self, _mock_grf):
        self._create_staging_assembly_folder()
        dest_folder = get_folder_path_for_species_assembly(
            self.release_properties.public_ftp_current_release_folder, SPECIES_FOLDER, ASSEMBLY_ACCESSION)
        os.makedirs(dest_folder, exist_ok=True)

        copy_current_assembly_data_to_ftp(self.assembly_info, self.release_properties, dest_folder)

        for filename in get_release_file_list_for_assembly(self.assembly_info):
            self.assertTrue(os.path.exists(os.path.join(dest_folder, filename)),
                            f"Expected {filename} in dest folder")
        md5file = os.path.join(dest_folder, 'md5checksums.txt')
        self.assertTrue(os.path.exists(md5file))
        self.assertGreater(os.path.getsize(md5file), 0)

    # ------------------------------------------------------------------
    # copy_current_unmapped_files
    # ------------------------------------------------------------------

    def test_copy_current_unmapped_files(self):
        source = os.path.join(self.tmp, 'source_species')
        dest = os.path.join(self.tmp, 'dest_species')
        os.makedirs(source)
        os.makedirs(dest)
        # Create the three expected species-level files
        unmapped_file = os.path.join(source, f'{TAXONOMY}_unmapped_ids.txt.gz')
        for path in [unmapped_file,
                     os.path.join(source, 'md5checksums.txt'),
                     os.path.join(source, 'README_unmapped_rs_ids_count.txt')]:
            open(path, 'w').close()

        copy_current_unmapped_files(source, dest)

        self.assertTrue(os.path.exists(os.path.join(dest, f'{TAXONOMY}_unmapped_ids.txt.gz')))

    def test_hardlink_previous_unmapped_files(self):
        prev_species = os.path.join(self.tmp, 'prev', SPECIES_FOLDER)
        curr_species = os.path.join(self.tmp, 'curr', SPECIES_FOLDER)
        os.makedirs(prev_species)
        os.makedirs(curr_species)
        unmapped_src = os.path.join(prev_species, f'{TAXONOMY}_unmapped_ids.txt.gz')
        with gzip.open(unmapped_src, 'wt') as f:
            f.write('RS123\t1\n')

        hardlink_previous_unmapped_files(prev_species, curr_species)

        expected_name = f'{SPECIES_FOLDER}_unmapped_ids.txt.gz'
        linked_file = os.path.join(curr_species, expected_name)
        self.assertTrue(os.path.exists(linked_file), "Hard-linked file should exist in current species folder")
        self.assertEqual(os.stat(unmapped_src).st_ino, os.stat(linked_file).st_ino,
                         "File should be a hard link (same inode)")
        self.assertTrue(os.path.exists(os.path.join(curr_species, 'md5checksums.txt')))
        self.assertTrue(os.path.exists(os.path.join(curr_species, 'README_unmapped_rs_ids_count.txt')))

    def test_hardlink_previous_unmapped_files_no_file(self):
        source = os.path.join(self.tmp, 'empty_source')
        dest = os.path.join(self.tmp, 'dest')
        os.makedirs(source)
        os.makedirs(dest)
        with self.assertRaises(AssertionError):
            hardlink_previous_unmapped_files(source, dest)

    # ------------------------------------------------------------------
    # publish_assembly_release_files_to_ftp
    # ------------------------------------------------------------------

    @patch('publish_release_to_ftp.publish_release_files_to_ftp.get_release_folder_name',
           return_value=SPECIES_FOLDER)
    def test_publish_assembly_release_files_to_ftp_released(self, _mock_grf):
        self._create_staging_assembly_folder()
        current_release = self.release_properties.public_ftp_current_release_folder
        assembly_folder = get_folder_path_for_assembly(current_release, ASSEMBLY_ACCESSION)
        species_folder = get_folder_path_for_species(current_release, SPECIES_FOLDER)
        os.makedirs(assembly_folder, exist_ok=True)
        os.makedirs(species_folder, exist_ok=True)

        publish_assembly_release_files_to_ftp(
            self.assembly_info, self.release_properties, assembly_folder, SPECIES_FOLDER)

        species_assembly_folder = get_folder_path_for_species_assembly(
            current_release, SPECIES_FOLDER, ASSEMBLY_ACCESSION)
        self.assertTrue(os.path.isdir(species_assembly_folder))
        for filename in get_release_file_list_for_assembly(self.assembly_info):
            self.assertTrue(os.path.exists(os.path.join(species_assembly_folder, filename)),
                            f"Expected {filename} in species-assembly folder")

    def test_publish_assembly_release_files_to_ftp_not_released(self):
        prev_folder = self._create_previous_assembly_folder()
        current_release = self.release_properties.public_ftp_current_release_folder
        assembly_folder = get_folder_path_for_assembly(current_release, ASSEMBLY_ACCESSION)
        assembly_info = make_assembly_info(should_be_released=False, num_rs_to_release=0)

        publish_assembly_release_files_to_ftp(
            assembly_info, self.release_properties, assembly_folder, SPECIES_FOLDER)

        curr_assembly_folder = get_folder_path_for_assembly(current_release, ASSEMBLY_ACCESSION)
        self.assertTrue(os.path.isdir(curr_assembly_folder))
        for filename in get_release_file_list_for_assembly(assembly_info):
            curr_file = os.path.join(curr_assembly_folder, filename)
            self.assertTrue(os.path.exists(curr_file), f"{filename} should exist via hardlink")
            self.assertGreater(os.stat(curr_file).st_nlink, 1,
                               f"{filename} should have nlink > 1 (hard link)")



    @patch('publish_release_to_ftp.publish_release_files_to_ftp.get_release_folder_name',
           return_value=SPECIES_FOLDER)
    @patch('publish_release_to_ftp.publish_release_files_to_ftp.create_assembly_name_symlinks')
    def test_full_publish_flow_released_assembly(self, _mock_symlinks, _mock_grf):
        self._create_staging_assembly_folder()
        self._create_previous_assembly_folder()

        current_release = self.release_properties.public_ftp_current_release_folder
        assembly_folder = get_folder_path_for_assembly(current_release, ASSEMBLY_ACCESSION)
        species_folder = get_folder_path_for_species(current_release, SPECIES_FOLDER)
        os.makedirs(assembly_folder, exist_ok=True)
        os.makedirs(species_folder, exist_ok=True)

        publish_assembly_release_files_to_ftp(
            self.assembly_info, self.release_properties, assembly_folder, SPECIES_FOLDER)

        # Species-assembly folder populated
        species_assembly_folder = get_folder_path_for_species_assembly(
            current_release, SPECIES_FOLDER, ASSEMBLY_ACCESSION)
        self.assertTrue(os.path.isdir(species_assembly_folder),
                        "by_species/<species>/<assembly> folder should exist")
        for filename in get_release_file_list_for_assembly(self.assembly_info):
            self.assertTrue(os.path.exists(os.path.join(species_assembly_folder, filename)),
                            f"Expected {filename} in species-assembly folder")
        md5file = os.path.join(species_assembly_folder, 'md5checksums.txt')
        self.assertTrue(os.path.exists(md5file), "md5checksums.txt should be written")
        self.assertGreater(os.path.getsize(md5file), 0, "md5checksums.txt should be non-empty")

        # Assembly-level folder created
        self.assertTrue(os.path.isdir(assembly_folder),
                        "by_assembly/<assembly> folder should exist")

    @patch('publish_release_to_ftp.publish_release_files_to_ftp.create_assembly_name_symlinks')
    def test_full_publish_flow_unchanged_assembly(self, _mock_symlinks):
        prev_folder = self._create_previous_assembly_folder()
        current_release = self.release_properties.public_ftp_current_release_folder
        assembly_folder = get_folder_path_for_assembly(current_release, ASSEMBLY_ACCESSION)
        assembly_info = make_assembly_info(should_be_released=False, num_rs_to_release=0)

        publish_assembly_release_files_to_ftp(
            assembly_info, self.release_properties, assembly_folder, SPECIES_FOLDER)

        curr_assembly_folder = get_folder_path_for_assembly(current_release, ASSEMBLY_ACCESSION)
        self.assertTrue(os.path.isdir(curr_assembly_folder))
        for filename in get_release_file_list_for_assembly(assembly_info):
            curr_file = os.path.join(curr_assembly_folder, filename)
            prev_file = os.path.join(prev_folder, filename)
            self.assertTrue(os.path.exists(curr_file), f"{filename} should be hard-linked")
            self.assertEqual(os.stat(curr_file).st_ino, os.stat(prev_file).st_ino,
                             f"{filename} should share inode with previous release (hard link)")
            self.assertGreater(os.stat(curr_file).st_nlink, 1,
                               f"{filename} nlink should be > 1")
