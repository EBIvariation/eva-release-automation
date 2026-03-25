import os
import re
import sqlite3
from datetime import date
from unittest import TestCase
from unittest.mock import patch, MagicMock

from release_automation.create_release_tracking_table import ReleaseTracker

resource_folder = os.path.join(os.path.dirname(__file__), 'resources')


def make_tracker(release_version=2):
    """Create a ReleaseTracker with cfg mocked out."""
    with patch('release_automation.create_release_tracking_table.cfg') as mock_cfg:
        mock_cfg.query.side_effect = lambda *args: {
            ('maven', 'settings_file'): 'settings.xml',
            ('maven', 'environment'): 'localhost',
            ('genome_downloader', 'output_directory'): '/tmp/ref',
        }.get(args)
        tracker = ReleaseTracker(release_version=release_version)
    # Inject a mock metadata connection so cached_property is never evaluated
    tracker.__dict__['metadata_conn'] = MagicMock()
    return tracker


def assert_no_multispace(expected, actual):
    expected = re.sub(r'\s+', ' ', expected)
    actual = re.sub(r'\s+', ' ', actual)
    assert expected == actual, '\n' + actual + '\n' + expected


class TestReleaseTracker(TestCase):

    def setUp(self):
        self.tracker = make_tracker(release_version = 2)

    @patch('release_automation.create_release_tracking_table.execute_query')
    @patch('release_automation.create_release_tracking_table.get_all_results_for_query', return_value=[])
    @patch('release_automation.create_release_tracking_table.NCBIAssembly')
    @patch('release_automation.create_release_tracking_table.get_scientific_name_from_ensembl', return_value='Oryza sativa')
    def test_insert_new_entry(self, mock_sci_name, mock_ncbi, mock_get_all, mock_execute):
        mock_ncbi.return_value.assembly_fasta_path = '/ref/GCA_000005425.2.fa'
        mock_ncbi.return_value.assembly_report_path = '/ref/GCA_000005425.2_report.txt'

        self.tracker._insert_entry_for_taxonomy_and_assembly(4530, 'GCA_000005425.2')

        mock_execute.assert_called_once()
        query = mock_execute.call_args[0][1]
        expected_query = """INSERT INTO eva_progress_tracker.clustering_release_tracker(
        taxonomy, scientific_name, assembly_accession, release_version, sources,
        fasta_path, report_path, tempmongo_instance, release_folder_name) 
        VALUES (4530, 'Oryza sativa', 'GCA_000005425.2', 2, 'EVA, DBSNP', 
        '/ref/GCA_000005425.2.fa', '/ref/GCA_000005425.2_report.txt', 'dummy', 'oryza_sativa') 
        ON CONFLICT DO NOTHING"""
        assert_no_multispace(query, expected_query)


    @patch('release_automation.create_release_tracking_table.execute_query')
    @patch('release_automation.create_release_tracking_table.get_all_results_for_query')
    def test_fill_from_previous_release(self, mock_get_all, mock_execute):
        # First call: previous-release query; second call: source-check inside _insert_entry
        mock_get_all.side_effect = [
            [(4530, 'Oryza sativa', 'GCA_000005425.2', '/p/fasta', '/p/report', 'oryza_sativa')],
            [],  # entry not yet in current release
        ]

        self.tracker._fill_from_previous_release()
        assert mock_get_all.call_count == 2
        query = mock_get_all.call_args_list[0][0][1]
        expected_query = """select taxonomy, scientific_name, assembly_accession, fasta_path, report_path,
                    release_folder_name from eva_progress_tracker.clustering_release_tracker
                    where release_version = 1"""
        assert_no_multispace(query, expected_query)
        expected_query = """SELECT sources from eva_progress_tracker.clustering_release_tracker 
        where taxonomy = 4530 and assembly_accession='GCA_000005425.2' and release_version = 2"""
        query = mock_get_all.call_args_list[1][0][1]
        assert_no_multispace(query, expected_query)

        mock_execute.assert_called_once()
        query = mock_execute.call_args[0][1]
        expected_query = """INSERT INTO eva_progress_tracker.clustering_release_tracker(
        taxonomy, scientific_name, assembly_accession, release_version, sources,
        fasta_path, report_path, tempmongo_instance, release_folder_name) 
        VALUES (4530, 'Oryza sativa', 'GCA_000005425.2', 2, 'EVA, DBSNP', 
        '/p/fasta', '/p/report', 'dummy', 'oryza_sativa') 
        ON CONFLICT DO NOTHING"""
        assert_no_multispace(query, expected_query)

    @patch('release_automation.create_release_tracking_table.execute_query')
    @patch('release_automation.create_release_tracking_table.get_all_results_for_query')
    @patch('release_automation.create_release_tracking_table.NCBIAssembly')
    @patch('release_automation.create_release_tracking_table.get_scientific_name_from_ensembl', return_value='Oryza sativa')
    def test_fill_from_eva_metadata_inserts_eva_entries(self, mock_sci_name, mock_ncbi, mock_get_all, mock_execute):
        mock_ncbi.return_value.assembly_fasta_path = '/ref/fasta'
        mock_ncbi.return_value.assembly_report_path = '/ref/report'
        mock_get_all.side_effect = [
            [(4530, 'GCA_000005425.2')],  # EVA metadata query
            [],  # source check: entry not yet present
        ]

        self.tracker._fill_from_eva_metadata()

        metadata_query = mock_get_all.call_args_list[0][0][1]
        expected_select = """select distinct  pt.taxonomy_id as taxonomy, asm.assembly_accession as assembly_accession
                    from evapro.project_taxonomy pt
                    join evapro.project_analysis pa on pt.project_accession = pa.project_accession
                    join evapro.analysis a on a.analysis_accession = pa.analysis_accession
                    join evapro.assembly asm on asm.assembly_set_id = a.assembly_set_id
                    and asm.assembly_accession is not null and assembly_accession like 'GCA%'"""
        assert_no_multispace(metadata_query, expected_select)

        insert_query = mock_execute.call_args[0][1]
        expected_insert = """INSERT INTO eva_progress_tracker.clustering_release_tracker(
                            taxonomy, scientific_name, assembly_accession, release_version, sources,
                            fasta_path, report_path, tempmongo_instance, release_folder_name)
                            VALUES (4530, 'Oryza sativa', 'GCA_000005425.2', 2, 'EVA, DBSNP',
                            '/ref/fasta', '/ref/report', 'dummy', 'oryza_sativa')
                            ON CONFLICT DO NOTHING"""
        assert_no_multispace(insert_query, expected_insert)

    @patch('release_automation.create_release_tracking_table.execute_query')
    @patch('release_automation.create_release_tracking_table.get_all_results_for_query')
    @patch('release_automation.create_release_tracking_table.NCBIAssembly')
    @patch('release_automation.create_release_tracking_table.get_scientific_name_from_ensembl', return_value='Oryza sativa')
    def test_fill_from_supported_assembly_tracker(self, mock_sci_name, mock_ncbi, mock_get_all, mock_execute):
        mock_ncbi.return_value.assembly_fasta_path = '/ref/fasta'
        mock_ncbi.return_value.assembly_report_path = '/ref/report'
        mock_get_all.side_effect = [
            [(4530, 'GCA_000005425.2')],  # supported_assembly_tracker query
            [],  # source check: entry not yet present
        ]

        self.tracker._fill_from_supported_assembly_tracker()

        tracker_query = mock_get_all.call_args_list[0][0][1]
        expected_select = """select distinct taxonomy_id as taxonomy, assembly_id as assembly_accession
                    from evapro.supported_assembly_tracker"""
        assert_no_multispace(tracker_query, expected_select)

        insert_query = mock_execute.call_args[0][1]
        expected_insert = """INSERT INTO eva_progress_tracker.clustering_release_tracker(
                            taxonomy, scientific_name, assembly_accession, release_version, sources,
                            fasta_path, report_path, tempmongo_instance, release_folder_name)
                            VALUES (4530, 'Oryza sativa', 'GCA_000005425.2', 2, 'EVA, DBSNP',
                            '/ref/fasta', '/ref/report', 'dummy', 'oryza_sativa')
                            ON CONFLICT DO NOTHING"""
        assert_no_multispace(insert_query, expected_insert)

    @patch('release_automation.create_release_tracking_table.execute_query')
    @patch('release_automation.create_release_tracking_table.get_all_results_for_query')
    def test_fill_should_be_released_from_clustered_variant_update_when_no_previous_release(self, mock_get_all, mock_execute):
        """With no prior release record, use a far-past date to capture all clustered updates."""
        mock_get_all.side_effect = [
            [],  # no rows in eva_stats.release_rs
            [(4530, 'GCA_000005425.2')],
        ]

        self.tracker.fill_should_be_released_from_clustered_variant_update()

        release_rs_query = mock_get_all.call_args_list[0][0][1]
        expected_release_rs_query = """select release_date from eva_stats.release_rs where release_version=1;"""
        assert_no_multispace(release_rs_query, expected_release_rs_query)

        clustered_update_query = mock_get_all.call_args_list[1][0][1]
        expected_clustered_query = """select distinct taxonomy_id, assembly_accession from clustered_variant_update
                where ingestion_time>'2000-01-01'; """
        assert_no_multispace(clustered_update_query, expected_clustered_query)

        executed_query = mock_execute.call_args[0][1]
        expected_update = """update eva_progress_tracker.clustering_release_tracker
                        set should_be_released=True, num_rs_to_release=1
                        where taxonomy=4530 and assembly_accession='GCA_000005425.2' and release_version=2"""
        assert_no_multispace(executed_query, expected_update)

    @patch('release_automation.create_release_tracking_table.execute_query')
    @patch('release_automation.create_release_tracking_table.get_all_results_for_query')
    def test_fill_should_be_released_from_clustered_variant_update_with_previous_release_available(self, mock_get_all, mock_execute):
        """Entries in clustered_variant_update after the previous release date should be marked."""
        mock_get_all.side_effect = [
            [(date(2024, 6, 15), 1)],  # previous release: version 1 on 2024-06-15
            [(9606, 'GCA_000001405.15')],
        ]

        self.tracker.fill_should_be_released_from_clustered_variant_update()

        clustered_update_query = mock_get_all.call_args_list[1][0][1]
        expected_clustered_query = """select distinct taxonomy_id, assembly_accession from clustered_variant_update
                where ingestion_time>'2024-06-15'; """
        assert_no_multispace(clustered_update_query, expected_clustered_query)

        executed_query = mock_execute.call_args[0][1]
        expected_update = """update eva_progress_tracker.clustering_release_tracker
                        set should_be_released=True, num_rs_to_release=1
                        where taxonomy=9606 and assembly_accession='GCA_000001405.15' and release_version=2"""
        assert_no_multispace(executed_query, expected_update)


    @patch('release_automation.create_release_tracking_table.execute_query')
    def test_fill_should_be_released_for_taxid_assembly(self, mock_execute):
        self.tracker.fill_should_be_released_for_taxid_assembly(4530, 'GCA_000005425.2')

        query = mock_execute.call_args[0][1]
        expected_update = """update eva_progress_tracker.clustering_release_tracker
                        set should_be_released=True, num_rs_to_release=1
                        where taxonomy=4530 and assembly_accession='GCA_000005425.2' and release_version=2"""
        assert_no_multispace(query, expected_update)


class SQLiteCompatibleCursor:
    """Wraps sqlite3.Cursor to add the context manager protocol expected by psycopg2 code."""

    def __init__(self, cursor):
        self._cursor = cursor

    def execute(self, query):
        self._cursor.execute(query)

    def fetchall(self):
        return self._cursor.fetchall()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


class SQLiteCompatibleConnection:
    """Wraps sqlite3 in-memory databases to simulate a psycopg2 PostgreSQL connection.

    Each PostgreSQL schema is represented by a separate attached in-memory SQLite database,
    making schema-qualified names (e.g. eva_progress_tracker.clustering_release_tracker) work.
    """

    def __init__(self):
        self._conn = sqlite3.connect(':memory:')
        self._conn.execute("ATTACH DATABASE ':memory:' AS eva_progress_tracker")
        self._conn.execute("ATTACH DATABASE ':memory:' AS evapro")
        self._conn.execute("ATTACH DATABASE ':memory:' AS eva_stats")
        self._conn.commit()
        self._create_tables()

    def _create_tables(self):
        statements = [
            """CREATE TABLE eva_progress_tracker.clustering_release_tracker (
                taxonomy INTEGER NOT NULL,
                scientific_name TEXT NOT NULL,
                assembly_accession TEXT NOT NULL,
                release_version INTEGER NOT NULL,
                sources TEXT NOT NULL,
                clustering_status TEXT,
                clustering_start TEXT,
                clustering_end TEXT,
                should_be_clustered INTEGER,
                fasta_path TEXT,
                report_path TEXT,
                tempmongo_instance TEXT,
                should_be_released INTEGER,
                num_rs_to_release INTEGER,
                total_num_variants INTEGER,
                release_folder_name TEXT,
                release_status TEXT,
                PRIMARY KEY (taxonomy, assembly_accession, release_version)
            )""",
            """CREATE TABLE evapro.project_taxonomy (
                project_accession TEXT,
                taxonomy_id INTEGER
            )""",
            """CREATE TABLE evapro.project_analysis (
                project_accession TEXT,
                analysis_accession TEXT
            )""",
            """CREATE TABLE evapro.analysis (
                analysis_accession TEXT,
                assembly_set_id INTEGER
            )""",
            """CREATE TABLE evapro.assembly (
                assembly_set_id INTEGER,
                assembly_accession TEXT
            )""",
            """CREATE TABLE evapro.supported_assembly_tracker (
                taxonomy_id INTEGER,
                assembly_id TEXT
            )""",
            """CREATE TABLE clustered_variant_update (
                taxonomy_id INTEGER,
                assembly_accession TEXT,
                ingestion_time TEXT
            )""",
            """CREATE TABLE eva_stats.release_rs (
                release_date TEXT,
                release_version INTEGER
            )""",
        ]
        for stmt in statements:
            self._conn.execute(stmt)
        self._conn.commit()

    def cursor(self):
        return SQLiteCompatibleCursor(self._conn.cursor())

    def commit(self):
        self._conn.commit()

    def execute(self, query, params=()):
        """Direct execute for test setup — bypasses the cursor wrapper."""
        self._conn.execute(query, params)

    def fetchall(self, query):
        """Direct query for test assertions."""
        return self._conn.execute(query).fetchall()


class TestReleaseTrackerEndToEnd(TestCase):
    """End-to-end tests using an SQLite in-memory database instead of a real PostgreSQL instance.

    Scenario:
      - Release version 2 is being prepared.
      - Version 1 already has one entry: Oryza sativa / GCA_000005425.2 (EVA-only).
      - EVA metadata adds a new species: Homo sapiens / GCA_000001405.15.
      - The supported assembly tracker lists both assemblies, upgrading both sources to DBSNP+EVA.
      - clustered_variant_update contains updates for both assemblies, so both are marked for release.
    """

    def setUp(self):
        self.conn = SQLiteCompatibleConnection()

        # Previous release row (version 1)
        self.conn.execute(
            """INSERT INTO eva_progress_tracker.clustering_release_tracker
               (taxonomy, scientific_name, assembly_accession, release_version, sources,
                fasta_path, report_path, tempmongo_instance, release_folder_name)
               VALUES (4530, 'Oryza sativa', 'GCA_000005425.2', 1, 'EVA',
                       '/old/fasta', '/old/report', 'dummy', 'oryza_sativa')"""
        )

        # Supported assembly tracker: both assemblies have DBSNP data
        self.conn.execute("INSERT INTO evapro.supported_assembly_tracker VALUES (4530, 'GCA_000005425.2')")
        self.conn.execute("INSERT INTO evapro.supported_assembly_tracker VALUES (9606, 'GCA_000001405.15')")

        # Clustered variant updates for both assemblies (well after epoch date)
        self.conn.execute("INSERT INTO clustered_variant_update VALUES (4530, 'GCA_000005425.2', '2024-01-01')")
        self.conn.execute("INSERT INTO clustered_variant_update VALUES (9606, 'GCA_000001405.15', '2024-06-01')")

        # No previous release record in eva_stats ('2000-01-01' will be used)
        self.conn.commit()

    @patch('release_automation.create_release_tracking_table.NCBIAssembly')
    @patch('release_automation.create_release_tracking_table.get_scientific_name_from_ensembl')
    def test_full_release_flow(self, mock_sci_name, mock_ncbi):
        mock_sci_name.side_effect = lambda tax: {4530: 'Oryza sativa', 9606: 'Homo sapiens'}.get(tax, f'Species {tax}')
        mock_ncbi.return_value.assembly_fasta_path = '/ref/fasta'
        mock_ncbi.return_value.assembly_report_path = '/ref/report'

        tracker = make_tracker(release_version=2)
        tracker.__dict__['metadata_conn'] = self.conn

        tracker.fill_release_entries()

        rows = {
            row[0]: row
            for row in self.conn.fetchall(
                """SELECT taxonomy, assembly_accession, sources, should_be_released
                   FROM eva_progress_tracker.clustering_release_tracker
                   WHERE release_version = 2"""
            )
        }

        assert len(rows) == 2

        oryza = rows[4530]
        assert oryza[1] == 'GCA_000005425.2'
        assert oryza[2] == 'EVA, DBSNP'
        assert oryza[3]  # should_be_released

        human = rows[9606]
        assert human[1] == 'GCA_000001405.15'
        assert human[2] == 'EVA, DBSNP'
        assert human[3]  # should_be_released

    @patch('release_automation.create_release_tracking_table.NCBIAssembly')
    @patch('release_automation.create_release_tracking_table.get_scientific_name_from_ensembl')
    def test_assembly_without_clustered_updates_not_marked_for_release(self, mock_sci_name, mock_ncbi):
        """An assembly present in the tracker but absent from clustered_variant_update should not be released."""
        mock_sci_name.side_effect = lambda tax: {4530: 'Oryza sativa', 9606: 'Homo sapiens',
                                                  10090: 'Mus musculus'}.get(tax, f'Species {tax}')
        mock_ncbi.return_value.assembly_fasta_path = '/ref/fasta'
        mock_ncbi.return_value.assembly_report_path = '/ref/report'

        # Add a third assembly with no clustered update
        self.conn.execute("INSERT INTO evapro.supported_assembly_tracker VALUES (10090, 'GCA_000001635.9')")
        self.conn.commit()

        tracker = make_tracker(release_version=2)
        tracker.__dict__['metadata_conn'] = self.conn

        tracker.fill_release_entries()

        rows = self.conn.fetchall(
            """SELECT taxonomy, should_be_released
               FROM eva_progress_tracker.clustering_release_tracker
               WHERE release_version = 2 ORDER BY taxonomy"""
        )
        released = {row[0]: row[1] for row in rows}

        assert released[4530]    # has clustered update → released
        assert released[9606]    # has clustered update → released
        assert not released[10090]  # no clustered update → not released
