CREATE USER metadata_user WITH PASSWORD 'metadata_pass';
CREATE DATABASE metadata;
-- Connect to the database
\c metadata

-----------------------------------------eva_progress_tracker-----------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------

CREATE SCHEMA eva_progress_tracker AUTHORIZATION metadata_user;

CREATE TABLE eva_progress_tracker.clustering_release_tracker (
	taxonomy int4 NOT NULL,
	scientific_name text NOT NULL,
	assembly_accession text NOT NULL,
	release_version int8 NOT NULL,
	sources text NOT NULL,
	release_start timestamp NULL,
	release_end timestamp NULL,
	fasta_path text NULL,
	report_path text NULL,
	tempmongo_instance text NULL,
	should_be_released bool NULL,
	num_rs_to_release int8 NULL,
	total_num_variants int8 NULL,
	release_folder_name text NULL,
	release_status text NULL,
	CONSTRAINT clustering_release_tracker_pkey PRIMARY KEY (taxonomy, assembly_accession, release_version)
);

ALTER TABLE eva_progress_tracker.clustering_release_tracker OWNER TO metadata_user;
GRANT ALL ON TABLE eva_progress_tracker.clustering_release_tracker TO metadata_user;

GRANT ALL ON SCHEMA eva_progress_tracker TO metadata_user;


------------------- Permission on Schema
GRANT ALL ON SCHEMA eva_progress_tracker TO metadata_user;
GRANT USAGE, SELECT ON ALL SEQUENCES IN SCHEMA eva_progress_tracker TO metadata_user;

ALTER DATABASE metadata SET search_path TO eva_progress_tracker, public, "$user";
------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------

