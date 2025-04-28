FROM postgres:9.6

COPY metadata.sql /docker-entrypoint-initdb.d/
COPY accession_jt.sql /docker-entrypoint-initdb.d/
COPY variant_jt.sql /docker-entrypoint-initdb.d/
