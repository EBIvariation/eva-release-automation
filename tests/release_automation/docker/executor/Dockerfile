FROM ubuntu:22.04

ENV JAVA_VERSION=8
ENV NXF_VER=21.10.5
ENV VCF_VALIDATOR_VERSION=0.9.4
ENV CONTIG_NAME_CONVERTER_VERSION=0.0.2

# install utility libs, python and java
RUN apt update
RUN apt-get install -y curl wget maven git bcftools python3-pip python3.10-venv openjdk-${JAVA_VERSION}-jdk
# Set Java Home
ENV JAVA_HOME=/usr/lib/jvm/java-${JAVA_VERSION}-openjdk-amd64
ENV PATH=$JAVA_HOME/bin:$PATH

# install nextflow
RUN curl -L "https://github.com/nextflow-io/nextflow/releases/download/v${NXF_VER}/nextflow-${NXF_VER}-all" | bash \
    && mv nextflow /usr/local/bin/

RUN mkdir /usr/local/software/
RUN mkdir -p /usr/local/test_eva_release/

# install vcf_validator
RUN curl -LJo /usr/local/bin/vcf_validator https://github.com/EBIvariation/vcf-validator/releases/download/v${VCF_VALIDATOR_VERSION}/vcf_validator_linux \
    && chmod 755 /usr/local/bin/vcf_validator

# install vcf_assembly_checker
RUN curl -LJo /usr/local/bin/vcf_assembly_checker https://github.com/EBIvariation/vcf-validator/releases/download/v${VCF_VALIDATOR_VERSION}/vcf_assembly_checker_linux \
    && chmod 755 /usr/local/bin/vcf_assembly_checker

# Install contig_name_converter
RUN pip install -q "git+https://github.com/EBIvariation/contig_name_converter.git@${CONTIG_NAME_CONVERTER_VERSION}"

COPY tests/release_automation/docker/executor/maven.xml /root
COPY tests/release_automation/docker/executor/release_config.yml /root/.release_config.yml
COPY tests/release_automation/docker/executor/nextflow.config /root/.nextflow/config

# Setup release pipeline
WORKDIR /usr/local/software
# TODO: change to Ebivariation-master once the eva-accesion-release PR is merged
# RUN git clone -b master https://github.com/EBIvariation/eva-accession.git
RUN git clone -b merge-deprecated-merged-jobs https://github.com/nitin-ebi/eva-accession.git
RUN mvn package -s /root/maven.xml -f eva-accession/pom.xml -DskipTests -P localhost
RUN cp eva-accession/eva-accession-release/target/eva-accession-release-*-exec.jar /usr/local/software/eva-accession-release.jar

COPY . /opt/

RUN cd /opt/ && python3 -m pip install .

WORKDIR /usr/local/test_eva_release/
