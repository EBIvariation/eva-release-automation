FROM ubuntu:22.04

ENV JAVA_VERSION=8
ENV NXF_VER=21.10.5
ENV VCF_VALIDATOR_VERSION=0.9.4
ENV CONTIG_NAME_CONVERTER_VERSION=v0.0.2

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

# Download release pipeline
RUN curl -O -J https://dl.cloudsmith.io/public/ebivariation/packages/maven/uk/ac/ebi/eva/eva-accession-release/0.6.47-SNAPSHOT/eva-accession-release-0.6.47-20250422.150441-1-exec.jar
RUN cp eva-accession-release-*exec.jar /usr/local/software/eva-accession-release.jar


# install vcf_validator
RUN curl -LJo /usr/local/bin/vcf_validator https://github.com/EBIvariation/vcf-validator/releases/download/v${VCF_VALIDATOR_VERSION}/vcf_validator_linux \
    && chmod 755 /usr/local/bin/vcf_validator

# install vcf_assembly_checker
RUN curl -LJo /usr/local/bin/vcf_assembly_checker https://github.com/EBIvariation/vcf-validator/releases/download/v${VCF_VALIDATOR_VERSION}/vcf_assembly_checker_linux \
    && chmod 755 /usr/local/bin/vcf_assembly_checker

# Install contig_name_converter
#RUN pip install -q "git+https://github.com/EBIvariation/contig_name_converter.git@${CONTIG_NAME_CONVERTER_VERSION}"
RUN pip install -q "git+https://github.com/tcezard/contig_name_converter.git@EVA3793-allow-accession-as-id"

COPY tests/release_automation/docker/executor/maven.xml /root
COPY tests/release_automation/docker/executor/release_config.yml /root/.release_config.yml
COPY tests/release_automation/docker/executor/nextflow.config /root/.nextflow/config
COPY . /opt/

RUN cd /opt/ && python3 -m pip install .

WORKDIR /usr/local/test_eva_release/
