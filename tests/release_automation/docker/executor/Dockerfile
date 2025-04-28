FROM ubuntu:22.04

ENV JAVA_VERSION=8
ENV NXF_VER=21.10.5
ENV EVA_ACCESSION_VERSION=v0.6.21

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

# setup eva-accession
#RUN git clone -b ${EVA_ACCESSION_VERSION} https://github.com/EBIvariation/eva-accession.git
#RUN mvn package -s /root/maven.xml -f eva-accession/pom.xml -DskipTests -P localhost
#RUN cp eva-accession/eva-accession-pipeline/target/eva-accession-release-*-exec*.jar /usr/local/software/eva-accession-release.jar
RUN curl -O -J https://dl.cloudsmith.io/public/ebivariation/packages/maven/uk/ac/ebi/eva/eva-accession-release/0.6.47-SNAPSHOT/eva-accession-release-0.6.47-20250422.150441-1-exec.jar
RUN cp eva-accession-release-*exec.jar /usr/local/software/eva-accession-release.jar


COPY tests/release_automation/docker/executor/maven.xml /root
COPY tests/release_automation/docker/executor/release_config.yml /root/.release_config.yml
COPY tests/release_automation/docker/executor/nextflow.config /root/.nextflow/config
COPY . /opt/

RUN cd /opt/ && ls && python3 -m pip install .

WORKDIR /usr/local/test_eva_release/
