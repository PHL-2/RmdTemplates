FROM amazoncorretto:17-al2023 as app

# List all software versions are ARGs near the top of the dockerfile
# 'ARG' sets environment variables during the build stage
# ARG variables are ONLY available during image build, they do not persist in the final image
ARG SOFTWARE_VER="24.10.3"

ENV NXF_VER="${SOFTWARE_VER}"
# set locale settings for singularity
ENV LC_ALL=C.UTF-8

# 'LABEL' instructions tag the image with metadata that might be important to the user
LABEL base.image="python:3.9.17-slim"
LABEL dockerfile.version="1"
LABEL software="nextflow"
LABEL software.version="${SOFTWARE_VER}"
LABEL description="Nextflow is a workflow manager that enables the deployment of complex, data-intensive pipelines in a portable and reproducible manner."
LABEL website="https://github.com/nextflow-io/nextflow"
LABEL license="https://github.com/nextflow-io/nextflow?tab=Apache-2.0-1-ov-file"
LABEL maintainer="Vincent Tu"
LABEL maintainer.email="vincent.tu@phila.gov"

# 'RUN' executes code during the build
# Install dependencies via apt-get or yum if using a centos or fedora base
# expect package provides the unbuffer program, which allows for commands to behave in a pty manner
RUN yum install -y \
  unzip \
  groff \
  ca-certificates \
  procps \
  expect \
  java-17-amazon-corretto

# Install and/or setup more things. Make /data for use as a working dir
# For readability, limit one install per 'RUN' statement.
RUN curl -s "https://get.nextflow.io" | bash && \
  chmod +x ./nextflow && \
  mv ./nextflow /usr/local/bin/ && \
  mkdir /data

# default command is to pull up help options
CMD nextflow help

# set working directory to /data
WORKDIR /data

# test layer
FROM app as test

# print version and grep for it; if not found, exit with error
RUN nextflow -version | grep "version ${SOFTWARE_VER}" || exit 1
