#The base image is the AMD64 version of centos:centos7.9.2009, which
#should correspond to the OS at MSI
#FROM amd64/centos:7.9.2009
#FROM ubuntu:latest
FROM python:3.9.16-slim-bullseye

# Prepare environment
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    apt-utils \
                    autoconf \
                    build-essential \
                    bzip2 \
                    ca-certificates \
                    curl \
                    gcc \
                    git \
                    gnupg \
                    libtool \
                    lsb-release \
                    pkg-config \
                    unzip \
                    wget \
                    xvfb \
                    default-jre \
                    zlib1g \
                    pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/

#Install relavent python packages
RUN python3 -m pip install --upgrade pip 
RUN python3 -m pip install numpy==1.19.5 #used to be 1.19.2
RUN python3 -m pip install scipy==1.8.0
RUN python3 -m pip install nibabel==3.2.2
RUN python3 -m pip install matplotlib==3.5.1


#Setup MCR - this grabs v910 of MCR that was downloaded from the matlab
#website, installed at MSI, and then zipped. If you want to use a
#different version of matlab then download the corresponding version
#of MCR, install it, zip it, and upload the new path to a public bucket
#on S3
RUN mkdir /mcr_path
RUN wget https://s3.msi.umn.edu/leex6144-public/v912.zip -O /mcr_path/mcr.zip \
    && cd /mcr_path && unzip -q ./mcr.zip \
    && rm /mcr_path/mcr.zip

#Download the unique code for this project
RUN mkdir /code
RUN wget https://s3.msi.umn.edu/leex6144-public/osprey_v2.4.0.zip -O /code/code.zip \
    && cd /code \
    && unzip -q ./code.zip \
    && rm /code/code.zip
RUN mkdir /python_code
COPY ./python_code/run.py /python_code
COPY ./python_code/localizer_alignment.py /python_code 
COPY hbcd_pilot_config.json /python_code

#Download the basis sets
RUN mkdir /HBCD_basissets
RUN wget https://s3.msi.umn.edu/leex6144-public/OSPREY_HBCD_BASISSETS.zip -O /HBCD_basissets/OSPREY_HBCD_BASISSETS.zip \
    && cd /HBCD_basissets \
    && unzip -q ./OSPREY_HBCD_BASISSETS.zip \
    && rm /HBCD_basissets/OSPREY_HBCD_BASISSETS.zip

#Export paths (make sure LD_LIBRARY_PATH is set to the correct version)
ENV BASIS_SETS_PATH=/HBCD_basissets
ENV MCR_PATH=/mcr_path/v912
ENV EXECUTABLE_PATH=/code/run_compiled.sh
ENV LD_LIBRARY_PATH ="${LD_LIBRARY_PATH}:/mcr_path/v912/runtime/glnxa64:/mcr_path/v912/bin/glnxa64:/mcr_path/v912/sys/os/glnxa64:/mcr_path/v912/extern/bin/glnxa64"

#Set permissions
RUN chmod 555 -R /mcr_path /code /python_code /HBCD_basissets

#Add code dir to path
ENV PATH="${PATH}:/python_code"
RUN pipeline_name=osprey && cp /python_code/run.py /python_code/$pipeline_name

ENTRYPOINT ["osprey"]
