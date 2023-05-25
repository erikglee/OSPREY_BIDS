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
                    zlib1g \
                    pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/

#Install relavent python packages
#RUN apt update 
#RUN apt install software-properties-common -y
#RUN add-apt-repository ppa:deadsnakes/ppa -y
#RUN apt update 
#ENV DEBIAN_FRONTEND=noninteractive
#RUN apt install python3.9 -y
#RUN apt-get install -y python3
#RUN apt-get install -y python3-dev
#RUN python3 -m pip install python-dev-tools
RUN python3 -m pip install numpy==1.19.2
RUN python3 -m pip install scipy==1.8.0
RUN python3 -m pip install nibabel==3.2.2
RUN python3 -m pip install matplotlib==3.5.1


#Setup MCR - this grabs v910 of MCR that was downloaded from the matlab
#website, installed at MSI, and then zipped. If you want to use a
#different version of matlab then download the corresponding version
#of MCR, install it, zip it, and upload the new path to a public bucket
#on S3
RUN mkdir /mcr_path
RUN wget https://s3.msi.umn.edu/leex6144-public/v910.zip -O /mcr_path/mcr.zip
RUN cd /mcr_path && unzip -q ./mcr.zip
RUN rm /mcr_path/mcr.zip

#Download the unique code for this project
RUN mkdir /code
RUN wget https://s3.msi.umn.edu/leex6144-public/osprey_v1.2.zip -O /code/code.zip
RUN cd /code && unzip -q ./code.zip
RUN mkdir /python_code
COPY ./code/run.py /python_code
COPY ./code/localizer_alignment.py /python_code 
COPY hbcd_pilot_config.json /python_code
RUN rm /code/code.zip

#Export paths
ENV MCR_PATH=/mcr_path
ENV EXECUTABLE_PATH=/code/HBCD/run_compiled.sh

#Set permissions
RUN chmod 555 -R /mcr_path /code /python_code

#Add code dir to path
ENV PATH="${PATH}:/python_code"
RUN pipeline_name=osprey && cp /python_code/run.py /python_code/$pipeline_name

ENTRYPOINT ["run.py"]
