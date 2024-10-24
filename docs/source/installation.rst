.. OSPREY_BIDS documentation master file, created by
   sphinx-quickstart on Wed Jun  5 10:48:12 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation
============

The intended use of this pipeline is through the use of a `Singularity <https://docs.sylabs.io/guides/3.7/user-guide/index.html>`_ or `Docker <https://docs.docker.com/get-started/>`_
image. The image can be built using the Dockerfile found in the `repository <https://github.com/erikglee/OSPREY_BIDS>`_,
or it can be pulled from `Dockerhub <https://hub.docker.com/r/dcanumn/osprey-bids/tags>`_ as a singularity using the following command: ::
    
        singularity pull docker://dcanumn/osprey:<version_num>

Where version_num denotes the specific version of the container. All available
versions of the container can be found `here <https://hub.docker.com/r/dcanumn/osprey-bids/tags>`_.

After downloading the container, singularity is the only other dependency needed
for processing. The full usage details can be seen under the :doc:`usage section <usage>`, but
the basic command to run the container is as follows: ::
    
        container_path=/path/to/container.sif
        bids_dir=/path/to/bids
        output_dir=/path/to/output
        bibsnet_dir=/path/to/bibsnet
        settings_file=/path/to/settings.json
        singularity run -B $bids_dir:/bids \
         -B $output_dir:/output \
         -B $bibsnet_dir:/bibsnet \
         -B $settings_file:/settings_file/file.json \
         $container_path /bids /output participant /settings_file/file.json

Where "singularity run" is followed by specific commands for singularity.
In this case it is a series of "bind" commands that will give singularity
access to the necessary directories for processing. Of note, you will need
to be sure that singularity is given access to the **BIDS directory**, the **output
directory** where you want pipeline outputs to be stored, and the **JSON file**
that will configure how processing is run.

The singularity specific arguments are then followed by the path to the
container. Finally, the user specifies the arguments that are unique to the current application,
such as input and output files, configuration files, and other processing settings.

Note: If you are not interested in interacting with this particular interface
of OSPREY, which includes BIDS and containerized functionality, consider using
OSPREY directly in Matlab following the instructions found `on the OSPREY github <https://github.com/schorschinho/osprey>`_.