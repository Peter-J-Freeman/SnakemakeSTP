# Use the Miniconda base image
FROM continuumio/miniconda3:latest

# Update the package list and install necessary packages using apt-get
# Note: I'm running a mac M1 and these packages do not seem to be available
# for ARM via conda, so using apt-get to install
RUN apt-get update && apt-get install -y \
    bwa \
    bcftools \
    samtools \
    wget

# Set the working directory
WORKDIR /snakemake_project

# Copy the environment.yml file to the container
COPY environment.yml .

# Create the Conda environment
RUN conda env create -f environment.yml

# Updrade pip
RUN pip install --upgrade pip

# Install any needed packages specified in requirements.txt
RUN pip install -r requirements.txt

# Activate the Conda environment
RUN echo "conda activate $(head -1 environment.yml | cut -d' ' -f2)" > ~/.bashrc

# Set the entry point to bash
CMD ["bash"]

