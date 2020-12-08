# using debian-based miniconda image
FROM continuumio/miniconda3:latest

#install commands
RUN apt-get update && apt-get install -y \
sudo \
wget \
vim \
git \
nginx \
unzip \
libyaml-dev &&\
apt-get clean && \
rm -rf /var/lib/apt/lists/*

#install commands(for analysis)
RUN conda install -c bioconda -y \
pandas \
numpy \
matplotlib \
tqdm \
fastp \
blast \
megahit \
spades \
skesa && \
python3 -m pip install --user --upgrade cutadapt ruamel.yaml && \
conda clean -a

#create working folder
RUN mkdir /ngs

#open port
EXPOSE 80

#set alias and create entry point
#activate nginx and execute whatever the container recieved
RUN echo 'alias ll="ls -l"' >> ~/.bashrc && \
echo 'alias la="ls -a"' >> ~/.bashrc &&\
echo "#!/bin/bash" >> /entry.sh && \
echo "nginx" >> /entry.sh && \
echo "exec \$@" >> /entry.sh && \
chmod 755 /entry.sh
ENTRYPOINT [ "/entry.sh" ]