# using debian-based miniconda image
FROM continuumio/miniconda3:latest

#install commands
RUN apt-get update && apt-get install -y \
sudo \
wget \
vim \
git \
nginx \
unzip

#install commands(for analysis)
RUN conda install -c bioconda -y \
pandas \
numpy \
matplotlib \
tqdm \
fastp \
blast \
skesa && \
python3 -m pip install --user --upgrade cutadapt

#create working folder
RUN mkdir /ngs


#open port
EXPOSE 80

#set alias
RUN echo 'alias ll="ls -l"' >> ~/.bashrc
RUN echo 'alias la="ls -a"' >> ~/.bashrc

#entry point
#activate nginx and execute whatever the container recieved
RUN echo "#!/bin/bash" >> /entry.sh
RUN echo "nginx" >> /entry.sh
RUN echo "exec \$@" >> /entry.sh
RUN chmod 755 /entry.sh
ENTRYPOINT [ "/entry.sh" ]