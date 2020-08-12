FROM ubuntu:20.04

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils build-essential sudo git wget python3 default-jre pigz unzip python3-pip samtools && apt-get clean

RUN mkdir /bio-tools
WORKDIR /bio-tools
RUN wget -O subread.tar.gz "https://sourceforge.net/projects/subread/files/subread-2.0.1/subread-2.0.1-Linux-x86_64.tar.gz/download" && tar -xzvf subread.tar.gz && cp subread-2.0.1-Linux-x86_64/bin/featureCounts . && rm -R subread-2.0.1-Linux-x86_64 subread.tar.gz
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && unzip fastqc_v0.11.9.zip && chmod +x FastQC/fastqc && sudo ln -sf /bio-tools/FastQC/fastqc /usr/local/bin/fastqc && rm fastqc_v0.11.9.zip
RUN wget https://github.com/broadinstitute/picard/releases/download/2.23.3/picard.jar
RUN wget -O HISAT2.zip "https://cloud.biohpc.swmed.edu/index.php/s/4pMgDq4oAF9QCfA/download" && unzip HISAT2.zip && rm HISAT2.zip
RUN wget -O BBmap.tar.gz "https://sourceforge.net/projects/bbmap/files/latest/download" && tar -xzvf BBmap.tar.gz && rm BBmap.tar.gz
RUN python3 -m pip install --system --upgrade multiqc
RUN wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip && unzip qualimap_v2.2.1.zip && rm qualimap_v2.2.1.zip
COPY config.yaml /bio-tools/config.yaml
ENV PATH="/bio-tools:${PATH}" 
ENV PATH="/bio-tools/hisat2-2.2.1/:${PATH}"
ENV PATH="/bio-tools/qualimap_v2.2.1/:${PATH}"
RUN ln -s /usr/bin/python3 /usr/bin/python

RUN mkdir /src
WORKDIR /src
