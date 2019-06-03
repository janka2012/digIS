FROM ubuntu:18.04

RUN apt-get update && \
    apt-get install -y hmmer=3.1b2+dfsg-5ubuntu1 \
    ncbi-blast+=2.6.0-1 \
    python3.7 \
    python3-pip

RUN pip3 install biopython==1.73

ENV INSTALL_PATH /digis_tool
RUN mkdir -p ${INSTALL_PATH}

WORKDIR ${INSTALL_PATH}

COPY . ${INSTALL_PATH}

# Run digIS_search.py when the container launches
CMD ["python3", "digIS_search.py"]
