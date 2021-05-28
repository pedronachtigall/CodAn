FROM ubuntu:18.04

MAINTAINER Pedro G Nachtigall <pedronachtigall[AT]gmail[DOT]com>

ENV DEBIAN_FRONTEND=noninteractive

# install required packages
RUN apt-get update
RUN apt-get install -y git python3 python3-biopython perl bioperl libmce-perl ncbi-blast+ unzip

# make an app folder
RUN mkdir /app

# download CodAn and install it
RUN git clone https://github.com/pedronachtigall/CodAn.git && mv CodAn /app && chmod +x /app/CodAn/bin/*

# decompress the models
RUN cd /app/CodAn/models && (for i in *.zip; do (echo $i; unzip $i); done;)

ENV LC_ALL=C
ENV PATH=/app/CodAn/bin:$PATH

VOLUME /project
WORKDIR /project

CMD ["/bin/bash"]
