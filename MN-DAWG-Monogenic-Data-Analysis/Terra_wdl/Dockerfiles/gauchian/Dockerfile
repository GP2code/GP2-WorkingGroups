FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -q && \
    apt-get install -q -y --no-install-recommends \
    bzip2 \
    ca-certificates \
    git \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender1 \
    mercurial \
    openssh-client \
    procps \
    subversion \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

#install miniconda
RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

#install Gauchain
ENV Gauchain_VERSION="1.0.2"
RUN wget https://github.com/Illumina/Gauchian/archive/refs/tags/v${Gauchain_VERSION}.tar.gz \
        && tar xvzf v${Gauchain_VERSION}.tar.gz \
        && rm v${Gauchain_VERSION}.tar.gz \
	&& cd Gauchian-${Gauchain_VERSION} \
	&& python setup.py install

