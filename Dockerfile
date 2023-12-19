FROM ncbi/blast
#PREAMBLE

ENV BAMTOOLS_VERSION 2.4.1
ENV BAMTOOLS_NAME bamtools
ENV BAMTOOLS_URL "https://github.com/pezmaster31/bamtools/archive/v${BAMTOOLS_VERSION}.tar.gz"

WORKDIR /home/genomics
COPY . /home/genomics
RUN cd /home/genomics

RUN apt-get --assume-yes update \
    && apt-get --assume-yes upgrade \
    && apt-get --assume-yes install wget cmake zlib1g-dev unzip

RUN wget -q -O - $BAMTOOLS_URL | tar -zxv \
    && cd ${BAMTOOLS_NAME}-${BAMTOOLS_VERSION} \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && cd ../.. \
    && cp ./${BAMTOOLS_NAME}-${BAMTOOLS_VERSION}/lib/libbamtools.so.${BAMTOOLS_VERSION} /usr/lib/ \
    && rm -rf ./${BAMTOOLS_NAME}-${BAMTOOLS_VERSION} \
    && strip /usr/local/bin/*; true

RUN wget https://github.com/sheinasim/HiFiAdapterFilt/archive/refs/heads/master.zip \
    && unzip master.zip \
    && cp HiFiAdapterFilt-master/pbadapterfilt.sh /usr/local/bin \
    && cp -R HiFiAdapterFilt-master/DB /usr/local/bin

RUN rm -rf *.tgz *.tar *.zip \
    && rm -rf /var/cache/apk/* \
    && rm -rf /tmp/*
