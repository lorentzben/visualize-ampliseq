FROM quay.io/qiime2/core:2020.8
LABEL maintainer="lorentzb"
RUN apt-get update --allow-releaseinfo-change -y --fix-missing && apt-get install -f -y  git procps nano && rm -rf /var/lib/apt/lists/*
#RUN conda install -c vitorheidrich q2_srs
RUN pip install git+https://github.com/vitorheidrich/q2-srs.git
#RUN QIIME2_RELEASE=2020.8 /bin/sh -c qiime dev refresh-cache
#RUN QIIME2_RELEASE=2020.8 /bin/sh -c /bin/bash -c "source activate qiime2-${QIIME2_RELEASE}"
#CMD [ "qiime info" ]
CMD [ "bash" ]