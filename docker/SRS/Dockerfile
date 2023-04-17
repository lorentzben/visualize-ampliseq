FROM quay.io/qiime2/core:2020.8
LABEL maintainer="lorentzb"
RUN apt-get update --allow-releaseinfo-change -y --fix-missing && apt-get install -f -y  git procps nano && rm -rf /var/lib/apt/lists/*
#RUN conda install -c vitorheidrich q2_srs
RUN pip install git+https://github.com/vitorheidrich/q2-srs.git
CMD [ "bash" ]