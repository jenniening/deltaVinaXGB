FROM ubuntu:16.04

RUN apt-get update && apt-get install --no-install-recommends -y g++ wget vim ca-certificates unzip build-essential git build-essential cmake libgtk-3-dev libboost-all-dev libcgal-dev apbs libcgal-demo scotch ptscotch \
                libeigen3-dev \
                libgfortran3 \
                libgmp-dev \
                libgmpxx4ldbl \
                libmpfr-dev \
                libboost-dev \
                libboost-thread-dev \
                libtbb-dev \
                python3-dev && rm -rf /var/lib/apt/lists/*

#Install MINICONDA
RUN wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda.sh && \
	/bin/bash Miniconda.sh -b -p /opt/conda && \
	rm Miniconda.sh
ENV PATH /opt/conda/bin:$PATH

#COPY Makefile requirements.txt setup.py test_environment.py  /app/
COPY . /app
WORKDIR /app
RUN make Makefile create_environment
SHELL ["conda", "run", "-n", "DXGB", "/bin/bash", "-c"]
RUN /bin/bash -c "source activate DXGB" && make Makefile requirements
RUN conda install -c conda-forge xgboost=0.80.0 && conda install -c rdkit rdkit=2019.03.1 && conda install -c openbabel openbabel && conda install moleculekit -c acellera && python setup.py install
RUN wget http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz
RUN tar -xvzf mgltools_x86_64Linux2_1.5.6.tar.gz && cd mgltools_x86_64Linux2_1.5.6/ && /bin/bash -c "source install.sh"
WORKDIR /app
RUN wget http://mgltools.scripps.edu/downloads/tars/releases/MSMSRELEASE/REL2.6.1/msms_i86_64Linux2_2.6.1.tar.gz
RUN mkdir msms && tar -xvzf msms_i86_64Linux2_2.6.1.tar.gz -C msms && cd msms && cp msms.x86_64Linux2.2.6.1 msms
WORKDIR /app/msms
RUN sed 's+./atmtypenumbers+/app/DXGB/atmtypenumbers+g' pdb_to_xyzr | tee pdb_to_xyzr
WORKDIR /app
RUN wget https://github.com/chengwang88/vina4dv/archive/master.zip && unzip master.zip

ENV PATH /opt/conda/envs/DXGB/bin:$PATH
RUN /bin/bash -c "source activate DXGB"
# path for MSMS
ENV PATH=$PATH:/app/msms/
# set mgltool variable (if mac, should change mgltools_x86_64Linux2_1.5.6 into your downloaded mac version)
ENV PATH=$PATH:/app/mgltools_x86_64Linux2_1.5.6/bin/
ENV MGL=/app/mgltools_x86_64Linux2_1.5.6/
ENV MGLPY=$MGL/bin/python
ENV MGLUTIL=$MGL/MGLToolsPckgs/AutoDockTools/Utilities24/
# set vina dir  (if mac, should use /mac/release/)
ENV VINADIR=/app/vina4dv-master/build/linux/release
# set DXGB dir (needed if run deltaVinaRF20 score)
ENV DXGB=/app/DXGB




