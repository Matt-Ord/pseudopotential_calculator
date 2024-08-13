#!/bin/bash

apt-get update

apt-get install -y gcc build-essential git python-is-python3 gfortran libopenblas-dev libfftw3-dev libopenmpi-dev python3-setuptools



CASTEP_VERSION="${VERSION:-"24.1"}"

CURRENT_DIR=$(dirname "${BASH_SOURCE[0]}")
cp "${CURRENT_DIR}/CASTEP-24.1.tar.gz" /tmp/CASTEP-${CASTEP_VERSION}.tar.gz

cd /tmp
tar xzf CASTEP-${CASTEP_VERSION}.tar.gz
cd CASTEP-${CASTEP_VERSION}
make -j 2

mkdir /opt/CASTEP/ 
cd /tmp/CASTEP-${CASTEP_VERSION} 
make install INSTALL_DIR="/opt/CASTEP/" 
cp /opt/CASTEP/castep.serial /usr/local/bin/
