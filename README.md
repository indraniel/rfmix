# RFMix _(v1.5.4)_


This is a `pip`-installable version of [RFMix](https://sites.google.com/site/rfmixlocalancestryinference/), Local Ancestry and Admixture Inference, version 1.  You also might be interested in [RFMix v2](https://github.com/slowkoni/rfmix) instead.

## Installation

RFMix assumes that you have the [`openmp` support](https://www.openmp.org/resources/openmp-compilers-tools/) enabled on your C++ compiler.

### MacOS

On Mac OS, you'll need `g++`, which can be installed via [homebrew](https://brew.sh) as:

    brew install gcc

Afterwards, you can do:

    CXX=/usr/local/bin/g++-9 pip install --no-cache-dir git+https://github.com/indraniel/rfmix.git

### Linux

#### Debian/Ubuntu based distributions

    sudo apt-get install libomp-dev
    pip install --no-cache-dir git+https://github.com/indraniel/rfmix.git
