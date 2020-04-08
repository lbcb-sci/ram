# Ram

[![Latest GitHub release](https://img.shields.io/github/release/lbcb-sci/ram.svg)](https://github.com/lbcb-sci/ram/releases/latest)
[![Build status for c++/clang++](https://travis-ci.org/lbcb-sci/ram.svg?branch=master)](https://travis-ci.org/lbcb-sci/ram)

Ram is a c++ implementation of [minimap](https://github.com/lh3/minimap) with few modifications.

# Usage

To build ram run the following commands:
```bash
git clone --recursive https://github.com/lbcb-sci/ram.git ram
cd ram && mkdir build && cd build
cmake -Dram_build_executable=ON -DCMAKE_BUILD_TYPE=Release .. && make
./bin/ram
```
which will display the following usage:
```bash
usage: ram [options ...] <target> [<sequences>]

  # default output is stdout
  <target>/<sequences>
    input file in FASTA/FASTQ format (can be compressed with gzip)

  options:
    -k, --kmer-length <int>
      default: 15
      length of minimizers
    -w, --window-length <int>
      default: 5
      length of sliding window from which minimizers are found
    -f, --frequency-threshold <float>
      default: 0.001
      threshold for ignoring most frequent minimizers
    -m, --micromize
      use only a portion of all minimizers
    -t, --threads <int>
      default: 1
      number of threads
    --version
      prints the version number
    -h, --help
      prints the usage
```

If you would like to add ram as a library to your project via CMake, add the following:
```cmake
if (NOT TARGET ram)
  add_subdirectory(<path_to_submodules>/ram EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(<your_exe> ram)
```

#### Dependencies

- gcc 4.8+ or clang 3.5+
- cmake 3.9+
- zlib (for binary only)

## Unit tests

To build ram unit tests run the following commands:

```bash
git clone https://github.com/lbcb-sci/ram.git ram
cd ram && mkdir build && cd build
cmake -Dram_build_tests=ON -DCMAKE_BUILD_TYPE=Release .. && make
./bin/ram_test
```

#### Dependencies
- gtest

## Acknowledgement

This work has been supported in part by the European Regional Development Fund under the grant KK.01.1.1.01.0009 (DATACROSS) and in part by the Croatian Science Foundation under the project Single genome and metagenome assembly (IP-2018-01-5886).
