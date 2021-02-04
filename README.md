# Ram

[![Latest GitHub release](https://img.shields.io/github/release/lbcb-sci/ram.svg)](https://github.com/lbcb-sci/ram/releases/latest)
[![Build status for c++/clang++](https://travis-ci.com/lbcb-sci/ram.svg?branch=master)](https://travis-ci.com/lbcb-sci/ram)

Ram is a c++ implementation of [minimap](https://github.com/lh3/minimap) with few modifications.

# Usage

To build ram run the following commands:

```bash
git clone https://github.com/lbcb-sci/ram && cd ram && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make
```

which will create ram library, executable and unit tests. Running the executable will display the following usage:

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
      length of sliding window from which minimizers are sampled
    -f, --frequency-threshold <float>
      default: 0.001
      threshold for ignoring most frequent minimizers
    --bandwidth <int>
      default: 500
      size of bandwidth in which minimizer hits can be chained
    --chain <int>
      default: 4
      minimal number of chained minimizer hits in overlap
    --matches <int>
      default: 100
      minimal number of matching bases in overlap
    --gap <int>
      default: 10000
      maximal gap between minimizer hits in a chain
    --minhash
      use only a portion of all minimizers
    -t, --threads <int>
      default: 1
      number of threads
    --version
      prints the version number
    -h, --help
      prints the usage
```

Running `make install` will install the executable. In order to install the library, both biosoup and thread_pool (see Dependencies) need to be installed beforehand, and option `ram_install` used while configuring the build. Once the library is installed, a package will be copied to your system that can be searched and linked with:

```cmake
find_package(ram)
target_link_libraries(<target> ram::ram)
```

On the other hand, you can include ram as a submodule and add it to your project with the following:

```cmake
if (NOT TARGET ram)
  add_subdirectory(<path_to_submodules>/ram EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(<target> ram)
```

#### Build options

- `ram_install`: generate library install target
- `ram_build_tests`: build unit tests

#### Dependencies

- gcc 4.8+ | clang 3.5+
- cmake 3.11+
- pthread
- (ram_exe) (ram_test) zlib 1.2.8+

###### Hidden

- biosoup 0.10.0
- thread_pool 3.0.3
- (ram_exe) (ram_test) bioparser 3.0.13
- (ram_test) googletest 1.10.0

## Acknowledgement

This work has been supported in part by the European Regional Development Fund under the grant KK.01.1.1.01.0009 (DATACROSS) and in part by the Croatian Science Foundation under the project Single genome and metagenome assembly (IP-2018-01-5886).
