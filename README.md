# Ram

[![Build status for c++/clang++](https://travis-ci.org/jmaricb/RAM.svg?branch=master)](https://travis-ci.org/jmaricb/RAM)

Mapping module for raw de novo genome assembly of long uncorrected reads.

## Description

Ram is a [Minimap](https://github.com/lh3/minimap) clone implemented as a C++ library.

## Dependencies

1. gcc 4.8+ or clang 3.4+
2. cmake 3.2+

## Installation

CmakeLists is provided in the project root folder. By running the following commands:

```bash
git clone --recursive https://github.com/jmaricb/ram ram
cd ram
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
a library named `libram.a` will appear in the `build/lib` directory. If you want the Ram executable, run the following two commands:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -Dram_build_executable=ON ..
make
```
which will place an executable named `ram` in `build/bin` directory.

Optionally, you can run `sudo make install` to install Ram library (and executable) to your machine.

***Note***: if you omitted `--recursive` from `git clone`, run `git submodule init` and `git submodule update` before proceeding with compilation.

## Usage

Usage of ram is as following:

```bash
usage: ram [options ...] <sequences> [<sequences> ...]

    <sequences>
        input file in FASTA/FASTQ format (can be compressed with gzip)
        containing sequences

    options:
        -k, --kmer-length <int>
            default: 15
            length of minimizers
        -w, --window-length <int>
            default: 5
            window length from which minimizers are found
        -f, --filter-threshold <float>
            default: 0.001
            threshold for ignoring most frequent minimizers
        -t, --threads <int>
            default: 1
            number of threads
        --version
            prints the version number
        -h, --help
            prints the usage
```

## Contact information

For additional information, help and bug reports please send an email to one of the following: josip.maric@fer.hr, robert.vaser@fer.hr.

## Acknowledgement

This work has been supported in part by the European Regional Development Fund under the grant KK.01.1.1.01.0009 (DATACROSS).
