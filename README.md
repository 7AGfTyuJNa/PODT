# PODTE

This codebase is part of the PODTE paper.

**What this codebase includes**: example and benchmark implementations in C++14 for some of the schemes in the PODTE paper.

**What this codebase is not**: it is not for production use; it is not extensively tested.

## Setup and Building Instructions

First, install the libraries [NTL](https://libntl.org/doc/tour-unix.html) required by PODTE. On several Ubuntu systems this can be done directly through links above.

To build library and executables:
```shell
mkdir build
cd build
cmake ..
make
```

To try an example, run e.g.:
```shell
build/PODT-code
```
