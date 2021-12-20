# ParGeo: A Library for Parallal Algorithms in Computational Geometry

## Tutorial

ParGeo uses the [parlaylib](https://github.com/cmuparlay/parlaylib) as a submodule for multi-core parallel programming. Before building the project, initialize the submodules:

```
git submodule init
git submodule update
```

ParGeo uses CMake as the build system. To build the project:

```
mkdir build
cd build
cmake ..
make -j // If you don't want to build everything, you could cd to the desired build/folder and then make
```

After the build completes, you will be able to find the executables in the nested directories of `build/`. Running `./program` in the terminal will return short usage instructions. A typical run will involve `./program <options> <data-path>`. For example, to run the HDBSCAN clustering

```
./hdbscan -m 3 data.txt
```

Example data sets can be found [here](https://github.com/wangyiqiu/pargeo/tree/main/test/datasets). Pargeo automatically parses csv-like point data files with or without header. The delimiters (space, comma etc) have to be one character each in length.
