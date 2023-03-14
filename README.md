# bumo: BUilding MOrphology 

CGAL-based reimplementation of https://github.com/tudelft3d/3d-building-metrics

It should be faster and robuster (to faulty input) and easier to install.

Most of the metrics described in the paper https://doi.org/10.7910/DVN/6QCRRF

- reads only CityJSON v1.1 files
- only Solid are processed
- made more-or-less for the <3dbag.nl>, but should work with any file

## Compilation

You first need to install the following free libraries:

  1. [CGAL v5.0+](http://www.cgal.org) 
  1. [Eigen library](http://eigen.tuxfamily.org)
  1. [CMake](http://www.cmake.org)

Under macOS, it's super easy, we suggest using [Homebrew](http://brew.sh/):

    $ brew install cgal
    $ brew install eigen
    $ brew install cmake

and then

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make


## Usage

  ```bash
  ./bumo myfile.city.json > metrics.csv
  ```