# bumo: BUilding MOrphology 

A C++ and purely CGAL-based reimplementation of the 3D building metrics described in the article:

> A. Labetski, S. Vitalis, F. Biljecki, K. A. Ohori, and J. Stoter. 3D building metrics for urban morphology. International Journal of Geographical Information Science, pages 1–32, 2022. [[DOI]](https://doi.org/10.1080/13658816.2022.2103818)

The original code is written in Python and relies mostly on [PyVista](https://pyvista.org): https://github.com/tudelft3d/3d-building-metrics


This implementation should be faster, robuster to faulty input, and easier to install and to maintain.

## Good to know

  - reads only CityJSON v1.1 files
  - only Solid are processed
  - made more-or-less for the [3dbag.nl](https://3dbag.nl), but should work with any file


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