# Unsteady-MODES
Adaptive Reduced Order Modelling for unsteady aerodynamics


NOTE: currently under heavy development.

Requirements
------

To install this project, please ensure that you have installed the following (installation guides are provided on the respective websites or github repositories):

SOFTWARES:
  - [Git](http://git-scm.com)
  - A C++ compiler, e.g., [GCC](https://gcc.gnu.org/), [clang](http://clang.llvm.org/) [Verified on GCC v6.4.0 or newer]
  - [CMake](http://www.cmake.org)

LIBRARIES:
  - Optional
         - [pagmo](https://esa.github.io/pagmo2/install.html)
         - [GSL](https://www.gnu.org/software/gsl/)

Compilation
------

Run the following commands to download, build, and compile this project.

    cd path/to/compilation/
    git clone https://github.com/strath-ace-labs/MODES.git
    cd Unsteady-MODES
    mkdir build && cd build
    cmake ..
    make


Build options
-------------

To be able to use optional capabilities you can change the options in the CMakelists.txt in MODES main folder before running cmake or run `ccmake ..` to bring up the interface that can be used to toggle options.


Running the code
-------------
For instruction on how to run the code, please refer to the tutorial folder


