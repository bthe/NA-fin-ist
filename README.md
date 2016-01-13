# NA-fin-ist
Informal code repository for North Atlantic fin whale simulation trials

# Setup 

To compile the naf-ist program one needs a working fortran compiler. These can be obtained from various source both commercial vendorsor as a part of the gnu compiler collection. The gnu fortran compiler (gfortran) is available for the common computer 
platforms. 

## Linux

Simply install the gnu compiler collection using the default package manager (yum/apt-get etc..)

## MacOsX

Follow the installation instructions on http://hpc.sourceforge.net/

## Windows

The gnu compiler collection can be downloaded as a part of the minimalistic GNU for Windows, see http://www.mingw.org for further details or alternatively rtools.

## Compiling

Assuming you have gfortran installed one simply has to navigate to the git directory (via terminal) and type:

> make

and to intall (assuming you are on linux/unix platforms):

> mv naf-ist ~/bin

# Running 

The naf-ist program can be run directly from the command line:

> naf-ist -main copyna 

