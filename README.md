---
output: html_document
---
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

The gnu compiler collection can be downloaded as a part of the minimalistic GNU for Windows, see http://www.mingw.org for further details or alternatively rtools (https://cran.r-project.org/bin/windows/Rtools/).

## Compiling

Assuming you have gfortran installed one simply has to navigate to the git directory (via terminal) and type:

> make

and to intall (assuming you are on linux/unix platforms):

> mv naf-ist ~/bin

# Running 

The naf-ist program can be run directly from the command line:

> naf-ist -main copyna 


# List of trials


| Trial No.    | Stock Hypothesis | MSYR 2 | No. of Stocks | Trial description                | Weight 1% | Weight 4% |
| ------------ | ---------------- | ------ | ------------- | -------------------------------- | --------- | --------- |
| **Baseline** |                  |        |               |                                  |           |           |
| NF-B1        | I                | 1, 4%  | 4             | 4 stocks, separate feeding areas |     M     |      H    |
| NF-B2        | II               | 1, 4%  | 4             | 4 stocks;  'W' & 'E' feed in central sub-areas | M | H |
| NF-B3 | III | 1, 4% | 4 | 4 stocks; 'C1' & 'C3' feed in adjacent sub-areas | M | H |
| NF-B5 | V | 1, 4% | 4 | 4 stocks as in Hypothesis I but stock 'S' in adjacent sub-areas | M | H |
| NF-B6 | VI | 4% | 3 | 3 stocks  (no 'E' stock) | n/a | H |
| **Sensitivity** |
| NF-H2 | II | 1, 4% | 4 | High historical catch series | M | M |
| NF-H3 | III | 1, 4% | 4 | High historical catch series | M | M |
| NF-Q3\* | III | 1, 4% | 4 | Future WI & EI/F surveys exc. strata S 60ºN | M | M |
| NF-A2 | II | 1, 4% | 4 | Pro-rate abundance data for conditioning | M | M |
| NF-A3 | III | 1, 4% | 4 | Pro-rate abundance data for conditioning | M | M |
| NF-U3 | III | 1, 4% | 4 | Selectivity decreases by 4%/yr for age 8+; _M_=0.04 | M | M |
| NF-G2 | II | 1, 4% | 4 | C2 sub-stock enters EG beginning yr 1985 (opt. a) | M | M |
| NF-G3 | III | 1, 4% | 4 | C2 sub-stock enters EG beginning yr 1985 (opt. a) | M | M |
| NF-F2 | II | 1, 4% | 4 | C2 sub-stock enters EG 1985-2025 (opt. b) | M | M |
| NF-F3 | III | 1, 4% | 4 | C2 sub-stock enters EG 1985-2025 (opt. b) | M | M |
| NF-S3 | III | 1, 4% | 4 | Selectivity estimated for pre and post 2000 & use all age data | M | M |
| NF-Y1\* | I | 1, 4% | 4 | 8 year future survey interval | M | H |
| NF-Y1\* | II | 1, 4% | 4 | 8 year future survey interval | M | H |
| NF-Y3\* | III | 1, 4% | 4 | 8 year future survey interval | M | H |
| NF-Y5\* | V | 1, 4% | 4 | 8 year future survey interval | M | H |
| NF-Y6\* | VI | 1, 4% | 3 | 8 year future survey interval | n/a | H |
| NF-E2 | II | 1, 4% | 4 | Exclude 1987/9 abundance in WI, EG & EI/F | M | M |
| NF-E3 | III | 1, 4% | 4 | Exclude 1987/9 abundance in WI, EG & EI/F | M | M |
| NF-D1 | I | 1% | 4 | Dispersal: max bound of 20% | M |   |
| NF-D3 | III | 1% | 4 | Dispersal: max bound of 20% | M |   |
| NF-J1 | II | 1, 4% | 4 | Assume _g_(0) = 0.8 (all estimates) | M | H |
| NF-J2 | III | 1, 4% | 4 | Assume _g_(0) = 0.8 (all estimates) | M | H |


