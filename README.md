# FEM Research

Optimizing deal.ii routines for parallel performance.

Notes:
- CMakeLists.txt has been altered to take any deal.ii code named "mark.cc"
- Project is run in release mode
- Output is saved to "out.txt"

load.sh does the following:
- Requests cluster resources
- Loads dependencies and the deal.ii library
- Makes and runs project using cluster resources
- Deletes CMake cache files
- Deletes all emacs backup files
- Deletes binary

Instructions:
- Copy target program to "mark.cc"
- Edit load.sh for preferred cluster resources
- Use "qsub load.sh" to submit a PBS batch
- Results are saved to "out.txt"
