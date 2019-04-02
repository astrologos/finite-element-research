#! /bin/bash
#PBS -l walltime=00:05:00
#PBS -l nodes=1:ppn=8
#PBS -W group_list=blueridge
#PBS -q open_q
#PBS -j oe

cd $PBS_O_WORKDIR

module purge
module load gcc
module load openmpi
module load atlas
module load glm
module load lua
module load python
module load boost
module load phdf5
module load p4est
module load trilinos
module load dealii
module load cmake

cmake -DDEAL_II_DIR=$DEAL_II_DIR .

make release run > out.txt

rm -f ./CMakeCache.txt
rm -f ./cmake_install.cmake
rm -f -r ./CMakeFiles
rm -f ./Makefile
rm -f ./*~
rm -f ./*#
rm -f ./mark