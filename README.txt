======Introduction======

MARS, Molecular Assembling and Representation Suite, is a program for
general purpose computer aided molecular design. This program uses an 
integer arry data structure to store the information of atom type, 
atom connectivity, and bond order of a molecule. There are also 
subroutines that can be used to modify the data structure in order to
generate new molecules. This program is espicially useful for
computer molecular design problems.

The program consists of 4 header files and 5 cpp files
ATOM.h  MATRIX.h  MOLECULE.h  OPT.h
MARS.cpp ATOM.cpp  MATRIX.cpp  MOLECULE.cpp  OPT.cpp

There are additional 4 text files for running examples molecular operations
on butane, methane, and phenol
butane.gjf  methane.gjf  phenol.gjf example.in

======Compile the program======

This program is written in C++. In Linux/Unix environment, run
./compile
or
make
to create the executable (MARS) for the program


======Usage======
1. basic use (generate new molecules using butane and phenol)
use the example files provided and execute
./MARS example.in

The input file, example.in, contain the filenames of two molecules,
butane and phenol. The program will expect the 3D coordinates of
these two molecules (butane.gjf and phenol.gjf in Gaussian input 
file format) to be in the same folder as the example.in file. 
The program will perform molecular operations on these two molecules.
The results will be produced in 6 output files
combination.txt crossover.txt subtraction.txt addition.txt exchange.txt ring.txt
These files contain all possible molecules that can be generated from 
each corresponding molecular operation. The results are shown in SMILES
format. The structure of the generated molecules can be visualized using 
other programs, such as Avogadro (https://avogadro.cc/).

Furthermore, the program also generates benzene from methane. The intermediate
molecules are produced in benzene.txt. 

A log file, named molecule.log, will be produced to note any warning and error
encountered during the molecular operations. 
A Warning message is produced when a certain operation cannot be performed
successfully.
A Error message is produced when the program needs to stop because of wrong
usage (e.g., missing input files)

2. Advanced usage 
The program accept keywords, gen3d and genmds (both optional), to request the
generation of 3D coordinates and MDS structure in the output files. For
example,
./MARS example.in gen3d
./MARS example.in genmds
./MARS example.in gen3d genmds

======Developers======
This program is developed by Hsuan Hao Hsu, Chen Hsuan Huang, and Shiang Tai Lin.
Computational Molecular Engineering Laboratory
Department of Chemical Engineering, National Taiwan University, Taipei, Taiwan
stlin@ntu.edu.tw
