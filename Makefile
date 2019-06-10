OBJS =  MOLECULE.o MARS.o ATOM.o MATRIX.o OPT.o
CC = g++

MARS : ${OBJS}
	${CC} ${OPTFLAGS} -o $@ ${OBJS} -limf
MATRIX.o : MATRIX.cpp MATRIX.h
ATOM_GA.o : ATOM.cpp ATOM.h
OPT.o : OPT.cpp OPT.h ATOM.h MATRIX.h
MOLECULE_GA.o : MOLECULE.cpp ATOM.h MOLECULE.h MATRIX.h OPT.h
MARS.o : MARS.cpp ATOM.h MOLECULE.h

%.o : %.cpp
	${CC} ${OPTFLAGS} ${CDEBUG} -o $@ -c $< ${LIBFLAGS} 

clean:
	rm -f *.o out core *.out

