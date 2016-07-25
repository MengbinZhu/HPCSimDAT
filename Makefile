FC=ifort
FCFLAGS=-g -traceback 

OBJS=Bdepart.o \
     Bdivers.o \
     Bimsl.o \
     Bnumadj.o \
     Bnumericas.o \
     Bsimplif.o \
     Bsouspec.o \
     Btransfo.o \
     Bvarie.o \
     Bprocess.o \
     Brandom.o \
     BInv.o \
     BEnKF.o \
     BEOF.o

all: qgcont_da.x

qgcont_da.x: $(OBJS)
	$(FC) $(FCFLAGS) -o qgcont_da.x qgcont_da.f90 $(OBJS) 

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	rm -f *.o qgcont_da.x
