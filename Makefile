FC = gfortran
FCFLAGS = -g -fcheck=all -Warray-bounds

# source files and objects
SRCS = limiteurs.f90 module.f90 main.f90

# dependances
#main.f90: Tools.o

# program name
PROGRAM = run

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod
