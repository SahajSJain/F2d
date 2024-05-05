CXX=g++-13
CFLAGS=-Wall -Wextra -O3 -Werror -pedantic -Ofast -Winline

HDIR=./headers
DEPSNAME = INPUT.h FUNCTIONS.h FIELDVARS.h CLASSES.h COEFFVARS.h FLOWVARS.h 
DEPS = $(patsubst %,$(HDIR)/%,$(DEPSNAME))

ODIR=./obj
OBJNAME = main_program.o\
 			read_input_data.o\
  			write_output_data.o\
   			setup_field_variables.o\
   			boundary_conditions.o\
	 		calculators.o \
	 		velocity_predictor.o\
	  		solvers.o\
	  		fractional_step.o
OBJ = $(patsubst %,$(ODIR)/%,$(OBJNAME))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

all: F2d

F2d: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	@echo "Cleaning Up"
	rm -f f2d $(ODIR)/*.o *~ core $(INCDIR)/*~ 
