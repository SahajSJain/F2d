CXX=icpx
CFLAGS=  -Wall -O3 -pedantic -Ofast -qopenmp
# -Wextra  -Wall -Wextra -O3 -Werror -pedantic -Ofast -Winline -std=c++2c
HDIR=./headers
DEPSNAME = INPUT.h FUNCTIONS.h FIELDVARS.h CLASSES.h COEFFVARS.h FLOWVARS.h IBMVARS.h 
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
	  		fractional_step.o\
	  		tecplot_printer.o\
			fractional_step_looper.o\
			debug_export.o\
			ssm_setup.o\
			ssm_bc.o
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
