PRE=PARALLEL

ifeq ($(PRE),SERIAL)
	FC=gfortran
	program = 3d_U1_serial
endif

ifeq ($(PRE),PARALLEL)
	FC=caf
	program = 3d_U1_parallel	
endif

src_files = indices pbc arrays parameters U1_functions lua dynamics main
obj_files = $(patsubst %, build/%.o, $(src_files) )

LIB = ~/Fortran/lib/
INC = ~/Fortran/include/

build/$(program): $(obj_files)
	$(FC) $^ -o $@ -L $(LIB) -lrandom -lconstants

build/main.o: app/main.F90
	$(FC) -D$(PRE) -O3 -c -I build -I $(INC) $< -o $@

build/%.o: src/%.F90
	$(FC) -D$(PRE) -O3 -c -J build -I build -I $(INC) $< -o $@

run:
	{ echo 1 2 2; echo input/input_parameters.nml; echo data/data_1x1x2.dat; } | LD_LIBRARY_PATH=$(LIB) cafrun -n 4 build/$(program) 


run_test:
	gfortran -J build src/indices.f90 test/test.f90 -o build/test
	build/test

run_serial:
	{ echo input/input_parameters.nml; echo data/data_serial.dat; } | LD_LIBRARY_PATH=$(LIB) build/$(program) 


clean:
	rm -rf build/*
help:
	echo $(obj_files)

documentation:
	pandoc -f gfm README.md -o doc/README.pdf
