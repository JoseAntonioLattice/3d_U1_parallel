
FC = gfortran #caf
program = 3d_U1
src_files = indices pbc arrays parameters U1_functions lua dynamics main
obj_files = $(patsubst %, build/%.o, $(src_files) )

LIB = ~/Fortran/lib/
INC = ~/Fortran/include/

build/$(program): $(obj_files)
	$(FC) $^ -o $@ -L $(LIB) -lrandom -lconstants

build/main.o: app/main.F90
	$(FC) -c -I build -I $(INC) $< -o $@

build/%.o: src/%.F90
	$(FC) -c -J build -I build -I $(INC) $< -o $@

run:
	{ echo 1 1 2; echo input/input_parameters.nml; echo data/data_1x1x2.dat; } | LD_LIBRARY_PATH=$(LIB) cafrun -n 2 build/3d_U1 


run_test:
	gfortran -J build src/indices.f90 test/test.f90 -o build/test
	build/test

run_serial:
	{ echo input/input_parameters.nml; echo data/data_1x1x2.dat; } | LD_LIBRARY_PATH=$(LIB) build/3d_U1 


clean:
	rm -rf build/*
help:
	echo $(obj_files)

