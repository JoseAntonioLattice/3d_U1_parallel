src_files = constants number2string check_files_directories random matlib statistics parallel_utils indices files pbc arrays parameters U1_functions hybridMC lua dynamics main
obj_files = $(patsubst %, build/%.o, $(src_files) )

LIB = ~/Fortran/lib/
INC = ~/Fortran/include/


ifeq ($(MODE), SERIAL)
	program=3d_U1_serial
	PRE=SERIAL
	FC=gfortran
endif


ifeq ($(PARALLEL),1)
	program=3d_U1_parallel
	PRE="PARALLEL=1"
	FC=caf
endif

ifeq ($(PARALLEL),2)
	program=3d_U1_parallel
	PRE="PARALLEL=2"
	FC=caf
endif


build/$(program): $(obj_files)
	$(FC) $^ -o $@ #-L $(LIB) -lrandom -lconstants -lfiles -lnum2str -lstats

build/main.o: app/main.F90
	$(FC) -D$(PRE) -O3 -c -I build -I $(INC) $< -o $@

build/%.o: src/%.F90
	$(FC) -D$(PRE) -O3 -c -J build -I build -I $(INC) $< -o $@

run:
	{ echo $(c1) $(c2) $(c3) ; echo input/input_parameters.nml; echo data_$(c1)x$(c2)x$(c3)_$(PARALLEL).dat; } | cafrun -n $$(( $(c1)*$(c2)*$(c3) )) build/$(program) 


run_test:
	gfortran -J build src/indices.f90 test/test.f90 -o build/test
	build/test

run_serial:
	{ echo input/input_parameters.nml; echo data_serial.dat; } | LD_LIBRARY_PATH=$(LIB) build/$(program) 


clean:
	rm -rf build/*
help:
	echo $(obj_files)

documentation:
	pandoc -f gfm README.md -o doc/README.pdf
