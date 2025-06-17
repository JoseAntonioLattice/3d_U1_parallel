
FC = caf
program = 3d_U1
src_files = indices arrays parameters U1_functions lua dynamics main
obj_files = $(patsubst %, build/%.o, $(src_files) )

build/$(program): $(obj_files)
	$(FC) $^ -o $@

build/main.o: app/main.f90
	$(FC) -c -I build $< -o $@

build/%.o: src/%.f90
	$(FC) -c -J build $< -o $@

run_test:
	gfortran -J build src/indices.f90 test/test.f90 -o build/test
	build/test


clean:
	rm -rf build/*
help:
	echo $(obj_files)

