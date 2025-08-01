cmake_minimum_required(VERSION 3.28)
include_directories(~/Fortran/include)
link_directories(~/Fortran/lib)
project(3d_U1 Fortran)
add_compile_options(-O3 -DSERIAL)
add_executable(${PROJECT_NAME}.exe app/main.F90 src/arrays.F90 src/dynamics.F90 
                                   src/files.F90 src/hybridMC.F90 src/indices.F90 src/lua.F90
                                   src/parameters.F90 src/pbc.F90 src/U1_functions.F90)
