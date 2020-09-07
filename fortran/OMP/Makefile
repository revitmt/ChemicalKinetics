# -w : disable warnings

ifort         = ifort
ifort_options = -c -O3 -w -auto
# /Qprof-use /Qipo
# ifort_options = -c /Qcov-gen
# ifort_options = -c -debug:full -Od -Zi -MD -check:none /libs:dll -threads/dbglibs

OMP_options = -openmp

mkl_inc_path  = D:\Windows\Programs\IntelParallelXE\compilers_and_libraries_2017\windows\mkl\include
mkl_lib_path  = D:\Windows\Programs\IntelParallelXE\compilers_and_libraries_2017\windows\mkl\lib\intel64_win

# link_mkl_lib     = /link mkl_intel_lp64.lib mkl_core.lib mkl_sequential.lib /LIBPATH:$(mkl_lib_path)
link_mkl_lib  = /link mkl_intel_lp64.lib mkl_core.lib mkl_intel_thread.lib /LIBPATH:$(mkl_lib_path)


# run_simulation
all: module_ChemicalSystem.obj module_MKL_wrappers.obj module_methods.obj run_simulation.obj toms611.obj toms778.obj
	$(ifort) -o result $(OMP_options) module_ChemicalSystem.obj module_methods.obj module_MKL_wrappers.obj toms611.obj toms778.obj toms778_timer.obj toms778_linpack.obj run_simulation.obj  $(link_mkl_lib)
	del *.obj *.mod

module_ChemicalSystem.obj: module_ChemicalSystem.f90
	$(ifort) $(ifort_options) $(OMP_options) module_ChemicalSystem.f90 -I $(mkl_inc_path)

module_methods.obj: module_methods.f90 SSA.f90 ThetaTauLeap.f90 ssTauLeap.f90 ssTauLeap2.f90 ssTauLeap_old.f90 estimate_parameters.f90 estimate_stationary_parameters.f90 estimate_parameters_bistable.f90 ssTauLeapBistable.f90 
	$(ifort) $(ifort_options) $(OMP_options) module_methods.f90 -I $(mkl_inc_path)

module_MKL_wrappers.obj: module_MKL_wrappers.f90 MKL_fsolve.f90 MKL_poissrnd.f90 linspace.f90 fsolve.f90 kron.f90 linsolve.f90 gradsearch.f90
	$(ifort) $(ifort_options) $(OMP_options) module_MKL_wrappers.f90 -I $(mkl_inc_path)

toms611.obj: toms611.f90
	$(ifort) $(ifort_options) $(OMP_options) toms611.f90

toms778.obj: toms778_timer.f toms778_linpack.f toms778.f
	$(ifort) $(ifort_options) $(OMP_options) toms778_timer.f toms778_linpack.f toms778.f $(link_mkl_lib)

run_simulation.obj: run_simulation.f90
	$(ifort) $(ifort_options) $(OMP_options) run_simulation.f90 -I $(mkl_inc_path)

clean:
	del *.obj *.exe *.mod 

#*.txt