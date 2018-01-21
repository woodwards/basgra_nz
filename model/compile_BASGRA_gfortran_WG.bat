gfortran -x f95-cpp-input -Dweathergen -O3 -c -fdefault-real-8 parameters_site.f90 parameters_plant.f90 environment.f90 resources.f90 soil.f90 plant.f90 set_params.f90 BASGRA.f90 
gfortran -shared -o BASGRA_WG.DLL parameters_site.o parameters_plant.o environment.o resources.o soil.o plant.o set_params.o BASGRA.o
del *.o
del *.mod
pause
