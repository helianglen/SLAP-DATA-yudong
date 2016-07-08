
#FFLAGS="-DDEBUG -g -u -traceback -fpe0 -nomixed_str_len_arg -names lowercase -convert big_endian -assume byterecl -I/home/ytian/proj-disk/libs/CRTM_Profile_Utility_intel_11_1_038/include"
FFLAGS="-DDEBUG -g -u -traceback -fpe0 -nomixed_str_len_arg -names lowercase -convert big_endian -assume byterecl" 

ifort -c $FFLAGS compute_sea_Tb.F90 sea_water.F90
ifort $FFLAGS -o compute_sea_Tb compute_sea_Tb.o sea_water.o /home/ytian/lswg/LIS-MEM/src/lis_mem.a

