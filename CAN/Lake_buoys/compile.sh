
#FFLAGS="-DDEBUG -g -u -traceback -fpe0 -nomixed_str_len_arg -names lowercase -convert big_endian -assume byterecl -I/home/ytian/proj-disk/libs/CRTM_Profile_Utility_intel_11_1_038/include"
FFLAGS="-DDEBUG -g -u -traceback -fpe0 -nomixed_str_len_arg -names lowercase -convert big_endian -assume byterecl" 

ifort -c $FFLAGS compute_buoy_Tb.F90
ifort $FFLAGS -o compute_buoy_Tb compute_buoy_Tb.F90 /home/ytian/lswg/LIS-MEM/src/lis_mem.a

