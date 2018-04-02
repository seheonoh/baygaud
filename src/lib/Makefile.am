# src/lib/Makefile.am
# by Se-Heon Oh (KASI/ICRAR)

# cfitsio
DIR_CFITSIO = /mnt/g0/research/packages/cfitsio/cfitsio
# gsl
DIR_GSL = /mnt/g0/research/packages/gsl/gsl-1.16/.libs

# additional include paths necessary to compile baygaud
AM_CFLAGS= -O3 -DMPI -fopenmp -m64 -g -I$(top_srcdir)/src/lib -I$(DIR_CFITSIO) -I$(DIR_GSL) -mcmodel=large

# the library names to build (note we are building static libs only)
lib_LIBRARIES = libbaygaud.a

# where to install the headers on the system
libbaygaud_adir = $(includedir)/lib

# the list of header files that belong to the library (to be installed later)
libbaygaud_a_HEADERS = baygaud.2dmaps.h baygaud.global_params.h baygaud.multinest.h baygaud.cfitsio.h baygaud.gsl.h baygaud.main.h baygaud.sort.h baygaud.memory.h baygaud.etc.h baygaud.mpi.h baygaud.gfit.h baygaud.mpi_datatype.h

# the sources to add to the library and to add to the source distribution
libbaygaud_a_SOURCES= $(libmylib_a_HEADERS) baygaud.2dmaps.c baygaud.etc.c baygaud.mpi.c baygaud.gfit.c baygaud.multinest.c baygaud.gsl.c baygaud.memory.c baygaud.sort.c 

# end of line
