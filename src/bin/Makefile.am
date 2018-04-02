# src/bin/Makefile.am
# by Se-Heon Oh (KASI/ICRAR)

#-----------------------------------------------------
# additional libraries linked to baygaud
#-----------------------------------------------------

# 1/3. multinest
DIR_MPIOFF_MULTINEST = /mnt/g0/research/packages/multinest/MultiNest_v3.7_MPIoff

#++++++++++
# 2/3. cfitsio
DIR_CFITSIO = /mnt/g0/research/packages/cfitsio/cfitsio

#++++++++++
# 3/3. gsl
DIR_GSL = /mnt/g0/research/packages/gsl/gsl-1.16/.libs


AM_CFLAGS= -O3 -DMPI -fopenmp -m64 -g -I$(top_srcdir)/src/lib -I$(DIR_CFITSIO) -I$(DIR_GSL) -mcmodel=large

# the program to build (the names of the final binaries)
bin_PROGRAMS= baygaud

# list of sources for the 'baygaud' binary
baygaud_SOURCES= baygaud.c

# the additional libraries needed to link baygaud
baygaud_LDADD= -L$(top_srcdir)/src/lib -lbaygaud -L$(DIR_MPIOFF_MULTINEST) -lnest3_mpioff -L$(DIR_CFITSIO) -lcfitsio -L$(DIR_GSL) -lgsl -lpthread -lm -lpthread

# end of line
