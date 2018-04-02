#!/bin/bash
# setup_baygaud.sh 
# by Se-Heon Oh (KASI/ICRAR)
# Version v1.2.0

#......................
# PRE-STEP FOR AUTOCONF
#......................

#--------------------------------------------------------
# EDIT ./src/bin/Makefile.am
#--------------------------------------------------------
# 1. MULTINEST DIRECTORY from the user
echo ""
echo "------------------------------------------------------------------------"
echo "... Setup the paths to the libraries required for BAYGAUD installation ..."
echo "------------------------------------------------------------------------"
echo ""
echo ""
default_path="/mnt/g0/research/packages/multinest/MultiNest_v3.7_MPIoff"
read -p "--> 1/3. Enter the absolute path to the directory where MULTINEST_MPIOFF is installed <--
    [Check with 'locate libnest3_mpioff.a']
    [Press 'Enter' to use the default directory '/mnt/g0/research/packages/multinest/MultiNest_v3.7_MPIoff' 
    or enter another directory and press 'Enter']: " dir_multinest_mpioff
dir_multinest_mpioff="${dir_multinest_mpioff:-$default_path}"
sed -i "s|^DIR_MPIOFF_MULTINEST.*|DIR_MPIOFF_MULTINEST = ${dir_multinest_mpioff}|g" ./src/bin/Makefile.am
echo "    >> DIR_MPIOFF_MULTINEST= ${dir_multinest_mpioff=}"


# 2. CFITSIO DIRECTORY from the user
echo ""
default_path="/mnt/g0/research/packages/cfitsio/cfitsio"
read -p "--> 2/3. Enter the absolute path to the directory where CFITSIO is installed <--
    [Check with 'locate libcfitsio.a']
    [Press 'Enter' to use the default directory '/mnt/g0/research/packages/cfitsio/cfitsio' 
    or enter another directory and press 'Enter']: " dir_cfitsio
dir_cfitsio="${dir_cfitsio:-$default_path}"
sed -i "s|^DIR_CFITSIO.*|DIR_CFITSIO = ${dir_cfitsio}|g" ./src/bin/Makefile.am
echo "    >> DIR_CFITSIO= ${dir_cfitsio=}"


# 3. GSL DIRECTORY from the user
echo ""
default_path="/mnt/g0/research/packages/gsl/gsl-1.16/.libs"
read -p "--> 3/3. Enter the absolute path to the directory where GSL is installed <--
    [Check with 'locate libgsl.a']
    [Press 'Enter' to use the default directory '/mnt/g0/research/packages/gsl/gsl-1.16/.libs' 
    or enter another directory and press 'Enter']: " dir_gsl
dir_gsl="${dir_gsl:-$default_path}"
sed -i "s|^DIR_GSL.*|DIR_GSL = ${dir_gsl}|g" ./src/bin/Makefile.am
echo "    >> DIR_GSL= ${dir_gsl=}"

#--------------------------------------------------------
# EDIT ./src/lib/Makefile.am
#--------------------------------------------------------
# Update the path of the CFITSIO DIRECTORY obtained from the user
sed -i "s|^DIR_CFITSIO.*|DIR_CFITSIO = ${dir_cfitsio}|g" ./src/lib/Makefile.am
# Update the path of the GSL DIRECTORY obtained from the user
sed -i "s|^DIR_GSL.*|DIR_GSL = ${dir_gsl}|g" ./src/lib/Makefile.am

echo ""
echo "-------------------------------------------------"
echo "... The paths to the libraries are configured ..."
echo "-------------------------------------------------"
echo ""

echo ""
echo "----------------------"
echo "... Start AUTOCONF ..."
echo "----------------------"
echo ""
#......................
# AUTOCONF STEP
#......................
aclocal
autoheader
autoconf
#touch NEWS README AUTHORS ChangeLog
#libtoolize --automake --copy --force
automake --foreign --copy --add-missing
./configure

echo ""
echo "---------------------------------------------------------------------------"
echo "... Makefile is configured... Now, type 'make' in the current directory ..."
echo "---------------------------------------------------------------------------"
echo ""


