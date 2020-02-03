#!/usr/bin/env bash

(
  echo "This is the old configuration file for Deepthought2."
  echo "If you are using CMake, use conf-dt2.sh instead."
  makefile_file="$(dirname "$0")/src/Makefile"
  sed -i -e 's:\(-I/usr/include/boost\) :\1 -I$(BOOST_INCLUDE) :g' \
         -e 's:\(-lboost_system\):\1 -L$(BOOST_LIB) -L$(GCC_ROOTDIR)/lib64 -Wl,-rpath,$(GCC_ROOTDIR)/lib64:g' \
         "$makefile_file"
)

