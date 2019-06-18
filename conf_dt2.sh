#!/usr/bin/env bash

(
  makefile_file="$(dirname "$0")/src/Makefile"
  sed -i -e 's:\(-I/usr/include/boost\) :\1 -I$(BOOST_INCLUDE) :g' \
         -e 's:\(-lboost_system\):\1 -L$(BOOST_LIB) -L$(GCC_ROOTDIR)/lib64 -Wl,-rpath,$(GCC_ROOTDIR)/lib64:g' \
         "$makefile_file"
)

