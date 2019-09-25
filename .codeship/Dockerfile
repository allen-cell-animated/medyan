FROM ubuntu:latest
COPY . /medyan
WORKDIR /medyan/src

ENV DEBIAN_FRONTEND noninteractive

# Install gcc
RUN apt-get -qq update && apt-get install -q -y software-properties-common
RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y
RUN apt-get -qq update && apt-get install -q -y gcc-7
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 60 --slave /usr/bin/g++ g++ /usr/bin/g++-7

# Install boost
RUN apt-get install -q -y libboost-all-dev

# Get google test source, and use CMake to build it
# RUN apt-get install -q -y libgtest-dev
# RUN apt-get install -q -y cmake
# RUN cmake /usr/src/gtest/CMakeLists.txt
# RUN make -C /usr/src/gtest

# RUN mkdir -p ../build/gtest/include ../build/gtest/lib
# RUN cp -a /usr/src/gtest/include/gtest/ ../build/gtest/include
# RUN cp /usr/src/gtest/*.a ../build/gtest/lib
 
CMD ["make"]
