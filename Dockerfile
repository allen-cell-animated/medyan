FROM ubuntu:20.04

WORKDIR /root/medyan
RUN mkdir input output

# Copy MEDYAN executable
COPY build/medyan .

# Install necessary libraries
RUN apt-get update && apt-get install -y \
    libxinerama-dev \
    libxcursor-dev \
    xorg-dev \
    libglu1-mesa-dev

# Run MEDYAN
ENTRYPOINT ["./medyan", "-i", "/home/input", "-o", "/home/output", "-s", "/home/input/systeminput.txt"]