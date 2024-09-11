# docker build -t urva:latest .
# docker run -it -v "$PWD:$PWD" -w "$PWD" urva:latest inp.inp
# docker run -it --entrypoint /bin/bash urva:latest

# Set build image
FROM debian:bookworm-slim AS build

# Install build tooling
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update &&\
 apt-get install -y\
 build-essential\
 gfortran

# Build lmodea executable
COPY lib/LModeA/src lmodea
RUN cd lmodea &&\
  make all

# Deploy image
FROM python:3.8.18-slim-bookworm
LABEL org.opencontainers.image.source="https://github.com/SouthernMethodistUniversity/URVA"

# Install Python dependencies
COPY requirements.txt /
RUN pip install -r /requirements.txt &&\
 rm /requirements.txt

# Copy Python scripts and LModeA executable
WORKDIR /urva
COPY --from=build lmodea/lmodea.exe lm90.test.exe
COPY src .

# Generate Python cache and remove source
RUN python3 -m compileall -b . &&\
 rm *.py

# Run URVA scripts
ENTRYPOINT ["python3", "/urva/main.pyc"]

