# docker build -t purva:latest .
# docker run -it -v "$PWD:$PWD" -w "$PWD" purva:latest inp.inp
# docker run -it --entrypoint /bin/bash purva:latest

# Label for GitHub Container Registry
LABEL org.opencontainers.image.source=https://github.com/SouthernMethodistUniversity/pURVA_release

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

# Install Python dependencies
COPY requirements.txt /
RUN pip install -r /requirements.txt &&\
 rm /requirements.txt

# Copy Python scripts and LModeA executable
WORKDIR /purva
COPY --from=build lmodea/lmodea.exe lm90.test.exe
COPY src .

# Run pURVA scripts
ENTRYPOINT ["python3", "/purva/main.py"]

