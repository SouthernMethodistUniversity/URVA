# docker build -t purva:latest .
# docker run -it -v "$PWD:$PWD" -w "$PWD" purva:latest main.py inp.inp
# docker run -it --entrypoint /bin/bash purva:latest
# docker save -o ~/Downloads/lasso.tar purva:latest

FROM --platform=linux/amd64 python:2.7.18-slim-buster

COPY requirements.txt /
RUN pip install -r /requirements.txt &&\
 rm /requirements.txt

ENTRYPOINT ["python"]

