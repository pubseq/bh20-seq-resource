# Dockerfile for containerizing the web interface
FROM python:3.6-jessie
WORKDIR /app

RUN pip3 install gunicorn

ADD LICENSE /app/
ADD gittaggers.py /app/
ADD setup.py /app/
ADD README.md /app/
ADD example /app/example
ADD bh20seqanalyzer /app/bh20simplewebuploader
ADD bh20sequploader /app/bh20sequploader
ADD bh20simplewebuploader /app/bh20simplewebuploader

RUN pip3 install -e .[web]

ENV PORT 8080
CMD ["gunicorn", "-w", "4", "-b", "0.0.0.0:8080", "bh20simplewebuploader.main:app"]
