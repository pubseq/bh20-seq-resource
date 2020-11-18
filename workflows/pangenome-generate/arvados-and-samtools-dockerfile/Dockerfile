FROM debian:10
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update -q
RUN apt-get install -yq --no-install-recommends gnupg
ADD 1078ECD7.key /tmp/
RUN cat /tmp/1078ECD7.key | apt-key add -
RUN echo 'deb http://apt.arvados.org/ buster main' > /etc/apt/sources.list.d/apt.arvados.org-stable.list
RUN apt-get update -q && apt-get install -yq --no-install-recommends samtools python3-python-client
RUN rm -f /usr/bin/python && ln -s /usr/share/python3/dist/python3-python-client/bin/python /usr/bin/python
RUN rm -f /usr/bin/python3 && ln -s /usr/share/python3/dist/python3-python-client/bin/python /usr/bin/python3
