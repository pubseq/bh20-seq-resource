# INSTALLATION

Other options for running this tool.

## GNU Guix

### Running the CLI uploader

Another way to install this tool is inside a [GNU Guix Environment](https://guix.gnu.org/manual/en/html_node/Invoking-guix-environment.html), which can handle installing dependencies for you

1. **Set up and enter a Guix environment with the necessary dependencies.** After installing Guix run:

```sh
guix environment -C guix --ad-hoc git python openssl python-pycurl nss-certs
```

2. **Install the tool.** From there you can follow the [user installation instructions](#installation-with-pip3---user). In brief:

```sh
pip3 install --user schema-salad  arvados-python-client
```

Pip installed the following modules

```
arvados-python-client-2.0.1 ciso8601-2.1.3 future-0.18.2 google-api-python-client-1.6.7 httplib2-0.17.1 oauth2client-4.1.3 pyasn1-0.4.8 pyasn1-modules-0.2.8 rsa-4.0 ruamel.yaml-0.15.77 six-1.14.0 uritemplate-3.0.1 ws4py-0.5.1
```

3. Run the tool directly with

```sh
guix environment guix --ad-hoc git python openssl python-pycurl nss-certs -- python3 bh20sequploader/main.py
```

### Using the Web Uploader

To run the web uploader in a GNU Guix environment

```
guix environment guix --ad-hoc git python python-flask python-pyyaml nss-certs --network openssl -- env FLASK_APP=bh20simplewebuploader/main.py flask run
```

The containerized version looks like

```
guix environment -C guix --ad-hoc git python python-flask python-pyyaml nss-certs --network openssl
```

and

```
env FLASK_APP=bh20simplewebuploader/main.py flask run
```

WIP: add gunicorn container
