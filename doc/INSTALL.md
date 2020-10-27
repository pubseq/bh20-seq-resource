#+OPTIONS: ^:nil

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
guix environment guix --ad-hoc git python openssl python-pycurl python-magic nss-certs python-pyshex -- python3 bh20sequploader/main.py example/sequence.fasta example/maximum_metadata_example.yaml
```

Note that python-pyshex is packaged in
http://git.genenetwork.org/guix-bioinformatics/guix-bioinformatics

so you'll need it to the GUIX_PACKAGE_PATH - see the README in that
repository. E.g.

```sh
env GUIX_PACKAGE_PATH=~/iwrk/opensource/guix/guix-bioinformatics/ ~/opt/guix/bin/guix environment -C guix --ad-hoc git python python-flask python-pyyaml python-pycurl python-magic  nss-certs python-pyshex python-pyyaml --network openssl python-pyshex python-pyshexc minimap2 python-schema-salad python-arvados-python-client --share=/export/tmp -- env TMPDIR=/export/tmp python3 bh20sequploader/main.py --help
```

### Using the Web Uploader

To run the web uploader in a GNU Guix environment/container run it with something like

```
guix environment guix --ad-hoc git python python-flask python-pyyaml python-pycurl python-magic  nss-certs --network openssl -- env FLASK_ENV=development PYTHONPATH=$PYTHONPATH:./bh20sequploader FLASK_APP=bh20simplewebuploader/main.py flask run
 * Serving Flask app "bh20simplewebuploader/main.py"
 * Environment: production
   WARNING: This is a development server. Do not use it in a production deployment.
   Use a production WSGI server instead.
 * Debug mode: off
 * Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)
```

WIP: add gunicorn container

Currently the full webserver container deploy command looks like

```
penguin2:~/iwrk/opensource/code/vg/bh20-seq-resource$  env GUIX_PACKAGE_PATH=~/iwrk/opensource/guix/guix-oinformatics/ ~/iwrk/opensource/guix/guix/pre-inst-env guix environment -C guix --ad-hoc git python python-flask python-pyyaml python-pycurl python-magic  nss-certs python-pyshex python-pyyaml --network openssl python-pyshex python-pyshexc minimap2 python-schema-salad python-arvados-python-client --share=/export/tmp -- env TMPDIR=/export/tmp FLASK_ENV=development FLASK_APP=bh20simplewebuploader/main.py flask run
```

Note: see above on GUIX_PACKAGE_PATH.
