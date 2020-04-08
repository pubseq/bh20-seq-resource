# Sequence uploader

This repository provides a sequence uploader for the COVID-19 Virtual Biohackathon's Public Sequence Resource project. You can use it to upload the genomes of SARS-CoV-2 samples to make them publicly and freely available to other researchers.

To get started, first [install the uploader](#installation), and use the `bh20-seq-uploader` command to [upload your data](#usage).

# Installation

There are several ways to install the uploader. The most portable is with a [virtualenv](#installation-with-virtualenv).

## Installation with `virtualenv`

1. **Prepare your system.** You need to make sure you have Python, and the ability to install modules such as `pycurl` and `pyopenssl`. On Ubuntu 18.04, you can run:

```sh
sudo apt update
sudo apt install -y virtualenv git libcurl4-openssl-dev build-essential python3-dev libssl-dev
```

2. **Create and enter your virtualenv.** Go to some memorable directory and make and enter a virtualenv:

```sh
virtualenv --python python3 venv
. venv/bin/activate
```

Note that you will need to repeat the `. venv/bin/activate` step from this directory to enter your virtualenv whenever you want to use the installed tool.

3. **Install the tool.** Once in your virtualenv, install this project:

```sh
pip3 install git+https://github.com/arvados/bh20-seq-resource.git@master
```

4. **Test the tool.** Try running:

```sh
bh20-seq-uploader --help
```

It should print some instructions about how to use the uploader.

**Make sure you are in your virtualenv whenever you run the tool!** If you ever can't run the tool, and your prompt doesn't say `(venv)`, try going to the directory where you put the virtualenv and running `. venv/bin/activate`. It only works for the current terminal window; you will need to run it again if you open a new terminal.

## Installation with `pip3 --user`

If you don't want to have to enter a virtualenv every time you use the uploader, you can use the `--user` feature of `pip3` to install the tool for your user.

1. **Prepare your system.** Just as for the `virtualenv` method, you need to install some dependencies. On Ubuntu 18.04, you can run:

```sh
sudo apt update
sudo apt install -y virtualenv git libcurl4-openssl-dev build-essential python3-dev libssl-dev
```

2. **Install the tool.** You can run:

```sh
pip3 install --user git+https://github.com/arvados/bh20-seq-resource.git@master
```

3. **Make sure the tool is on your `PATH`.** THe `pip3` command will install the uploader in `.local/bin` inside your home directory. Your shell may not know to look for commands there by default. To fix this for the terminal you currently have open, run:

```sh
export PATH=$PATH:$HOME/.local/bin
```

To make this change permanent, assuming your shell is Bash, run:

```sh
echo 'export PATH=$PATH:$HOME/.local/bin' >>~/.bashrc
```

4. **Test the tool.** Try running:

```sh
bh20-seq-uploader --help
```

It should print some instructions about how to use the uploader.

## Installation from Source for Development

If you plan to contribute to the project, you may want to install an editable copy from source. With this method, changes to the source code are automatically reflected in the installed copy of the tool.

1. **Prepare your system.** On Ubuntu 18.04, you can run:

```sh
sudo apt update
sudo apt install -y virtualenv git libcurl4-openssl-dev build-essential python3-dev libssl-dev
```

2. **Clone and enter the repository.** You can run:

```sh
git clone https://github.com/arvados/bh20-seq-resource.git
cd bh20-seq-resource
```

3. **Create and enter a virtualenv.** Go to some memorable directory and make and enter a virtualenv:

```sh
virtualenv --python python3 venv
. venv/bin/activate
```

Note that you will need to repeat the `. venv/bin/activate` step from this directory to enter your virtualenv whenever you want to use the installed tool.

4. **Install the checked-out repository in editable mode.** Once in your virtualenv, install with this special pip command:

```sh
pip3 install -e .
```

5. **Test the tool.** Try running:

```sh
bh20-seq-uploader --help
```

It should print some instructions about how to use the uploader.

## Installation with GNU Guix

Another way to install this tool is inside a [GNU Guix Environment](https://guix.gnu.org/manual/en/html_node/Invoking-guix-environment.html), which can handle installing dependencies for you even when you don't have root access on an Ubuntu system.

1. **Set up and enter a container with the necessary dependencies.** After installing Guix as `~/opt/guix/bin/guix`, run:

```sh
~/opt/guix/bin/guix environment -C guix --ad-hoc git python openssl python-pycurl nss-certs
```
   
2. **Install the tool.** From there you can follow the [user installation instructions](#installation-with-pip3---user). In brief:

```sh
pip3 install --user git+https://github.com/arvados/bh20-seq-resource.git@master
```

# Usage

Run the uploader with a FASTA file and accompanying metadata file in [JSON-LD format](https://json-ld.org/):

```sh
bh20-seq-uploader example/sequence.fasta example/metadata.json
```

## Workflow for Generating a Pangenome

All these uploaded sequences are being fed into a workflow to generate a [pangenome](https://academic.oup.com/bib/article/19/1/118/2566735) for the virus. You can replicate this workflow yourself.

Get your SARS-CoV-2 sequences from GenBank in `seqs.fa`, and then run:

```sh
minimap2 -cx asm20 -X seqs.fa seqs.fa >seqs.paf
seqwish -s seqs.fa -p seqs.paf -g seqs.gfa
odgi build -g seqs.gfa -s -o seqs.odgi
odgi viz -i seqs.odgi -o seqs.png -x 4000 -y 500 -R -P 5
```

For more information on building pangenome models, [see this wiki page](https://github.com/virtual-biohackathons/covid-19-bh20/wiki/Pangenome#pangenome-model-from-available-genomes).


