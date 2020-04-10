---
title: 'CPSR: COVID-19 Public Sequence Resource'
title_short: 'CPSR: COVID-19 Public Sequence Resource'
tags:
  - Sequencing
  - COVID-19
authors:
  - name: Pjotr Prins
    orcid: 0000-0002-8021-9162
    affiliation: 1
  - name: Peter Amstutz
    orcid: 0000
    affiliation: 2
  - name: Tazro Ohta
    orcid: 0000
    affiliation: 3
  - name: Thomas Liener
    orcid: 0000
    affiliation: 4
  - name: Erik Garrison
    orcid: 0000
    affiliation: 5
  - name: Michael Crusoe
    orcid: 0000
    affiliation: 6
  - name: Rutger Vos
    orcid: 0000
    affiliation: 7
  - name: Michael Heuer
    orcid: 0000
    affiliation: 8
  - name: Adam Novak
    orcid: 0000
    affiliation: 9
  - name: Alex Kanitz
    orcid: 0000
    affiliation: 10
  - name: Jerven Bolleman
    orcid: 0000
    affiliation: 11
  - name: Joep de Ligt
    orcid: 0000
    affiliation: 12
affiliations:
  - name: Department of Genetics, Genomics and Informatics, The University of Tennessee Health Science Center, Memphis, TN, USA.
    index: 1
  - name: Curii, Boston, USA
    index: 2
date: 11 April 2020
event: COVID2020
group: Public Sequence Uploader
authors_short: Pjotr Prins & Peter Amstutz \emph{et al.}
bibliography: paper.bib
---

<!--

The paper.md, bibtex and figure file can be found in this repo:

  https://github.com/arvados/bh20-seq-resource

To modify, please clone the repo. You can generate PDF of the paper by
pasting above link (or yours) with

  https://github.com/biohackrxiv/bhxiv-gen-pdf

Note that author order will change!

-->

# Introduction

As part of the COVID-19 Biohackathion 2020 we formed a working
group to create a COVID-19 Public Sequence Resource (CPSR) for
Corona virus sequences. The general idea was to create a
repository that has a low barrier to entry for uploading sequence
data using best practices. I.e., data published with a creative
commons 4.0 (CC-4.0) license with metadata using state-of-the art
standards and, perhaps most importantly, providing standardized
workflows that get triggered on upload, so that results are
immediately available in standardized data formats.

Existing data repositories for viral data include GISAID, EBI ENA
and NCBI. These repositories allow for free sharing of data, but
do not add value in terms of running immediate
computations. Also, GISAID, at this point, has the most complete
collection of genetic sequence data of influenza viruses and
related clinical and epidemiological data through its
database. But, due to a restricted license, data submitted to
GISAID can not be used for online web services and on-the-fly
computation. In addition GISAID registration which can take weeks
and, painfully, forces users to download sequences one at a time
to do any type of analysis. In our opinion this does not fit a
pandemic scenario where fast turnaround times are key and data
analysis has to be agile.

We managed to create a useful sequence uploader utility within
one week by leveraging existing technologies, such as the Arvados
Cloud platform [@Arvados], the Common Workflow Langauge (CWL)
[@CWL], Docker images built with Debian packages, and the many
free and open source software packages that are available for
bioinformatics.

The source code for the CLI uploader and web uploader can be
found [here](https://github.com/arvados/bh20-seq-resource)
(FIXME: we'll have a full page). The CWL workflow definitions can
be found [here](https://github.com/hpobio-lab/viral-analysis) and
on CWL hub (FIXME).

<!--

    RESULTS!

    For each section below

    State the problem you worked on
    Give the state-of-the art/plan
    Describe what you have done/results starting with The working group created...
    Write a conclusion
    Write up any future work

-->

## Cloud computing backend

The development of CPSR was accelerated by using the Arvados
Cloud platform. Arvados is an open source platform for managing,
processing, and sharing genomic and other large scientific and
biomedical data. The Arvados instance was deployed on Amazon AWS
for testing and development and a project was created that
allows for uploading data.

## Sequence uploader

We wrote a Python-based uploader that authenticates with Arvados
using a token. Data gets validated for being a FASTA sequence,
FASTQ raw data and/or metadata in the form of JSON LD that gets
validated against a schema. The uploader can be used
from a command line or using a simple web interface.

## Creating a Pangenome

### FASTA to GFA workflow

The first workflow (1) we implemented was a FASTA to Graphical
Fragment Assembly (GFA) Format conversion. When someone uploads a
sequence in FASTA format it gets combined with all known viral
sequences in our storage to generate a pangenome or variation
graph (VG). The full pangenome is made available as a
downloadable GFA file together with a visualisation (Figure 1).

### FASTQ to GFA workflow

In the next step we introduced a workflow (2) that takes raw
sequence data in fastq format and converts that into FASTA.
This FASTA file, in turn, gets fed to workflow (1) to generate
the pangenome.

## Creating linked data workflow

We created a workflow (3) that takes GFA and turns that into
RDF. Together with the metadata at upload time a single RDF
resource is compiled that can be linked against external
resources such as Uniprot and Wikidata. The generated RDF file
can be hosted in any triple store and queried using SPARQL.

## Creating a Phylogeny workflow

WIP

## Other workflows?

# Discussion

CPSR is a data repository with computational pipelines that will
persist during pandemics.  Unlike other data repositories for
Sars-COV-2 we created a repository that immediately computes the
pangenome of all available data and presents that in useful
formats for futher analysis, including visualisations, GFA and
RDF. Code and data are available and written using best practises
and state-of-the-art standards. CPSR can be deployed by anyone,
anywhere.

CPSR is designed to abide by FAIR data principles (expand...)

CPSR is primed with viral data coming from repositories that have
no sharing restrictions. The metadata includes relevant
attribution to uploaders. Some institutes have already committed
to uploading their data to CPSR first so as to warrant sharing
for computation.

CPSR is currently running on an Arvados cluster in the cloud. To
ascertain the service remains running we will source money from
project during pandemics. The workflows are written in CWL which
means they can be deployed on any infrastructure that runs
CWL. One of the advantages of the CC-4.0 license is that we make
available all uploaded sequence and meta data, as well as
results, online to anyone. So the data can be mirrored by any
party. This guarantees the data will live on.

<!-- Future work... -->

We aim to add more workflows to CPSR, for example to prepare
sequence data for submitting in other public repositories, such
as EBI ENA and GISAID. This will allow researchers to share data
in multiple systems without pain, circumventing current sharing
restrictions.

# Acknowledgements

We thank the COVID-19 BioHackathon 2020 and ELIXIR for creating a
unique event that triggered many collaborations. We thank Curii
Corporation for their financial support for creating and running
Arvados instances.  We thank Amazon AWS for their financial
support to run COVID-19 workflows. We also want to thank the
other working groups in the BioHackathon who generously
contributed onthologies, workflows and software.


# References
