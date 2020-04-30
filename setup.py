#!/usr/bin/env python3
import os
import sys

import setuptools.command.egg_info as egg_info_cmd
from setuptools import setup

SETUP_DIR = os.path.dirname(__file__)
README = os.path.join(SETUP_DIR, "README.md")

try:
    import gittaggers

    tagger = gittaggers.EggInfoFromGit
except ImportError:
    tagger = egg_info_cmd.egg_info

install_requires = ["arvados-python-client", "schema-salad",
                    "python-magic", "pyshex", "py-dateutil"]
web_requires = ["flask", "pyyaml"]

needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest < 6", "pytest-runner < 5"] if needs_pytest else []

setup(
    name="bh20-seq-uploader",
    version="1.0",
    description="Biohackathon sequence uploader",
    long_description=open(README).read(),
    long_description_content_type="text/markdown",
    author="Peter Amstutz",
    author_email="peter.amstutz@curii.com",
    license="Apache 2.0",
    packages=["bh20sequploader", "bh20seqanalyzer", "bh20simplewebuploader"],
    package_data={"bh20sequploader": ["bh20seq-schema.yml",
                                      "bh20seq-options.yml",
                                      "bh20seq-shex.rdf",
                                      "validation/formats",
                                      "SARS-CoV-2-reference.fasta",],
    },
    install_requires=install_requires,
    extras_require={
        'web': web_requires
    },
    setup_requires=[] + pytest_runner,
    tests_require=["pytest<5"],
    entry_points={
        "console_scripts": [
            "bh20-seq-uploader=bh20sequploader.main:main",
            "bh20-seq-analyzer=bh20seqanalyzer.main:main"
        ]
    },
    zip_safe=True,
    cmdclass={"egg_info": tagger},
    python_requires=">=3.5, <4",
)
