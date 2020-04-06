#!/usr/bin/env python3
import os
import sys

import setuptools.command.egg_info as egg_info_cmd
from setuptools import setup

SETUP_DIR = os.path.dirname(__file__)
README = os.path.join(SETUP_DIR, "README.rst")

try:
    import gittaggers

    tagger = gittaggers.EggInfoFromGit
except ImportError:
    tagger = egg_info_cmd.egg_info

install_requires = []

needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest < 6", "pytest-runner < 5"] if needs_pytest else []

setup(
    name="bh20-seq-uploader",
    version="1.0",
    description="Biohackathon sequence uploader",
    long_description=open(README).read(),
    long_description_content_type="text/x-rst",
    author="Peter Amstutz",
    author_email="peter.amstutz@curii.com",
    license="Apache 2.0",
    packages=["bh20sequploader"],
    install_requires=install_requires,
    setup_requires=[] + pytest_runner,
    tests_require=["pytest<5"],
    entry_points={
        "console_scripts": [
            "bh20sequploader=bh20sequploader.main:main"
        ]
    },
    zip_safe=True,
    cmdclass={"egg_info": tagger},
    python_requires=">=3.5, <4",
)
