import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "psfam",
    version = "0.0.15",
    author = "Ben Reggio",
    author_email = "benjreggio@gmail.com",
    description = ("Pauli string organizer for dense operators."),
    license = "MIT",
    packages=['psfam'],
    long_description=read('README.md'),
    long_description_content_type = "text/markdown",
    install_requires = ['scikit-build>= 0.1','markdown','numpy>=1.22','qiskit>=0.35','galois>=0.1','scipy>=1.8'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
    ],
)