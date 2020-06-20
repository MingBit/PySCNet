import os
import setuptools

with open('requirements.txt') as f:
    required = f.read().splitlines()

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyscnet",
    version="0.0.1",
    author="Ming Wu",
    license='MIT',
    author_email="ming.wu@tum.de",
    description="A tool for reconstructing and analyzing gene regulatory network from single-cell RNA-Seq data",
    url="https://github.com/MingBit/PySCNet",
    download_url="https://github.com/MingBit/PySCNet/archive/v0.0.1.tar.gz",
    packages=setuptools.find_packages(),
    install_requires=required,
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)
