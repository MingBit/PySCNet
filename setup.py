import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PySCNet",
    version="0.0.1",
    author="Ming Wu",
    author_email="ming.wu@tum.de",
    description="A tool for reconstructing and analyzing gene regulatory network from single-cell RNA-Seq data",
    url="https://github.com/MingBit/PySCNet",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)