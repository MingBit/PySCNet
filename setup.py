
import setuptools

with open("requirements.txt", "r") as f:
    required = f.read().splitlines()

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyscnet",
    version="0.0.3",
    author="Ming Wu",
    license='MIT',
    author_email="ming.wu@tum.de",
    description="A tool for reconstructing and analyzing gene regulatory network from single-cell RNA-Seq data",
    url="https://github.com/MingBit/PySCNet",
    packages=setuptools.find_packages(),
    install_requires=required,
    include_package_data=True,
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)
