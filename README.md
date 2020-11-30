[![Codacy Badge](https://api.codacy.com/project/badge/Grade/d3c17aac77e14f6bb17b33f875ff7471)](https://app.codacy.com/manual/MingBit/PySCNet?utm_source=github.com&utm_medium=referral&utm_content=MingBit/PySCNet&utm_campaign=Badge_Grade_Dashboard)
[![License](https://img.shields.io/github/license/MingBit/PySCNet)](https://github.com/MingBit/PySCNet/blob/master/LICENSE)
[![Build Status](https://travis-ci.org/MingBit/PySCNet.svg?branch=master)](https://travis-ci.org/MingBit/PySCNet)
[![Documentation Status](https://readthedocs.org/projects/pyscnet/badge/?version=latest)](https://pyscnet.readthedocs.io/en/latest/?badge=latest)

# PySCNet: A tool for reconstructing and analyzing gene regulatory network from single-cell RNA-Seq data
There are four modules:
1) **Pro-precessing**: initialize a gnetData object consisting of Expression Matrix, Cell Attributes, Gene Attributes and Network Attributes;
2) **BuildNet**: reconstruct GRNs by various methods implemented in docker;
3) **NetEnrich**: network analysis including consensus network detection, gene module identification and trigger path prediction as well as network fusion;
4) **Visulization**: network illustration.

![Overview](https://github.com/MingBit/PySCNet/blob/master/images/workflow_update.eps)

# :tada: :confetti_ball: Create your own GRNs
[Dashboard](https://github.com/MingBit/PySCNet/blob/master/images/pyscnet_dashboard.gif) is available now for creating your own GRNs.
Cell specific GRNs and network analysis results can be saved as a pickle object and upload onto PySCNet-Dashboard.
It provides parameter settings and allows for parameter adjustment and GRNs customization. 

# Installation
1) clone from github:
`git clone https://github.com/MingBit/PySCNet`
2) create a new folder and set up:
`mkdir dist | python setup.py sdist`
3) install pyscnet:
`pip install dist/pyscnet-0.0.3.tar.gz`

# Tutorial
Mouse HSC data preprocessed and analyzed by stream as explained in this 
[tutorial](https://github.com/MingBit/PySCNet/blob/master/tutorial/pyscnet_stream.ipynb). 

# TO-DO
1) Add an Auto-ML based pipeline to Pre-Processing module;
2) Collect more GRN methods to BuildNet module;
3) Update other network-based algorithms to NetEnrich module;
4) Build R-Shiny app for visualization
5) Test with integrated sc RNA-seq data.

# Cite
- :smile_cat: This work has been presented at [ISMB/ECCB 2019](https://www.iscb.org/ismbeccb2019);
- :paw_prints: Go to [my poster](https://f1000research.com/posters/8-1359);
- :page_with_curl: Reference: *Wu M, Kanev K, Roelli P and Zehn D. PySCNet:
A tool for reconstructing and analyzing gene regulatory network from single-cell RNA-Seq data [version 1; not peer reviewed]. F1000Research 2019, 8(ISCB Comm J):1359 (poster) (doi: 10.7490/f1000research.1117280.1)*
