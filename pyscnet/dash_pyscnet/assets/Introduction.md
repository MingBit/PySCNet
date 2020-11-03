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

### GnetData object contains the following parts:  
1) **ExpMatrix**: Raw Gene count matrix;  
2) **CellAttrs**: Cell annotation;  
3) **GeneAttrs**: Gene annotation;  
4) **NetAttrs**: multiple unstructured annotation (e.g. linkage table, graph, gene module, gene centrality).   


# Create your own GRNs
[Dashboard](https://github.com/MingBit/PySCNet/blob/master/images/ShinyApp.gif) is available now for creating your own GRNs.
Once the cells are grouped into several clusters and linkage tables are generated for each/all clusters, you can export the results
as pickle object and uplaod onto Shinyapp. Cell attributes, Gene attributes and Network attributes are illustrated here.
As shown belows, you can set your own thresholds to build each/all cluster-specific GRNs.

# Tutorial
PBMC data preprocessed and analyzed by scanpy as explained in this [tutorial](https://github.com/MingBit/PySCNet/blob/master/tutorial/pyscnet_pbmc.ipynb). 

# Cite
-  This work has been presented at [ISMB/ECCB 2019](https://www.iscb.org/ismbeccb2019);
-  Go to [my poster](https://f1000research.com/posters/8-1359);
-  Reference: *Wu M, Kanev K, Roelli P and Zehn D. PySCNet:
A tool for reconstructing and analyzing gene regulatory network from single-cell RNA-Seq data [version 1; not peer reviewed]. F1000Research 2019, 8(ISCB Comm J):1359 (poster) (doi: 10.7490/f1000research.1117280.1)*
