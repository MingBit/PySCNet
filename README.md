# PySCNet: A tool for reconstructing and analyzing gene regulatory network from single-cell RNA-Seq data
There are four modules:
1) **Pro-precessing**: initialize a gnetData object consisting of Expression Matrix, Cell Attributes, Gene Attributes and Network Attributes;
2) **BuildNet**: reconstruct GRNs by various methods implemented in docker;
3) **NetEnrich**: network analysis including consensus network detection, gene module identification and trigger path prediction as well as network fusion;
4) **Visulization**: network illustration.

![Overview](https://github.com/MingBit/PySCNet/blob/master/Overview.png)

# TO-DO:
1) Add an Auto-ML based pipeline to Pre-Processing module;
2) Collect more GRN methods to BuildNet module;
3) Update other network-based algorithms to NetEnrich module;
4) Build R-Shiny app for visualization
5) Test with integrated sc RNA-seq data.