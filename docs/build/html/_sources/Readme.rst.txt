Introduction
=============================================================================================================
PySCNet: A tool for reconstructing and analyzing gene regulatory network from single-cell RNA-Seq data

Modules
----------------

There are four modules:

1)  **Pro-precessing**: initialize a gnetData object consisting of Expression Matrix, Cell Attributes, Gene Attributes and Network Attributes;
2) **BuildNet**: reconstruct GRNs by various methods implemented in docker;
3) **NetEnrich**: network analysis including consensus network detection, gene module identification and trigger path prediction as well as network fusion;
4) **Visulization**: network illustration.

.. figure:: ../../images/Overview.png
    :width: 600

Features
----------------
`Shinyapp`_ is available now for creating your own GRNs.
Once the cells are grouped into several clusters and linkage tables are generated for each/all clusters, you can export the results
as pickle object and uplaod onto Shinyapp. Cell attributes, Gene attributes and Network attributes are illustrated here.
As shown belows, you can set your own thresholds to build each/all cluster-specific GRNs.

.. figure:: ../../images/ShinyApp.gif
    :width: 600

Cite
----------------
- This work has been presented at `ISMB/ECCB 2019`_
- Go to `my poster`_

.. _my poster: https://f1000research.com/posters/8-1359
.. _ISMB/ECCB 2019: https://www.iscb.org/ismbeccb2019
.. _Shinyapp: https://github.com/MingBit/PySCNet/blob/master/ShinyApp.gif