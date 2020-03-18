#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 12:09:41 2019

@author: mwu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:51:43 2019

@author: mwu
"""


from PySCNet.BuildNet import gne_dockercaller as gdocker
from PySCNet.BuildNet.Models import model_node2vec


def call_node2vec(gnetdata, method, p, q, dim_list, walk_list, num_walks_list, workers = 8, cell_clusterid = None):
    
    if cell_clusterid is not None:

         cell_info = gnetdata.CellAttrs
         Expr = gnetdata.ExpMatrix[list(cell_info.loc[cell_info.cluster_id.isin([cell_clusterid])].index)]
    else:
         Expr = gnetdata.ExpMatrix
    
    for dimensions in dim_list:
        for walk_length in walk_list:
            for num_walks in num_walks_list:
                
                links = model_node2vec.run_node2vec(Expr, method = method, p = p, q = q, walk_len=walk_length, 
                                            num_walks=num_walks, size=dimensions, workers = workers)
    
#                links = gdocker._remove_duplicate(links)
    gnetdata._add_netattr('links', links)
    
    return(gnetdata)
                
    
