
import pandas as pd
import os
import itertools

def run_scode(Expr):
    
    Expr.to_csv('Expr_1.txt', sep='\t', header=False, index=False)
#        Rscript SCODE.R data/Expr.txt data/time_2.txt out/ 114 4 241 114
    cmd = 'ruby run_R.rb Expr_1.txt time_point.txt out ' + str(Expr.shape[0]) + ' 4 ' + str(Expr.shape[1]) + ' ' + str(Expr.shape[0]) + ' 10' 
    os.system(cmd)
    

def convert_to_links(Expr):
    
    res_A = pd.read_csv('out/meanA.txt', sep = '\t', header = -1).abs()
    res_A.index = Expr.index
    res_A.columns = Expr.index
    
    
    comb = list(itertools.combinations(list(res_A.index), 2))
    link = pd.DataFrame(comb, columns = ['source', 'target'])
    weights = list()
    for i in range(len(link)):
       weights.append((res_A[link['source'][i]][link['target'][i]] + res_A[link['target'][i]][link['source'][i]]) / 2)
    
    link['weight'] = weights
    link.to_csv('links.txt', sep = '\t', index = False)
    

Expr = pd.read_csv('Expr.txt', sep = '\t', header = 0, index_col = 0)
run_scode(Expr)
convert_to_links(Expr)
    
    
