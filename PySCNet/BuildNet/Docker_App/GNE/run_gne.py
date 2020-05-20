from __future__ import absolute_import
import os
import sys
import numpy as np
sys.path.append(os.getenv("HOME") + '/GNE/')

from sklearn.linear_model import LogisticRegression

import LoadData as data
from evaluation import *
from GNE import GNE
from convertdata import *
import networkx as nx
from sklearn import metrics
import matplotlib.pyplot as plt
import seaborn as sns

#################################### Define dataset and files ##################################################
# Define path
input_path = os.getenv("HOME") + '/MING_V9T/PhD_Pro/Test/Simulation/BoolODE_Data/'
output_path = '/home/mwu/MING_V9T/PhD_Pro/Test/Simulation/GNE_BoolODE_Res/'

filename = sys.argv[1]
if not os.path.exists(output_path + filename):
    os.makedirs(output_path + filename)

new_path = output_path + filename + '/'   
raw_Expr = pd.read_csv(input_path + filename + '/ExpressionData.csv', index_col = 0)

# adj = np.array(raw_Expr.T.corr())
# adj[:]=np.where(adj < 0.5, 0, 1)
# np.fill_diagonal(adj, 0)

refNet = pd.read_csv(input_path + filename + '/refNetwork.csv').drop(columns = 'Type')
refNet.columns = ['source', 'target']

G = nx.from_pandas_edgelist(refNet.sample(np.int(refNet.shape[0] * 0.5)))
adj = nx.to_numpy_array(G)


data_standard = np.log(raw_Expr + 1)
data_standard.to_csv(new_path + '/data_standard.txt', index = False, header = False)

raw_Expr.index = [('G' + str(i + 1)) for i in range(raw_Expr.shape[0])]
raw_Expr.T.to_csv(new_path + '/ExpressionData_2.csv', index = False, sep = '\t')



# Define the input to GNE model
feature_file = new_path + '/ExpressionData_2.csv'

################################# Define parameters to train GNE model #######################################
parameters = {}
# Dimension of topological structure embeddings
parameters['id_embedding_size'] = 5
# Dimension of expression data embeddings
parameters['attr_embedding_size'] = 5
# Dimension of final representation after transformation of concatenated topological properties and expression data representation
parameters['representation_size'] = 5
# Importance of gene expression relative to topological properties
parameters['alpha'] = 1

# Number of negative samples for Negative Sampling
parameters['n_neg_samples'] = 10
# Number of epochs to run the model
parameters['epoch'] = 20
# Number of sample to consider in the batch
parameters['batch_size'] = 10
# Learning rate
parameters['learning_rate'] = 0.05

################################# Load network and split to train and test######################################

num_genes = raw_Expr.shape[0]

def run_gne(path, adj, feature_file, parameters):
    
    dataset = create_train_test_split(path + '/', adj, test_size=0.2, validation_size=0.2)
    
    train_edges = dataset['train_pos']
    train_edges_false = dataset['train_neg']
    val_edges = dataset['val_pos']
    val_edges_false = dataset['val_neg']
    test_edges = dataset['test_pos']
    test_edges_false = dataset['test_neg']
        
    ###################### Combine positive and negative interactions for valdiation and test ######################
    # Create validation edges and labels
    validation_edges = np.concatenate([val_edges, val_edges_false])
    val_edge_labels = np.concatenate([np.ones(len(val_edges)), np.zeros(len(val_edges_false))])
    
    # Create test edges and labels
    test_edges_data = np.concatenate([test_edges, test_edges_false])
    test_edge_labels = np.concatenate([np.ones(len(test_edges)), np.zeros(len(test_edges_false))])
    
    ################################################################################################################
    
    
    ################## load interaction and expression data to fit GNE model and learn embeddings ##################
    # load dataset to fit GNE model
    Data = data.LoadData(path, train_links=train_edges, features_file=feature_file)
    
    # Define GNE model with data and parameters
    model = GNE(path, Data, 1, parameters)
    
    # learn embeddings
    embeddings, attr_embeddings = model.train(validation_edges, val_edge_labels)
     
    ################## Create feature matrix and true labels for training and randomize the rows  ##################
    # Train-set edge embeddings
    pos_train_edge_embs = get_edge_embeddings(embeddings, train_edges)
    neg_train_edge_embs = get_edge_embeddings(embeddings, train_edges_false)
    train_edge_embs = np.concatenate([pos_train_edge_embs, neg_train_edge_embs])
    
    # Create train-set edge labels: 1 = real edge, 0 = false edge
    train_edge_labels = np.concatenate([np.ones(len(train_edges)), np.zeros(len(train_edges_false))])
    
    # Randomize train edges and labels
    index = np.random.permutation([i for i in range(len(train_edge_labels))])
    train_data = train_edge_embs[index, :]
    train_labels = train_edge_labels[index]
    
    ################## Train the logistic regression on training data and predict on test dataset ##################
    # Train logistic regression on train-set edge embeddings
    edge_classifier = LogisticRegression(random_state=0)
    edge_classifier.fit(train_data, train_labels)
    
    # Test-set edge embeddings, labels
    pos_test_edge_embs = get_edge_embeddings(embeddings, test_edges)
    neg_test_edge_embs = get_edge_embeddings(embeddings, test_edges_false)
    test_edge_embs = np.concatenate([pos_test_edge_embs, neg_test_edge_embs])
    
    # Randomize test edges and labels
    index = np.random.permutation([i for i in range(len(test_edge_labels))])
    test_data = test_edge_embs[index, :]
    test_labels = test_edge_labels[index]
    
    # Predict the probabilty for test edges by trained classifier
    test_preds = edge_classifier.predict_proba(test_data)[:, 1]
    test_roc = roc_auc_score(test_labels, test_preds)
    test_ap = average_precision_score(test_labels, test_preds)
    
    fpr, tpr, threshold_1 = metrics.roc_curve(test_labels, test_preds)
    pre, recall, threshold_2 = metrics.precision_recall_curve(test_labels, test_preds)
    
    return([fpr, tpr, pre, recall, test_roc, test_ap])

    

# msg = "Alpha: {0:>6}, GNE Test ROC Score: {1:.9f}, GNE Test AP score: {2:.9f}"
# print(msg.format(parameters['alpha'], test_roc, test_ap))

def build_curves(ax, node_dict_list, curve):

    keywords = ['fpr', 'tpr'] if curve == 'ROC' else ['recall', 'pre']
    fill = 1 if curve == 'ROC' else 0
    
    node_fpr_recall_df = pd.DataFrame(node_dict_list[i][keywords[0]] for i in range(len(node_dict_list))).fillna(fill)
    node_tpr_pre_df = pd.DataFrame(node_dict_list[i][keywords[1]] for i in range(len(node_dict_list))).fillna(fill)
    
    node_auc_list = list(node_dict_list[i]['test_roc'] for i in range(len(node_dict_list)))
    node_avgpre_list = list(node_dict_list[i]['test_ap'] for i in range(len(node_dict_list)))
    
    colors = sns.color_palette().as_hex() + sns.color_palette('hls', 8).as_hex()
    
    for i in range(len(node_dict_list)):
        
        ax.plot(node_dict_list[i][keywords[0]], node_dict_list[i][keywords[1]], color = 'grey', alpha = 0.4)    
    
    fpr_dict = {'GNE': node_fpr_recall_df.mean(0)}
    tpr_dict = {'GNE': node_tpr_pre_df.mean(0)}
    auc_dict = {'GNE': np.mean(node_auc_list)}
    avgpre_dict = {'GNE': np.mean(node_avgpre_list)}
    
    for i, color in zip(list(auc_dict.keys()), colors):
            ax.plot(fpr_dict[i], tpr_dict[i], 
                     label='ROC curve {0} (area = {1:0.2f})'.format(i, auc_dict[i]) if curve == 'ROC' else 
                     'PR curve {0} (area = {1:0.2f})'.format(i, avgpre_dict[i]),
                     color = color)
    if curve == 'ROC':
        
        ax.plot([0, 1], [0, 1], 'k--') 
        
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.0])
    ax.set_xlabel('False Positive Rate', fontsize = 20)
    ax.set_ylabel('True Positive Rate', fontsize = 20)
    # ax.set_title(filename + '_p_' + str(p) + '_q_' + str(q) + ':Algorithms performance', fontsize = 12)
    ax.legend(loc="lower right")


repeat_dict_list = list()
for i in range(10):
       
    repeat_dict = dict(zip(['fpr', 'tpr', 'pre', 'recall', 'test_roc', 'test_ap'], 
                         run_gne(new_path, adj, feature_file, parameters)))
    repeat_dict_list.append(repeat_dict)
 


fig = plt.figure(figsize = (20,10))
grid = plt.GridSpec(1, 2, wspace=0.2, hspace=0.2)
ax1 = fig.add_subplot(grid[0,0])
ax2 = fig.add_subplot(grid[0,1])
fig.suptitle(filename + ':Algorithms performance', fontsize=30)
build_curves(ax1, repeat_dict_list, 'ROC')
build_curves(ax2, repeat_dict_list, 'PCR')

fig.savefig(newpath + filename + '_GNE.pdf')




# embeddings_file = open(path + "embeddings_trainsize_alpha_"+str(parameters['alpha'])+".pkl", 'wb')
# pickle.dump(embeddings, embeddings_file)
# embeddings_file.close()

################################################################################################################

