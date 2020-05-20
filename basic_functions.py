from __future__ import absolute_import
import pandas as pd
import numpy as np
import itertools
import copy
import seaborn as sns


def read_data(file_path, dropRate=0):
    Sim = pd.read_csv(file_path + '/ExpressionData.csv', sep=',', index_col=0)
    Sim_Ref_raw = pd.read_csv(file_path + '/refNetwork.csv', sep=',', index_col=None)

    Sim_Ref_raw['weight'] = Sim_Ref_raw['Type'].map({'+': 1, '-': -1})
    Sim_Ref_raw = Sim_Ref_raw.drop(columns='Type')
    Sim_Ref_raw.columns = ['source', 'target', 'weight']

    Full_Ref = pd.DataFrame(itertools.permutations(Sim.index, 2), columns=['source', 'target'])
    Sim_Ref = pd.merge(Full_Ref, Sim_Ref_raw, how='outer').fillna(0)

    Sim_Ref_dropout = copy.deepcopy(Sim_Ref)
    pos = list(Sim_Ref_dropout[abs(Sim_Ref_dropout.weight) == 1].index)
    replaceNum = np.random.choice(pos, int(len(pos) * dropRate), replace=False)
    Sim_Ref_dropout.loc[replaceNum, 'weight'] = 0

    return [Sim, Sim_Ref, Sim_Ref_dropout, Sim_Ref_raw.shape[0]]


def build_curves(ax, node_dict_list, GENIE3_dict, PIDC_dict, CORR_dict,
                 curve, filename, p, q, dropRate):

    keywords = ['fpr', 'tpr'] if curve == 'ROC' else ['recall_list', 'pre']
    fill = 1 if curve == 'ROC' else 0
    node_fpr_pre_df = pd.DataFrame(node_dict_list[i][keywords[0]] for i in range(len(node_dict_list))).fillna(fill)
    node_tpr_recall_df = pd.DataFrame(node_dict_list[i][keywords[1]] for i in range(len(node_dict_list))).fillna(fill)
    node_auc_list = list(node_dict_list[i]['auc'] for i in range(len(node_dict_list)))
    node_avgpre_list = list(node_dict_list[i]['avg_pre'] for i in range(len(node_dict_list)))

    colors = sns.color_palette().as_hex() + sns.color_palette('hls', 8).as_hex()

    for i in range(len(node_dict_list)):
        ax.plot(node_dict_list[i][keywords[0]], node_dict_list[i][keywords[1]], color='grey', alpha=0.2)

    fpr_dict = {'GENIE3': GENIE3_dict[keywords[0]], 'PIDC': PIDC_dict[keywords[0]], 'CORR': CORR_dict[keywords[0]],
                'Node2Vec': node_fpr_pre_df.mean(0)}
    tpr_dict = {'GENIE3': GENIE3_dict[keywords[1]], 'PIDC': PIDC_dict[keywords[1]], 'CORR': CORR_dict[keywords[1]],
                'Node2Vec': node_tpr_recall_df.mean(0)}
    auc_dict = {'GENIE3': GENIE3_dict['auc'], 'PIDC': PIDC_dict['auc'], 'CORR': CORR_dict['auc'],
                'Node2Vec': np.mean(node_auc_list)}
    avgpre_dict = {'GENIE3': GENIE3_dict['avg_pre'], 'PIDC': PIDC_dict['avg_pre'], 'CORR': CORR_dict['avg_pre'],
                   'Node2Vec': np.mean(node_avgpre_list)}

    for i, color in zip(list(auc_dict.keys()), colors):
        ax.plot(fpr_dict[i], tpr_dict[i],
                label='ROC curve {0} (area = {1:0.2f})'.format(i, auc_dict[i]) if curve == 'ROC' else
                'PR curve {0} (area = {1:0.2f})'.format(i, avgpre_dict[i]),
                color=color)
    if curve == 'ROC':
        ax.plot([0, 1], [0, 1], 'k--')

    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.0])
    ax.set_xlabel('False Positive Rate' if curve == 'ROC' else 'Recall', fontsize=20)
    ax.set_ylabel('True Positive Rate' if curve == 'ROC' else 'Precision', fontsize=20)
    ax.set_title(filename + '_p_' + str(p) + '_q_' + str(q) + '_dropRate_' + str(dropRate), fontsize=12)
    ax.legend(loc="lower right")


def build_plot(ax, node_dict_list, GENIE3_dict, PIDC_dict, CORR_dict):
    node_precision_list = list(node_dict_list[i]['precision'] for i in range(len(node_dict_list)))
    node_recall_list = list(node_dict_list[i]['recall'] for i in range(len(node_dict_list)))
    node_f1score_list = list(node_dict_list[i]['f1_score'] for i in range(len(node_dict_list)))

    Algorithms = ['GENIE3', 'PIDC', 'CORR', 'Node2Vec']
    Eva_Methods = ['Precision', 'Recall', 'F1_Score']

    data = [[GENIE3_dict['precision'], PIDC_dict['precision'], CORR_dict['precision'],
             np.mean(node_precision_list)],
            [GENIE3_dict['recall'], PIDC_dict['recall'], CORR_dict['recall'],
             np.mean(node_recall_list)],
            [GENIE3_dict['f1_score'], PIDC_dict['f1_score'], CORR_dict['f1_score'],
             np.mean(node_f1score_list)]]

    colors = sns.color_palette().as_hex() + sns.color_palette('hls', 8).as_hex()
    X = np.arange(4)
    for i in range(len(data)):
        print(data[i])
        ax.bar(X + 0.25 * i, data[i], color=colors[i], width=0.25)

    ax.set_xlabel('Algorithms', fontsize=25)
    ax.set_ylabel('Performance', fontsize=25)
    ax.set_xticks(X + 0.25)
    ax.set_xticklabels(Algorithms, fontsize=20)
    ax.legend(Eva_Methods, loc='upper left')
