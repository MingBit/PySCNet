# map edges and calculate roc score

from sklearn import metrics
import seaborn as sns
import numpy as np


def mapping_edges(df_1_link, df_2_link):

    df_1_link['edge'] = df_1_link[['source', 'target']].apply(lambda x: '_'.join(sorted(x)), axis=1)
    df_2_link['edge'] = df_2_link[['source', 'target']].apply(lambda x: '_'.join(sorted(x)), axis=1)

    return len(set(df_1_link['edge']).intersection(df_2_link['edge']))

def evaluation(links, Full_Ref):
    
    ##ToDo: Test
    
    Detected = links.shape[0]
    Ref_links = Full_Ref[Full_Ref.weight == 1].reset_index(drop=True)
    TP = mapping_edges(links, Ref_links)
    FN = Ref_links.shape[0] - TP
    FP = Detected - TP
    TN = Full_Ref.shape[0] - Ref_links.shape[0] - Detected + TP

    Precision = TP / (TP + FP)
    Recall = TP / (TP + FN)
    FDR = FP / (TN + FP)
    F1_Score = (2 * Precision * Recall) / (Precision + Recall)

    links['name'] = links['source'] + '_' + links['target']

    name_exists = Full_Ref['name'].isin(links['name'])
    tmp = links.loc[links['name'].isin(Full_Ref.loc[name_exists, 'name']), 'weight']
    Full_Ref.loc[name_exists, 'weight_2'] = tmp.astype(float).apply(lambda x: format(x, '.2f'))

    name_exists_reversed = Full_Ref['name'].str.split('_').apply(lambda x: x[1] + '_' + x[0]).isin(links['name'])
    tmp_reversed = links.loc[links['name'].isin(Full_Ref.loc[name_exists_reversed, 'name']), 'weight']
    Full_Ref.loc[name_exists_reversed, 'weight_2'] = tmp_reversed.astype(float).apply(lambda x: format(x, '.2f'))

    Full_Ref['weight_2'] = Full_Ref['weight_2'].fillna(0)
    auc = metrics.roc_auc_score(np.abs(Full_Ref['weight']), np.abs(Full_Ref['weight_2']))
    fpr, tpr, _ = metrics.roc_curve(np.abs(Full_Ref['weight']), np.abs(Full_Ref['weight_2']))
    pre, recall, _ = metrics.precision_recall_curve(np.abs(Full_Ref['weight']), np.abs(Full_Ref['weight_2']))

    avg_pre_auc = metrics.auc(recall, pre)

    return {'fpr': fpr, 'tpr': tpr, 'pre': pre, 'recall': recall, 'auc': auc,
            'avg_pre': avg_pre_auc, 'Precision': Precision, 'Recall': Recall, 'F1_Score': F1_Score}


# Create ROC and PR curves
def build_curves(ax, dict_of_algs, curve, filename, colors=None, **kwarg):
    names = dict_of_algs.keys()

    keywords = ['fpr', 'tpr'] if curve == 'ROC' else ['recall', 'pre']
    fill = 1 if curve == 'ROC' else 0

    colors = sns.color_palette().as_hex() + sns.color_palette('hls', 8).as_hex() if colors is None else colors
    section_dict = {sec: {} for sec in ['fpr_dict', 'tpr_dict', 'auc_dict', 'avgpre_dict']}

    for name in names:
        
        section_dict['fpr_dict'][name] = dict_of_algs[name][keywords[0]]
        section_dict['tpr_dict'][name] = dict_of_algs[name][keywords[1]]
        section_dict['auc_dict'][name] = dict_of_algs[name]['auc']
        section_dict['avgpre_dict'][name] = dict_of_algs[name]['avg_pre']

    for i, color in zip(list(section_dict['auc_dict'].keys()), colors):
        ax.plot(section_dict['fpr_dict'][i], section_dict['tpr_dict'][i],
                label='ROC curve {0} (area = {1:0.2f})'.format(i, section_dict['auc_dict'][i]) if curve == 'ROC' else
                'PR curve {0} (area = {1:0.2f})'.format(i, section_dict['avgpre_dict'][i]),
                color=color, **kwarg)

    if curve == 'ROC':
        ax.plot([0, 1], [0, 1], 'k--')

    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.0])
    ax.set_xlabel('False Positive Rate' if curve == 'ROC' else 'Recall', fontsize=20)
    ax.set_ylabel('True Positive Rate' if curve == 'ROC' else 'Precision', fontsize=20)
    ax.legend(loc="lower right")

