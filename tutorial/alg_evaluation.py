# map edges and calculate roc score
from sklearn import metrics
import seaborn as sns
import numpy as np


def mapping_edges(df_1, df_2, df_1_col_1, df_1_col_2, df_2_col_1, df_2_col_2):
    df_1['tmp1'] = df_1[df_1_col_1] + '_' + df_1[df_1_col_2]
    df_2['tmp1'] = df_2[df_2_col_1] + '_' + df_2[df_2_col_2]

    df_2['tmp2'] = df_2[df_2_col_2] + '_' + df_2[df_2_col_1]

    return len(set(df_1['tmp1']) & set(df_2['tmp1'])) + len(set(df_1['tmp1']) & set(df_2['tmp2']))


def evaluation(links, Full_Ref):
    Detected = links.shape[0]
    Ref_links = Full_Ref[Full_Ref.weight == 1]
    TP = mapping_edges(links, Ref_links, 'source', 'target', 'source', 'target')
    FN = Ref_links.shape[0] - TP
    FP = Detected - TP
    TN = Full_Ref.shape[0] - Ref_links.shape[0] - Detected + TP

    Precision = TP / (TP + FP)
    Recall = TP / (TP + FN)
    FDR = FP / (TN + FP)
    F1_Score = (2 * Precision * Recall) / (Precision + Recall)

    # print('TP:', TP, '\n', 'FN:', FN, '\n', 'FP:', FP, '\n', 'TN:', TN)

    links['name'] = links['source'] + '_' + links['target']

    for name in links.name:
        tmp = list(name.split('_'))
        if name in list(Full_Ref['name']):
            Full_Ref.loc[Full_Ref.name == name, 'weight_2'] = float(
                format(float(links.loc[links.name == name, 'weight']), '.2f'))

        elif tmp[1] + '_' + tmp[0] in list(Full_Ref['name']):
            Full_Ref.loc[Full_Ref.name == tmp[1] + '_' + tmp[0], 'weight_2'] = float(
                format(float(links.loc[links.name == name, 'weight']), '.2f'))

        else:
            continue

    Full_Ref = Full_Ref.fillna(0)
    auc = metrics.roc_auc_score(np.array(Full_Ref['weight'].abs()), np.array(Full_Ref['weight_2'].abs()))
    fpr, tpr, threshold_1 = metrics.roc_curve(Full_Ref['weight'].abs(), Full_Ref['weight_2'].abs())
    pre, recall, threshold_2 = metrics.precision_recall_curve(Full_Ref['weight'].abs(),
                                                              Full_Ref['weight_2'].abs())

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

