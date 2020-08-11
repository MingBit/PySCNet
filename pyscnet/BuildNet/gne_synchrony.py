#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:51:34 2019

@author: mwu
"""

import numpy as np
import pandas as pd
from scipy.signal import hilbert, butter, filtfilt
from itertools import combinations


def __butter_bandpass(low_cut, high_cut, fs, order=5):
    nyq = 0.5 * fs
    low = low_cut / nyq
    high = high_cut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def __butter_bandpass_filter(data, low_cut, high_cut, fs, order=5):
    b, a = __butter_bandpass(low_cut, high_cut, fs, order=order)
    y = filtfilt(b, a, data, axis=0)
    return y


def __synchrony_cal(df, low_cut, high_cut, fs, order):
    d1 = df.iloc[:, :1].interpolate().values
    d2 = df.iloc[:, :2].interpolate().values
    y1 = __butter_bandpass_filter(d1, low_cut, high_cut, fs, order)
    y2 = __butter_bandpass_filter(d2, low_cut, high_cut, fs, order)

    al1 = np.angle(hilbert(y1), deg=False)
    al2 = np.angle(hilbert(y2), deg=False)
    phase_synchrony = (1 - np.sin(np.abs(al1 - al2) / 2)).mean()

    return phase_synchrony


def __window_rolling(df, window_size):
    single_rank = list()
    for window in window_size:
        rolling_r = df['source'].rolling(window=window, center=True).corr(
            df['target'].rolling(window=window, center=True)).fillna(0)
        single_rank.append(rolling_r)
    single_rank = sum(single_rank)/len(single_rank)
    rank = abs(len(single_rank[single_rank > 0]) - len(single_rank[single_rank < 0]))

    return rank


def get_synchrony(gnetdata, method, cell=None, feature=None, cell_clusterid=None, select_by=None,
                  fs=50., low_cut=2, high_cut=15, order=5, window_size=[50]):
    """
    :param gnetdata: Gnetdata object, default None.
    :param method: str, default None. methods: [window_rolling, phase_synchrony]
    :param cell: list, default None. a list of cell names
    :param feature: list, default None. a list of gene names
    :param cell_clusterid: str, default None. cell with cell_clusterid will be selected
    :param select_by: str, default None. key of filtering cells
    :param fs: float, default 50. the sampling frequency of the digital system
    :param low_cut: int, default 2. low threshold for butter filtering
    :param high_cut: int, default 15. high threshold for butter filtering
    :param order: int, default 5. the order of the butter filter
    :param window_size: list, default [50]. window size for window rolling correlation test
    :return: Gnetdata object with links saved in NetAttrs
    """

    assert method in ['window_rolling', 'phase_synchrony'], 'only window_rolling and phase_synchrony are valid!'

    feature = gnetdata.ExpMatrix.index if feature is None else feature

    if cell_clusterid is None:
        cell = gnetdata.ExpMatrix.columns if cell is None else cell
    else:
        cell_info = gnetdata.CellAttrs['CellInfo']
        cell = list(cell_info.loc[cell_info[select_by].isin([cell_clusterid])].index)

    subExpr = gnetdata.ExpMatrix.loc[feature, cell]
    link = pd.DataFrame(combinations(set(feature), 2), columns=['source', 'target'])
    phase_synchrony = list()

    for i in range(link.shape[0]):
        source = link.loc[i]['source']
        target = link.loc[i]['target']

        a = subExpr.loc[source]
        b = subExpr.loc[target]

        synchrony = __synchrony_cal(pd.DataFrame({'source': a, 'target': b}), low_cut, high_cut,
                                    fs, order) if method == 'phase_synchrony' else \
            __window_rolling(pd.DataFrame({'source': a, 'target': b}), window_size)

        phase_synchrony.append(synchrony)

    link['weight'] = phase_synchrony
    gnetdata._add_netattr(method + '_links', link)
    print(method + '_links added into NetAttrs')

    return gnetdata
