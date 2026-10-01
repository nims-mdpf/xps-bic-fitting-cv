# ------------------------------------------------------------
# Copyright (c) 2024, National Institute for Materials Science
#
# This software is released under the MIT License.
#
# Contributor:
#     Hiroshi Shinotsuka
# ------------------------------------------------------------
# coding: utf-8

# Active Shirleyの出力ファイルからAIC, BICを計算

import math
import subprocess

import numpy as np
import pandas as pd
from scipy.special import voigt_profile
from scipy.optimize import fmin

from glob import glob


def voigt_fwhm_approx(sigma, gamma):
    """
    Calculate the approximated Full Width at Half Maximum (FWHM) of a Voigt profile.

    Parameters:
    - sigma (float): The standard deviation (Gaussian component) of the Voigt profile.
    - gamma (float): The half-width at half-maximum (Lorentzian component) of the Voigt profile.

    Returns:
    float: The approximated FWHM of the Voigt profile.
    """
    fwhm_gauss = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sigma
    fwhm_lorentz = 2.0 * gamma
    fwhm_approx = np.sqrt(fwhm_gauss**2 + fwhm_lorentz**2)
    return fwhm_approx


def voigt_fwhm(sigma, gamma):
    """
    Calculate the Full Width at Half Maximum (FWHM) of a Voigt profile.

    Parameters:
    - sigma (float): The standard deviation (Gaussian component) of the Voigt profile.
    - gamma (float): The half-width at half-maximum (Lorentzian component) of the Voigt profile.

    Returns:
    float: The FWHM of the Voigt profile.
    """

    if sigma == 0 and gamma == 0:
        return 0.0

    # Voigt関数を定義
    def voigt_localfunction(x):
        return voigt_profile(x, sigma, gamma)

    # 半値幅を求めるための最適化関数
    def half_max_width(x, halfheight):
        return abs(voigt_localfunction(x) - halfheight)

    # 最適化を行い、Voigt関数の半値幅を取得
    halfheight = voigt_localfunction(0) / 2
    decimal = 5
    # 半値半幅の初期値
    # x0 = 0.0
    x0 = voigt_fwhm_approx(sigma, gamma) / 2
    result = fmin(half_max_width, x0, args=(halfheight,), xtol=10 ** (-decimal), disp=False, full_output=True)
    # print(result)
    hwhm = result[0][0]
    fwhm = np.round(hwhm * 2.0, decimal)
    return fwhm


def calc_AIC_BIC(fitting_result_file, fitting_parameter_file, idHigherS=-1, idLowerS=-1):
    # 高BEと低BEのインデックスを読み込む
    df = pd.read_csv(fitting_result_file, delimiter="=", header=None, nrows=2)
    idHigherOrg = df[1][0]
    idLowerOrg = df[1][1]

    # 入力スペクトルおよびフィッティングスペクトル読み込み
    df = pd.read_csv(fitting_result_file, comment="#")

    # 解析範囲以外を削除
    if idHigherS >= 0 and idLowerS >= 0:
        idHigher = idHigherS
        idLower = idLowerS
    else:
        # df1 = df[np.isnan(df.background) == False]
        df1 = df[pd.isnull(df.background) == False]

        df2 = df1
        if type(df1.background[0]) == object:
            df2 = df1[df1.background != "-nan(ind)"]
        else:
            df2 = df1

        numData = len(df2.index)
        if numData == 0:
            print("Error. All background data is NaN.")
            return

        # idHigher, idLower を読み込み、データ数を取得
        mycmd = 'grep "# idHigher" %s | cut -d" " -f4' % (fitting_result_file)
        idHigher = subprocess.call(mycmd, shell=True)
        mycmd = 'grep "# idLower" %s | cut -d" " -f4' % (fitting_result_file)
        idLower = subprocess.call(mycmd, shell=True)

    numData = idLower - idHigher + 1

    # 始点と終点のエネルギー
    energyHigher = df["Binding_Energy(eV)"][idHigher]
    energyLower = df["Binding_Energy(eV)"][idLower]

    energyHigherOrg = df["Binding_Energy(eV)"][idHigherOrg]
    energyLowerOrg = df["Binding_Energy(eV)"][idLowerOrg]

    df1 = df.iloc[idHigher : (idLower + 1), :]

    # バックグラウンドに負の値があれば評価から外す
    if np.any(df1.background < 0):
        print("Error. Background has negative value.")
        return

    # フィッティングパラメータ読み込み（ピーク本数を取得）
    dfp = pd.read_csv(fitting_parameter_file)

    # ピーク本数およびパラメータ数
    numPeak = len(dfp.index)
    numParam = 4 * numPeak + 2

    # 最小半値半幅
    # min_HWHM = dfp.width[dfp.width.idxmin()]
    dfp["fwhm"] = dfp.apply(lambda row: voigt_fwhm(row["sigma"], row["gamma"]), axis=1)
    min_HWHM = dfp["fwhm"].min() / 2

    # 一様なガウスノイズを仮定して、ノイズの分散および逸脱度を計算
    variance_gauss = ((df1.spectrum - df1.fitting) ** 2).mean()
    deviance_gauss = numData * (1.0 + np.log(2.0 * np.pi * variance_gauss))

    # そのときのAIC, BIC
    aic_gauss = deviance_gauss + 2.0 * numParam
    bic_gauss = deviance_gauss + np.log(numData) * numParam

    # 強度に依存したガウスノイズを仮定して，分散の最尤推定。
    # スケーリング因子を計算。逸脱度を計算
    param_b = ((df1.spectrum - df1.fitting) ** 2 / df1.fitting).mean()
    scaling_pgauss = 1.0 / param_b
    #    deviance_pgauss = (np.log(df1.fitting)).sum() + numData * (1.0 + np.log(2.0 * np.pi * scaling_pgauss))
    # そのときのAIC, BIC
    #    aic_pgauss = deviance_pgauss + 2.0 * numParam
    #    bic_pgauss = deviance_pgauss + np.log(numData) * numParam

    # 強度に依存したガウスノイズを仮定して，分散の最尤推定。
    # スケーリングを行うとデータそのものが変わってしまうため，AIC, BIC計算には不適切
    # スケーリングは行わずに逸脱度を計算
    deviance_pgauss = (np.log(df1.fitting)).sum() + numData * (1.0 + np.log(2.0 * np.pi * param_b))
    # そのときのAIC, BIC
    aic_pgauss = deviance_pgauss + 2.0 * numParam
    bic_pgauss = deviance_pgauss + np.log(numData) * numParam

    # ポアソンノイズを仮定して、逸脱度を計算
    # 縦軸の単位がカウント数の場合にのみ意味がある
    # 対数ガンマ関数は配列計算できないので tmp 使用
    tmp = [math.lgamma(y + 1.0) for y in df1.spectrum]
    deviance_poisson = 2.0 * sum(df1.fitting - df1.spectrum * np.log(df1.fitting) + tmp)

    # そのときのAIC, BIC
    aic_poisson = deviance_poisson + 2.0 * numParam
    bic_poisson = deviance_poisson + np.log(numData) * numParam

    # カイ二乗
    chi2 = sum((df1.spectrum - df1.fitting) ** 2 / df1.fitting)

    # Reduced カイ二乗
    reduced_chi2 = chi2 / (numData - numParam)

    return [
        numData,
        numPeak,
        numParam,
        variance_gauss,
        deviance_gauss,
        aic_gauss,
        bic_gauss,
        deviance_poisson,
        aic_poisson,
        bic_poisson,
        scaling_pgauss,
        deviance_pgauss,
        aic_pgauss,
        bic_pgauss,
        chi2,
        reduced_chi2,
        min_HWHM,
        idHigher,
        idLower,
        energyHigher,
        energyLower,
        idHigherOrg,
        idLowerOrg,
        energyHigherOrg,
        energyLowerOrg,
    ]


def get_settings_as_dict(globkey="./*/settings.inp"):
    dict_settings = {}
    for f in glob(globkey):
        df_settings = pd.read_csv(f, sep="\s+", comment="#", header=None, usecols=[0, 1], names=["key", "value"])
        dict_settings = df_settings.set_index("key")["value"].to_dict()
    # print(dict_settings)
    return dict_settings


def get_value_from_settings(key, default_value=None, globkey="./*/settings.inp"):
    value = default_value
    dict_settings = get_settings_as_dict(globkey=globkey)
    if key in dict_settings.keys():
        value = dict_settings[key]
    return value
