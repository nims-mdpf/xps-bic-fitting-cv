# coding: utf-8
import os
import shutil
import sys
import csv
import subprocess
import time
from glob import glob
import numpy as np
import pandas as pd
import time, datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MaxNLocator

from auto_simplex.myfunc import get_value_from_settings

startTime = time.time()  ## start time

PYTHON = "python3"
bindir = os.path.abspath(os.path.dirname(__file__))
PY_AUTO_SIMPLEX = f"{bindir}/auto_simplex/automatic_xps_peak_separation_single.py"
EXE_GA = f"{bindir}/ga/decompositionMultipleSpectra.exe"
PY_PARAM_SIMPLEX = f"{bindir}/param_simplex/automatic_xps_peak_separation_single_after_GA.py"
PY_PPT = f"{bindir}/make_result_figures_ppt2.py"
K_max = 6


def make_config(fnm, noiseType="-g"):
    data = np.loadtxt(fnm, delimiter=",")
    # print(data.shape)
    ymin = data[:, 1].min()
    ymax = data[:, 1].max()
    intensity_max = ymax - ymin
    # print(f"intensity_max = {intensity_max}")

    fwhm = 1.6
    gw_alpha = 4.0
    gw_beta = fwhm / 4.0 / gw_alpha
    lw_alpha = gw_alpha
    lw_beta = gw_beta
    intensity_alpha = 2.0
    intensity_beta = intensity_max

    bgh_mean = data[:10, 1].mean()
    bgh_std = data[:10, 1].std()
    bgl_mean = data[-10:, 1].mean()
    bgl_std = data[-10:, 1].std()

    # bgh_alpha = 2.0
    # bgh_beta = bgh_mean
    # bgl_alpha = 2.0
    # bgl_beta = bgl_mean
    ## 端点近傍の平均値と標準偏差を使ってガンマ分布の形状母数と尺度母数を求める
    ## 値がそれなりに大きければ平均値は最頻値とおおよそ一致する
    bgh_alpha = (bgh_mean / bgh_std) ** 2
    bgh_beta = bgh_mean / bgh_alpha
    bgl_alpha = (bgl_mean / bgl_std) ** 2
    bgl_beta = bgl_mean / bgl_alpha
    ## 端点近傍の平均値をガンマ分布の平均値，
    ## 端点近傍の標準偏差をガンマ分布の標準偏差の2倍とするように
    ## ガンマ分布の形状母数と尺度母数を求める
    bgh_alpha = (bgh_mean / (bgh_std * 2.0)) ** 2
    bgh_beta = bgh_mean / bgh_alpha
    bgl_alpha = (bgl_mean / (bgl_std * 2.0)) ** 2
    bgl_beta = bgl_mean / bgl_alpha

    bg_stddev = intensity_max / 10.0

    noiseTypeText = "gauss"
    if noiseType == "-p":
        noiseTypeText = "poisson"

    configfile = "config.txt"
    print("Output to ", configfile)

    lines = [
        f"MaxOfLoop :  10\n",
        f"generation :  50\n",
        # f"MaxOfLoop :  100\n",
        # f"generation :  50\n",
        f"individual :  30\n",
        f"num_island :   16\n",
        f"beta_log_upper :  1.25\n",
        f"beta_log_lower :  -2\n",
        f"gw_alpha :  {gw_alpha}\n",
        f"gw_beta :   {gw_beta}\n",
        f"lw_alpha :  {lw_alpha}\n",
        f"lw_beta :   {lw_beta}\n",
        f"bgh_alpha :  {bgh_alpha}\n",
        f"bgh_beta :   {bgh_beta}\n",
        f"bgl_alpha :  {bgl_alpha}\n",
        f"bgl_beta :   {bgl_beta}\n",
        f"intensity_alpha :  {intensity_alpha}\n",
        f"intensity_beta :  {intensity_beta}\n",
        f"intensity_max :  {intensity_max}\n",
        f"#bg_stddev :  {bg_stddev}\n",
        f"#shift_stddev : 0.2\n",
        f"#prob_exchange :  0.3\n",
        f"#prob_cross :  0.2\n",
        f"#prob_mutation : 0.2\n",
        f"noiseType : {noiseTypeText}",
    ]
    # print(lines)

    with open(configfile, mode="w") as f:
        f.writelines(lines)


#
def get_values_from_lines(lines, key):
    values = [line.strip().replace(key, "") for line in lines if line.startswith(key)]
    if len(values) == 0:
        print(f"Warning. Not match {key}")
    values = [int(v) for v in values]
    return values


#
def get_idH_idL_from_file(targetdir):
    # targetdir = './'
    idH = []
    idL = []
    logfile = targetdir + "/log_as_shirley.txt"
    if os.path.exists(logfile):
        with open(logfile) as f:
            lines = f.readlines()
            idH.extend(get_values_from_lines(lines, "idHigherGlobal = "))
            idL.extend(get_values_from_lines(lines, "idLowerGlobal  = "))
    logfile = targetdir + "/log_as_linear.txt"
    if os.path.exists(logfile):
        with open(logfile) as f:
            lines = f.readlines()
            idH.extend(get_values_from_lines(lines, "idHigherGlobal = "))
            idL.extend(get_values_from_lines(lines, "idLowerGlobal  = "))

    idHigherGlobal = min(idH) if len(idH) > 0 else None
    idLowerGlobal = max(idL) if len(idL) > 0 else None

    if (idHigherGlobal is None) or (idLowerGlobal is None):
        print("Error. idHigherGlobal or idLowerGlobal is None")
        sys.exit(1)

    print(idHigherGlobal, idLowerGlobal)
    return idHigherGlobal, idLowerGlobal


def get_idH_idL_from_file_2(targetdir):
    idHigherGlobal = idLowerGlobal = None

    target_file_path_list = glob(f"{targetdir}/global_energy_range.csv")
    if len(target_file_path_list) == 0:
        raise FileNotFoundError(f"{targetdir}/global_energy_range.csv")

    for target_file_path in target_file_path_list:
        with open(target_file_path, "r") as file:
            reader = csv.reader(file)
            for row in reader:
                idHigherGlobal, idLowerGlobal = map(int, row)

    if (idHigherGlobal is None) or (idLowerGlobal is None):
        raise ValueError("idHigherGlobal or idLowerGlobal is None")

    print(f"idHigherGlobal = {idHigherGlobal}")
    print(f"idLowerGlobal = {idLowerGlobal}")
    return idHigherGlobal, idLowerGlobal


def trim_spectrum_data_file(ifn, ofn, idH, idL):
    # print(ifn)
    df = pd.read_csv(ifn, header=None)
    if idH is None:
        idH = 0
    if idL is None:
        idL = len(df) - 1
    # print(df.head())
    # print(len(df), idH, idL)
    # print(type(len(df)), type(max(idH, idL)))
    if len(df) < (max(idH, idL) + 1):
        print("Error. Mismatch between file length and (idH, idL)", len(df), idH, idL)
        sys.exit(1)

    df = df.iloc[idH : idL + 1, :]
    # ofn = ifn.replace('.csv', '_trim.csv')
    df.to_csv(ofn, index=False, header=None, float_format="%g")
    print(f"Output to {ofn}")
    return


# def link_to_voigtDB(targetdir):
#     ## ターゲットディレクトリに VoigtDB.csv へのシンボリックリンクを作成
#     voigtDBsrc = f"{bindir}/VoigtDB.csv"
#     voigtDBdst = f"{targetdir}/VoigtDB.csv"
#     if os.path.islink(voigtDBdst):
#         os.unlink(voigtDBdst)
#     os.symlink(voigtDBsrc, voigtDBdst)


def jobAutoSimplex(spectrumfile, noiseType="-g", bgType="", devType=""):
    print("# -------------------------------------------------------")
    print("#        Auto Simplex")
    print("# -------------------------------------------------------")
    print("", flush=True)
    workDir = "auto"
    os.makedirs(workDir, exist_ok=True)

    # link_to_voigtDB(workDir)

    cmd = [PYTHON, PY_AUTO_SIMPLEX, spectrumfile, noiseType, bgType, devType]
    targetdir = os.path.join(os.getcwd(), workDir)
    if os.path.exists(spectrumfile):
        shutil.copy2(spectrumfile, targetdir)
    status = None
    status = subprocess.run(
        cmd,
        # shell = True,
        cwd=targetdir,
        # stdout = subprocess.PIPE,
        # stderr = subprocess.PIPE,
        check=False,
    )

    # status = exec_command(cmd, targetdir=targetdir)

    print(status)

    # idHigherGlobal, idLowerGlobal = get_idH_idL_from_file(targetdir)
    idHigherGlobal, idLowerGlobal = get_idH_idL_from_file_2(targetdir)

    trimmedSpectrumFile = spectrumfile.replace(".csv", "_trim.csv")
    trim_spectrum_data_file(spectrumfile, trimmedSpectrumFile, idHigherGlobal, idLowerGlobal)

    targetKey = "bic_gauss"
    if noiseType == "-p":
        targetKey = "bic_poisson"

    fnmBIC = targetdir + "/summary_BIC.csv"
    dfbic = pd.read_csv(fnmBIC)
    dfbic2 = dfbic.loc[dfbic.groupby("numPeak")[targetKey].idxmin()]
    numPeakOfMinBIC = dfbic["numPeak"][dfbic[targetKey].idxmin()]
    numPeakMax = numPeakOfMinBIC + 3

    print(f"numPeakOfMinBIC = {numPeakOfMinBIC}")
    print(f"numPeakMax = {numPeakMax}")

    return status, trimmedSpectrumFile, numPeakMax


def jobParamSimplex(spectrumfile, optparamfile, noiseType, bgType="", devType=""):
    print("# -------------------------------------------------------")
    print("#        Param Simplex")
    print("# -------------------------------------------------------")
    print("", flush=True)
    workDir = "param"
    os.makedirs(workDir, exist_ok=True)

    # link_to_voigtDB(workDir)

    cmd = [PYTHON, PY_PARAM_SIMPLEX, spectrumfile, optparamfile, noiseType, bgType, devType]
    targetdir = os.path.join(os.getcwd(), workDir)
    print(targetdir)
    if os.path.exists(optparamfile):
        shutil.copy2(optparamfile, targetdir)
    status = subprocess.run(
        cmd,
        cwd=targetdir,
        # stdout = subprocess.PIPE,
        # stderr = subprocess.PIPE,
        check=False,
    )
    return status


#
def collect_optparamfile(optparamfile):
    globkey = "./GA/num_peaks_*/opt_parameters_all_island.csv"
    print(f"Output to {optparamfile}")
    read_files = glob(globkey)
    if len(read_files) >= 1:
        with open(optparamfile, "wb") as outfile:
            for f in read_files:
                with open(f, "rb") as infile:
                    outfile.write(infile.read())


#
def jobGA(k, spectrumfile):
    # ndir = 'num_peaks_{:02}'.format(k)
    ndir = f"GA/num_peaks_{k:02d}"
    os.makedirs(ndir, exist_ok=True)

    # link_to_voigtDB(ndir)

    cmd = [
        EXE_GA,
        "123",
        str(k),
        spectrumfile,
    ]
    targetdir = os.path.join(os.getcwd(), ndir)
    print(targetdir)
    configfile = "config.txt"
    if os.path.exists(configfile):
        shutil.copy2(configfile, targetdir)
    if os.path.exists(spectrumfile):
        shutil.copy2(spectrumfile, targetdir)
    status = subprocess.run(
        cmd,
        cwd=targetdir,
        # stdout = subprocess.PIPE,
        # stderr = subprocess.PIPE,
        check=False,
    )
    return status


#
def plot_BIC_vs_numPeak(df3, targetKey="bic_gauss"):
    plt.rcParams["font.size"] = 18
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax.tick_params(axis="both", which="both")
    ax.grid(ls=":")

    if "method" in df3.columns:
        xoffset = 0.0
        df3S = df3[df3["method"] == "auto_shirley"]
        if len(df3S) > 0:
            ax.plot(df3S["numPeak"] + xoffset, df3S[targetKey], "o", mfc="None", c="C0", label="auto_simplex_shirley")
            xoffset += 0.2
        df3L = df3[df3["method"] == "auto_linear"]
        if len(df3L) > 0:
            ax.plot(df3L["numPeak"] + xoffset, df3L[targetKey], "s", mfc="None", c="C1", label="auto_simplex_linear")
            xoffset += 0.2
        df3GAS = df3[df3["method"] == "param_shirley"]
        if len(df3GAS) > 0:
            ax.plot(df3GAS["numPeak"] + xoffset, df3GAS[targetKey], "^", mfc="None", c="C2", label="param_simplex_shirley")
            xoffset += 0.2
        df3GAL = df3[df3["method"] == "param_linear"]
        if len(df3GAL) > 0:
            ax.plot(df3GAL["numPeak"] + xoffset, df3GAL[targetKey], "v", mfc="None", c="C3", label="param_simplex_linear")
            xoffset += 0.2
        ax.legend(loc="upper right", fontsize=18)
    else:
        ax.plot(df3["numPeak"], df3[targetKey], "o", mfc="None", c="C0", label="BIC")

    ax.set_xlabel("Number of peaks")
    ax.set_ylabel("BIC")
    # ax.grid(ls=":")
    # ax.legend(fontsize=14)
    # ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    # ax.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
    ax.minorticks_on()
    ax.tick_params(axis="x", which="minor", bottom=False, top=False)

    ## y軸の範囲を自動調整．各ピーク本数におけるBIC最小値に注目．
    df3tmp = df3.loc[df3.groupby("numPeak")[targetKey].idxmin()]
    ##
    print("Table of minimum BIC values vs number of peaks")
    print(df3tmp)
    ymin = df3tmp[targetKey].min()
    ymax = df3tmp[targetKey].max()
    yrange = ymax - ymin
    ax.set_ylim(ymin - 0.1 * yrange, ymax + 0.1 * yrange)

    fig.tight_layout()
    pngfile = f"BIC_vs_NumPeak.png"
    print("Output to ", pngfile)
    fig.savefig(pngfile)


# plot_result_core
def plot_result_core(csvfile, idHigherS, idLowerS, paper=False, ylim=None, figLevel=2, title=None):
    if figLevel <= 1:
        return
    if not os.path.exists(csvfile):
        print("Warning:", csvfile, "is not exists")
        return
    df = pd.read_csv(csvfile, comment="#")
    df["residue"] = df["spectrum"] - df["fitting"]

    if idHigherS < 0:
        idHigherS = 0
    if idLowerS < 0:
        idLowerS = len(df) - 1

    peakList = [s for s in df.columns if "peak[" in s]
    numPeakList = len(peakList)
    if numPeakList == 0:
        print("Warning: %s has no peak", csvfile)
        return

    peakLabel = "%d peaks" % numPeakList if (numPeakList > 1) else "1 peak"

    pngfile = csvfile.replace(".csv", ".png")
    print("Read data from ", csvfile, "  Output to", pngfile)

    plt.rcParams["font.size"] = 18
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    ax.plot(df["Binding_Energy(eV)"], df["spectrum"], "o", mfc="None", ms=5, label="spectrum")

    dfc = df[idHigherS : idLowerS + 1]
    ax.plot(dfc["Binding_Energy(eV)"], dfc["fitting"], lw=3, label="fitted", alpha=0.8)
    ax.plot(dfc["Binding_Energy(eV)"], dfc[peakList[0]], color="#333333", lw=1, label=peakLabel)
    for c in peakList[1:]:
        ax.plot(dfc["Binding_Energy(eV)"], dfc[c], color="#333333", lw=1, label="_nolegend_")
    ax.plot(dfc["Binding_Energy(eV)"], dfc["background"], c="C2", label="background")
    # ax.plot(dfc["Binding_Energy(eV)"], dfc["fitting"] - dfc["spectrum"], ".", c="C4", ms=4, label="residue")

    # Residue
    # ymax = dfc["spectrum"].max()
    # residueMax = dfc["residue"].max()
    # residueMin = dfc["residue"].min()
    datarange = dfc["spectrum"].max() - dfc["spectrum"].min()
    dataBaseline = dfc[["spectrum", "background"]].min().min()
    residueBaseline = dataBaseline - datarange * 0.15
    ylimMin = dataBaseline - datarange * 0.25
    ylimMin = min(ylimMin, dfc["background"].min())
    ax.axhline(residueBaseline, ls="--", c="gray")
    ax.plot(dfc["Binding_Energy(eV)"], dfc["residue"] + residueBaseline, "-", c="black", ms=2, lw=1, label="residue")

    if ylim is not None:
        ax.set_ylim(ylim)

    plot_xlabel = get_value_from_settings("plot_xlabel", "Binding energy (eV)")
    plot_xlabel = plot_xlabel.replace("_", " ")
    ax.set_xlabel(plot_xlabel)
    ax.set_ylabel("Intensity")

    ax.set_xlim(df["Binding_Energy(eV)"].max(), df["Binding_Energy(eV)"].min())
    ax.grid(ls=":")
    ax.legend(fontsize=14)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax.minorticks_on()
    ax.tick_params(axis="both", which="both")

    if title is not None:
        ax.set_title(title, fontsize=20, loc="right")

    plt.tight_layout()
    fig.savefig(pngfile, transparent=False)

    if paper:
        svgfile = pngfile.replace(".png", ".svg")
        print(f"Output to {svgfile}")
        fig.savefig(svgfile, transparent=False, dpi=100)

        # figdpi = 1000
        # figfile2 = pngfile.replace(".png", ".tif")
        # savefig_by_PIL(fig, figfile2, figdpi=figdpi)

    plt.clf()
    plt.close()


#
def copy_and_rename(fname1, fname2):
    """
    Copy a file from fname1 to fname2 after checking the file exists or not
    """
    if os.path.exists(fname1):
        print("Copy from %s to %s" % (fname1, fname2))
        shutil.copy2(fname1, fname2)
    else:
        print("Not found %s" % fname1)


#
def make_ranking(noiseType="-g"):
    targetKey = "bic_gauss"
    if noiseType == "-p":
        targetKey = "bic_poisson"

    fnm1 = "./auto/summary_BIC.csv"
    df1 = pd.DataFrame()
    if os.path.exists(fnm1):
        df1 = pd.read_csv(fnm1)
        df1["method"] = "auto_" + df1["method"]
        df1["ID"] = None
    fnm2 = "./param/summary_BIC.csv"
    df2 = pd.DataFrame()
    if os.path.exists(fnm2):
        df2 = pd.read_csv(fnm2)
        df2["method"] = "param_" + df2["method"]
        df2["smoothing_point"] = None
    dfv = pd.concat([df1, df2], axis=0, sort=True)
    if len(dfv) == 0:
        print("Error: No data found")
        sys.exit(1)
    df3p = dfv.sort_values(targetKey)
    df3p.reset_index(inplace=True, drop=True)

    # df3p = df3p.reindex(columns=['numPeak', targetKey, 'method', 'ID', 'smoothing_point'])
    # df3p = df3p.reindex(columns=["numPeak", "deviance_gauss", "bic_gauss", "deviance_poisson", "bic_poisson", "method", "ID", "smoothing_point"])
    df3p = df3p.reindex(columns=["numPeak", targetKey, "method", "smoothing_point", "ID"])
    df3p.rename(columns={"ID": "param_ID", "smoothing_point": "auto_ID"}, inplace=True)

    df3p.to_csv("./summary_BIC.csv", index=False)

    # BIC vs. number of peaks
    plot_BIC_vs_numPeak(df3p, targetKey=targetKey)

    # --------------------------------------------------------------

    # numPeakでくくり，それぞれ targetKey が最小となる行だけをピックアップする
    df3p = df3p.loc[df3p.groupby("numPeak")[targetKey].idxmin()]

    # データをBIC順にソートし，上位5件までのうちの最大のピーク数をnumPeakPickupMaxとして取得する．
    numPeakPickupMax = df3p.sort_values(targetKey)[:5]["numPeak"].max()

    # データをピーク本数順にソートする
    df3p.sort_values("numPeak", inplace=True)

    # ## numPeak <= numPeakMax の範囲に絞り込む
    # if numPeakMax is not None:
    #     df3p = df3p[df3p['numPeak'] <= numPeakMax]

    # numPeak <= numPeakPickupMax の範囲に絞り込む
    df3p = df3p[df3p["numPeak"] <= numPeakPickupMax]

    # インデックスを振り直す
    df3p.reset_index(drop=True, inplace=True)

    # print(df3p)

    # ランキング付けする最大順位
    # rankMax = min(5, len(df3p))

    # 全部取得する場合
    rankMax = len(df3p)

    # rankMaxまでで切り取る
    df3p = df3p.iloc[:rankMax, :]

    # targetKey の小さい順に並び替える。
    df3p = df3p.sort_values(targetKey, ascending=True)

    # print(df3p)
    # print("rankMax = ", rankMax)

    idHigherGlobal, idLowerGlobal = get_idH_idL_from_file_2("auto")

    # 第1位から第rankMax位までの図，パラメータファイル，スペクトルデータファイルを名前を変えてコピーする
    for i in np.arange(rankMax):
        rank = i + 1
        method = df3p["method"].values[i]
        numPeak = df3p["numPeak"].values[i]
        param_ID = df3p["param_ID"].values[i]
        auto_ID = df3p["auto_ID"].values[i]
        bicvalue = df3p[targetKey].values[i]
        conditionText = f"rank={rank}, method={method}, numPeak={numPeak}, param_ID={param_ID}, auto_ID={auto_ID}"
        # print(rank, method, numPeak, param_ID, auto_ID)
        # print(conditionText)

        fnms = []
        globkey = "file_not_found"
        if method.startswith("auto"):
            globkey = f"./auto/gbp_rank*_numPeak*_*_mSG{auto_ID:04d}_parameters.csv"
        elif method.startswith("param"):
            globkey = f"./param/gbp_rank*_numPeak*_*_ll{param_ID:04d}_parameters.csv"

        fnms = glob(globkey)
        if len(fnms) == 0:
            print(f"File not found. globkey = {globkey}")
            continue

        fnameBase1 = fnms[0]
        # fnameBase2 = f"gbp_rank{rank}_numPeak{numPeak}_{method}_{param_ID}_{auto_ID}_parameters.csv"
        fnameBase2 = f"gbp_rank{rank}_numPeak{numPeak}_{method}_parameters.csv"

        fname1 = fnameBase1
        fname2 = fnameBase2
        copy_and_rename(fname1, fname2)

        fname1 = fnameBase1.replace("_parameters.csv", "_result.csv")
        fname2 = fnameBase2.replace("_parameters.csv", "_result.csv")
        copy_and_rename(fname1, fname2)

        csvfile = fname2
        # title = f"{method}, BIC={bicvalue:.2f}"
        title = None
        plot_result_core(csvfile, idHigherGlobal, idLowerGlobal, title=title)

        # fname1 = fnameBase1.replace('_parameters.csv', '_result.png')
        # fname2 = fnameBase2.replace('_parameters.csv', '_result.png')
        # copy_and_rename(fname1, fname2)

        continue

    copy_and_rename("./auto/input_spectrum2.png", "input_spectrum2.png")

    return


#
if __name__ == "__main__":
    print("# -------------------------------------------------------")
    print("#        ", os.path.basename(__file__))
    print("# -------------------------------------------------------")
    print("", flush=True)

    noiseType = "-g"
    bgType = ""

    # if (len(sys.argv) == 1) or (len(sys.argv) > 4):
    if len(sys.argv) == 1:
        print(f"Error! Usage: {sys.argv[0]} <spectrum data file> [-s|-l] [-g|-p] [-x] [-d]")
        sys.exit(1)

    if "-s" in sys.argv:
        bgType = "-s"
    elif "-l" in sys.argv:
        bgType = "-l"
    print(f"bgType = {bgType}")

    if "-g" in sys.argv:
        noiseType = "-g"
    elif "-p" in sys.argv:
        noiseType = "-p"
    print(f"noiseType = {noiseType}")

    full_mode_flag = "--full" in sys.argv
    ppt_flag = "-x" in sys.argv
    devType = "-d" if "-d" in sys.argv else ""

    # link_to_voigtDB("./")

    spectrumfile = os.path.abspath(sys.argv[1])
    # print(f"spectrumfile = {spectrumfile}, K_max = {K_max}")
    print(f"spectrumfile = {spectrumfile}")

    trimmedSpectrumFile = spectrumfile.replace(".csv", "_trim.csv")
    numPeakMax = 6

    # 解析範囲の自動探索，初期モデル自動生成，バックグラウンド2種，Simplex法で最適化
    startTimeTmp = time.time()  # start time
    # print("debug20240709. skip jobAutoSimplex")
    status, trimmedSpectrumFile, numPeakMax = jobAutoSimplex(spectrumfile, noiseType, bgType, devType)
    td = time.time() - startTimeTmp
    print(f"AutoSimplex run time is {td:.2f} (sec) =  {datetime.timedelta(seconds=td)}", flush=True)

    if full_mode_flag:
        print("# -------------------------------------------------------")
        print("#        Genetic Algorithm")
        print("# -------------------------------------------------------")
        print("", flush=True)

        # GA計算を行うピーク本数の最大値
        K_max = numPeakMax
        print(f"K_max = {K_max}")

        # GA用のconfigファイルを生成
        make_config(spectrumfile, noiseType)

        # GAによる初期モデル探索
        startTimeTmp = time.time()  # start time
        for k in range(1, K_max + 1):
            print("< GA for k = {} >".format(k))
            status = jobGA(k, trimmedSpectrumFile)
            print(status)
        td = time.time() - startTimeTmp
        print(f"GA run time is {td:.2f} (sec) =  {datetime.timedelta(seconds=td)}", flush=True)

        # GAで得た初期モデルを用いてSimplex法で最適化を行う
        optparamfile = "opt_parameters_all_island.csv"
        collect_optparamfile(optparamfile)

        startTimeTmp = time.time()  # start time
        status = jobParamSimplex(trimmedSpectrumFile, optparamfile, noiseType, bgType, devType)
        print(status)
        # print(status.stdout.decode("utf8"))
        td = time.time() - startTimeTmp
        print(f"ParamSimplex run time is {td:.2f} (sec) =  {datetime.timedelta(seconds=td)}", flush=True)

    startTimeTmp = time.time()  # start time
    make_ranking(noiseType)
    td = time.time() - startTimeTmp
    print(f"make_ranking run time is {td:.2f} (sec) =  {datetime.timedelta(seconds=td)}", flush=True)

    # 無用の中間ファイルを削除
    if devType != "-d":
        for path in ["auto/", "GA/", "param/", "config.txt", "opt_parameters_all_island.csv", trimmedSpectrumFile]:
            if os.path.exists(path):
                if os.path.isfile(path):
                    os.remove(path)
                elif os.path.isdir(path):
                    shutil.rmtree(path)

    if ppt_flag:
        # 結果をpptにまとめる
        print("# -------------------------------------------------------")
        print("#        ", os.path.basename(PY_PPT))
        print("# -------------------------------------------------------")
        print("", flush=True)
        status = subprocess.run([PYTHON, PY_PPT], check=False)
        print(status)

    td = time.time() - startTime
    print(f"{sys.argv[0]} run time is {td:.2f} (sec) =  {datetime.timedelta(seconds=td)}", flush=True)
