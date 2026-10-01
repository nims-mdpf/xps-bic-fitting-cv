# import csv
# import pandas as pd
import os
import shutil
import pathlib
import subprocess
import sys
from glob import glob
import time, datetime
from multiprocessing import Pool


# Specify the directory which has "activeshirley_gcc.exe" and "plot_result.py"
bindir = os.path.abspath(os.path.dirname(__file__))
activeshirley_exe_file = f"{bindir}/active_shirley_from_params.exe"
plot_py_file = f"{bindir}/plot_result_shirley_and_linear_after_GA.py"
settingsfile = f"{bindir}/../settings.inp"
# settingsShirleyBG = f"{bindir}/settings_Shirley.inp"
# settingsLinearBG = f"{bindir}/settings_Linear.inp"
PYTHON = "python3"

print("bindir                 = ", bindir)
print("activeshirley_exe_file = ", activeshirley_exe_file)
print("plot_py_file           = ", plot_py_file)

if not os.path.exists(activeshirley_exe_file):
    print(f"Error. File not found {activeshirley_exe_file}")
    sys.exit(1)
if not os.path.exists(plot_py_file):
    print(f"Error. File not found {plot_py_file}")
    sys.exit(1)
# if not os.path.exists(settingsShirleyBG):
#     print(f"Error. File not found {settingsShirleyBG}")
#     sys.exit(1)
# if not os.path.exists(settingsLinearBG):
#     print(f"Error. File not found {settingsLinearBG}")
#     sys.exit(1)


def exec_command(cmd, targetdir="./", logfile=""):
    """
    Execute a command cmd assuming a Linux system.
    If specified, stdout and stderr will be output to logfile.
    Otherwise, the log is output to the standard console.

    >>> exec_command('echo hello')
    0
    """
    if logfile == "":
        p = subprocess.Popen(cmd, cwd=targetdir, shell=True)
    else:
        p = subprocess.Popen(cmd.split(), cwd=targetdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        with open(logfile, mode="w") as f:
            for line in p.stdout:
                line = line.decode("utf-8").rstrip()
                f.write(line + "\n")
    p.wait()
    return p.returncode


def main_core(fnm, optparamfile, kwd, noiseType="-g", devFlag=False):
    # pwd = os.getcwd()
    print(f"bgType = {kwd}, noiseType = {noiseType}")

    # kwdはshirleyまたはlinear
    # kwdディレクトリを作る
    targetdir = kwd
    os.makedirs(targetdir, exist_ok=True)

    # スペクトルデータファイルをコピーする
    shutil.copy2(fnm, targetdir)
    # パラメータファイルをコピーする
    shutil.copy2(optparamfile, targetdir)

    # デフォルトの settings.inp を targetdir にコピーする．
    # もしデフォルトの  settings.inp が存在しなければ， targetdir に空のファイルを作る．
    target_settings_file = f"{targetdir}/settings.inp"
    if os.path.exists(settingsfile):
        shutil.copy2(settingsfile, target_settings_file)
    else:
        pathlib.Path(target_settings_file).touch()

    # shirley の場合は settings.inp をそのまま使う．
    # linear の場合は settings.inp に backgroundType と allowedNegativeRatio の設定を追記する．
    if kwd == "shirley":
        os.path.exists(target_settings_file)
    elif kwd == "linear":
        with open(target_settings_file, "a") as f:
            print("backgroundType               Linear", file=f)
            print("allowedNegativeRatio         0.05", file=f)
    else:
        raise ValueError("kwd must be either shirley or linear")

    if noiseType == "-g":
        with open(target_settings_file, "a") as f:
            print("noiseType    Gauss", file=f)
    elif noiseType == "-p":
        with open(target_settings_file, "a") as f:
            print("noiseType    Poisson", file=f)
    else:
        raise ValueError("noiseType must be -g or -p.")

    # 自動ピーク分離を実行し，標準出力をlog_as_{kwd}.txtに書き出す
    cmd = f"{activeshirley_exe_file} {fnm} {optparamfile}"
    print(cmd)
    # exec_command(cmd, targetdir, f"log_as_{kwd}.txt")
    exec_command(cmd, targetdir)


def main(fnm, optparamfile, bgType="auto", noiseType="-g", figType="-minfig", devFlag=False):

    if (bgType == "auto") or (bgType == "shirley"):
        startTimeTmp = time.time()
        main_core(fnm, optparamfile, "shirley", noiseType=noiseType, devFlag=devFlag)
        td = time.time() - startTimeTmp
        print(f"shirley run time is {td:.2f} (sec) =  {datetime.timedelta(seconds=td)}")
    if (bgType == "auto") or (bgType == "linear"):
        startTimeTmp = time.time()
        main_core(fnm, optparamfile, "linear", noiseType=noiseType, devFlag=devFlag)
        td = time.time() - startTimeTmp
        print(f"linear run time is {td:.2f} (sec) =  {datetime.timedelta(seconds=td)}")

    # 事後処理を実行し，標準出力をlog_plot.txtに書き出す
    cmd = f"{PYTHON} {plot_py_file} {fnm} {figType} {noiseType}"
    print(cmd)
    startTimeTmp = time.time()
    # exec_command(cmd, logfile="log_plot.txt")
    exec_command(cmd)
    td = time.time() - startTimeTmp
    print(f"plot run time is {td:.2f} (sec) =  {datetime.timedelta(seconds=td)}")

    # 無用の中間ファイルを削除
    if not devFlag:
        # cmd = "rm -f out_*.csv parameters_*.csv prefered_*.png py_*.png result_* summary_evalFunc_minvalue.csv summary_evalFunc_vs_mSG.csv"
        cmd = "rm -rf linear shirley prefered_*.png py_*.png summary_evalFunc_minvalue.csv summary_evalFunc_vs_mSG.csv"
        print(cmd)
        exec_command(cmd)
    return


def show_usage():
    print(f"Usage:\n        {args[0]} <csv file> <optparamfile> [-d] [-s|-l] [-allfig|-minfig|-nofig] [-g|-p]")
    print(f"Example:\n        {args[0]} spectrumdata.csv")
    print("        <csv file> is a XPS data file in CSV-format.")
    print("        <optparamfile> is a parameter file in CSV-format.")
    print('        "-d" for developer mode')
    print('        "-s" for only Shirley background')
    print('        "-l" for only Linear background')
    print('        "-g" for Gauss noise')
    print('        "-p" for Poisson noise')
    print('        "-allfig" will create all figures')
    print('        "-minfig" will create only minimum necessary figures (default)')
    print('        "-nofig" will not create any figures')
    sys.exit(1)


#
if __name__ == "__main__":
    startTime = time.time()  # start time

    args = sys.argv

    bgType = "auto"
    noiseType = "-g"
    figType = "-minfig"

    # 引数に '-d' があれば開発者モード．全ての図を作る．中間ファイルを削除せずに残す．
    devFlag = "-d" in args
    if devFlag:
        figType = "-allfig"

    # 引数に '-l' も '-s' もなければデフォルトで auto とする．直線法とShirley法の両方を実行し，両者の中からモデル選択をする
    # 引数に '-l' があればバックグラウンドを直線法で近似する
    # 引数に '-s' があればバックグラウンドをShirley法で近似する
    # 引数に '-nofig', '-minfig', '-allfig' のいずれもなければデフォルトで '-minfig' とみなす
    # 引数に '-nofig' があれば図を作らない．時間短縮用
    # 引数に '-minfig' があれば必要最小限の図を作る．時間短縮用
    # 引数に '-allfig' があれば全ての図を作る．
    # 引数に '-g' も '-p' もなければデフォルトで '-g' とみなす
    # 引数に '-g' があればガウスノイズで扱う
    # 引数に '-p' があればポアソンノイズで扱う

    if "-s" in args:
        bgType = "shirley"
    elif "-l" in args:
        bgType = "linear"
    #
    if "-nofig" in args:
        figType = "-nofig"
    elif "-minfig" in args:
        figType = "-minfig"
    elif "-allfig" in args:
        figType = "-allfig"

    if "-g" in args:
        noiseType = "-g"
    elif "-p" in args:
        noiseType = "-p"

    args = [a for a in args if (not a.startswith("-"))]
    print(f"bgType = {bgType}")
    print(f"noiseType = {noiseType}")
    print(f"figType = {figType}")
    print(f"devFlag = {devFlag}")
    print(args)

    if len(args) <= 2:
        show_usage()

    main(args[1], args[2], bgType, noiseType, figType, devFlag)

    # print("Total script run time is --- %s" %(time.time()-startTime))
    td = time.time() - startTime
    print(f"{os.path.basename(sys.argv[0])} run time is {td:.2f} (sec) =  {datetime.timedelta(seconds=td)}")
