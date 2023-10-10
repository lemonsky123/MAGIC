import os
import sh
import re
import csv
import time
import random
import warnings
import numpy as np
from tqdm import tqdm
import multiprocessing
import pandas as pd
from tqdm import tqdm
import scipy.stats as stats
from anndata import read_h5ad
from scipy import sparse, io
from itertools import groupby
from operator import itemgetter
from scipy.stats import pearsonr
import statsmodels.formula.api as smf
import statsmodels.stats.multitest as multi
from collections import Counter, defaultdict, OrderedDict
from colorama import init, Fore, Back
from magic_pipe import tools

init(autoreset=True)


def set_highlighted_excepthook():
    """
    this code is from https://newbedev.com/coloring-exceptions-from-python-on-a-terminal
    to color the exceptions.
    """
    import sys
    import traceback
    from pygments import highlight
    from pygments.lexers import get_lexer_by_name
    from pygments.formatters import TerminalFormatter

    lexer = get_lexer_by_name("pytb" if sys.version_info.major < 3 else "py3tb")
    formatter = TerminalFormatter()

    def myexcepthook(type, value, tb):
        tbtext = ''.join(traceback.format_exception(type, value, tb))
        sys.stderr.write(highlight(tbtext, lexer, formatter))

    sys.excepthook = myexcepthook


set_highlighted_excepthook()
path_tools = tools.__path__[0]


def red(s):
    return Fore.RED + s + Fore.RESET


def green(s):
    return Fore.GREEN + s + Fore.RESET


def yellow(s):
    return Fore.YELLOW + s + Fore.RESET


def blue(s):
    return Fore.BLUE + s + Fore.RESET


def magenta(s):
    return Fore.MAGENTA + s + Fore.RESET


def cyan(s):
    return Fore.CYAN + s + Fore.RESET


def white(s):
    return Fore.WHITE + s + Fore.RESET


def black(s):
    return Fore.BLACK


def white_green(s):
    return Fore.WHITE + Back.GREEN + s


class InvalidSampleTypeError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class GenesParsingError(BaseException):
    def __init__(self, genes):
        self.genes = genes

    def __str__(self):
        return "Unexpected Error occurred when parsing the following genes {}".format(",".join(self.genes))


class GenesPermNumberNotEqualError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class EmptyResultError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class SamplesNotFoundError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class SamplesNotEqualError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


def seurat_to_h5ad(path_rds, path_h5ad, path_rscript):
    if path_rscript == "Rscript" and not sh.which("Rscript"):
        raise FileNotFoundError("Could not find Rscript, please add Rscript in your path and try again.")
    path_seurat_to_h5ad = f"{path_tools}/scPBS_seurat_to_h5ad.R"
    rscript = sh.Command(path_rscript)
    rscript(path_seurat_to_h5ad, path_rds, path_h5ad)


def seurat_plot_cluster(path_rds, path_cor, path_pdf, path_rscript):
    if path_rscript == "Rscript" and not sh.which("Rscript"):
        raise FileNotFoundError("Could not find Rscript, please add Rscript in your path and try again.")
    path_seurat_plot_cluster = f"{path_tools}/scPBS_plot_cluster.R"
    rscript = sh.Command(path_rscript)
    rscript(path_seurat_plot_cluster, path_rds, path_cor, path_pdf)
