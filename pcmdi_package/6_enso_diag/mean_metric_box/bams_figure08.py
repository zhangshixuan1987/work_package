# -*- coding:UTF-8 -*- 
#---------------------------------------------------#
# Plot Figure 8 of Planton et al. (2021) https://doi.org/10.1175/BAMS-D-19-0337.1
#      - Read json file with all data (adapt to your json structure)
#           * 'reference' dataset is defined in the ENSO_metrics code but you can define your own
#      - Plot CMIP data as boxplot in parallel coordinates plot
#      - Plot members of selected model ('my_model') as markers and ensemble mean as line
#---------------------------------------------------#


#---------------------------------------------------#
# Import the right packages
#---------------------------------------------------#
from copy import deepcopy
from glob import iglob as GLOBiglob
import json
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from numpy.ma import array as NUMPYma__array
from numpy.ma import masked_invalid as NUMPYma__masked_invalid
from numpy.ma import masked_where as NUMPYmasked_where
from os.path import join as OSpath__join
from re import search as REsearch
from sys import exit as SYSexit
import time

# ENSO_metrics functions
from EnsoPlotLib import plot_param


#---------------------------------------------------#
# colors for printing
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
#---------------------------------------------------#


#---------------------------------------------------#
# Variables initialization
#---------------------------------------------------#
# arguments
metric_collection = ["ENSO_perf", "ENSO_proc", "ENSO_tel"]
experiment = "historical"
my_model = "IPSL-CM6A-LR"
my_project = "cmip6"
reduced_set = True  # False  #
# path
path_cmip = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2020_05_report/Data"
path_plot = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2019_10_ENSO_evaluation/Review/r01/r05"

list_project = ["cmip5", "cmip6"]

expe = "hist" if experiment == "historical" else "pi"

met_o1 = ["BiasPrLatRmse", "BiasPrLonRmse", "BiasSshLatRmse", "BiasSshLonRmse", "BiasSstLatRmse", "BiasSstLonRmse",
          "BiasTauxLatRmse", "BiasTauxLonRmse", "SeasonalPrLatRmse", "SeasonalPrLonRmse", "SeasonalSshLatRmse",
          "SeasonalSshLonRmse", "SeasonalSstLatRmse", "SeasonalSstLonRmse", "SeasonalTauxLatRmse",
          "SeasonalTauxLonRmse"]
met_o2 = ["EnsoSstLonRmse", "EnsoPrTsRmse", "EnsoSstTsRmse", "EnsoTauxTsRmse", "EnsoAmpl", "EnsoSeasonality",
          "EnsoSstSkew", "EnsoDuration", "EnsoSstDiversity_1", "EnsoSstDiversity_2", "NinoSstDiversity_1",
          "NinoSstDiversity_2"]
met_o3 = ["EnsoPrMapCorr", "EnsoPrMapRmse", "EnsoPrMapStd", "EnsoPrMapDjfCorr", "EnsoPrMapDjfRmse", "EnsoPrMapDjfStd",
          "EnsoPrMapJjaCorr", "EnsoPrMapJjaRmse", "EnsoPrMapJjaStd", "EnsoSlpMapCorr", "EnsoSlpMapRmse",
          "EnsoSlpMapStd", "EnsoSlpMapDjfCorr", "EnsoSlpMapDjfRmse", "EnsoSlpMapDjfStd", "EnsoSlpMapJjaCorr",
          "EnsoSlpMapJjaRmse", "EnsoSlpMapJjaStd", "EnsoSstMapCorr", "EnsoSstMapRmse", "EnsoSstMapStd",
          "EnsoSstMapDjfCorr", "EnsoSstMapDjfRmse", "EnsoSstMapDjfStd", "EnsoSstMapJjaCorr", "EnsoSstMapJjaRmse",
          "EnsoSstMapJjaStd"]
met_o4 = ["EnsoFbSstTaux", "EnsoFbTauxSsh", "EnsoFbSshSst",
          "EnsoFbSstThf", "EnsoFbSstSwr", "EnsoFbSstLhf", "EnsoFbSstLwr", "EnsoFbSstShf",
          "EnsodSstOce_1", "EnsodSstOce_2"]
met_order = met_o1 + met_o2 + met_o3 + met_o4

met_names = {
    "BiasPrLatRmse": "double_ITCZ_bias", "BiasPrLonRmse": "eq_PR_bias",
    "BiasSstLonRmse": "eq_SST_bias", "BiasTauxLonRmse": "eq_Taux_bias",
    "SeasonalPrLatRmse": "double_ITCZ_sea_cycle", "SeasonalPrLonRmse": "eq_PR_sea_cycle",
    "SeasonalSstLonRmse": "eq_SST_sea_cycle", "SeasonalTauxLonRmse": "eq_Taux_sea_cycle",
    "EnsoSstLonRmse": "ENSO_pattern", "EnsoSstTsRmse": "ENSO_lifecycle",
    "EnsoAmpl": "ENSO_amplitude", "EnsoSeasonality": "ENSO_seasonality",
    "EnsoSstSkew": "ENSO_asymmetry", "EnsoDuration": "ENSO_duration",
    "EnsoSstDiversity": "ENSO_diversity", "EnsoSstDiversity_1": "ENSO_diversity", "EnsoSstDiversity_2": "ENSO_diversity",
    "EnsoPrMapDjfRmse": "DJF_PR_teleconnection", "EnsoPrMapJjaRmse": "JJA_PR_teleconnection",
    "EnsoSstMapDjfRmse": "DJF_TS_teleconnection", "EnsoSstMapJjaRmse": "JJA_TS_teleconnection",
    "EnsoFbSstTaux": "SST-Taux_feedback", "EnsoFbTauxSsh": "Taux-SSH_feedback", "EnsoFbSshSst": "SSH-SST_feedback",
    "EnsoFbSstThf": "SST-NHF_feedback", 
    "EnsodSstOce": "ocean_driven_SST", "EnsodSstOce_1": "ocean_driven_SST", "EnsodSstOce_2": "ocean_driven_SST"}


#---------------------------------------------------#
print bcolors.OKGREEN + '%%%%%     -----     %%%%%'
print str().ljust(5)+"Parallel coordinate plot will be plotted"
print str().ljust(10)+'metric collection: '+str(metric_collection)
print str().ljust(10)+'experiment:        '+str(experiment)
print str().ljust(10)+'models:            '+str(my_project).upper()+' ' +str(my_model)
print '%%%%%     -----     %%%%%' + bcolors.ENDC
for ii in range(3): print ''


#---------------------------------------------------#
# Functions
#---------------------------------------------------#
def get_reference(metric_collection, metric):
    """
    Gets main reference for the given metric_collection / metric from EnsoPlotLib.plot_param

    Inputs:
    ------
    :param metric_collection: string
        Name of a metric collection.
    :param metric: string
        Name of a metric.

    Output:
    ------
    :return reference: string
        Name of the main reference for the given metric_collection / metric
    """
    if metric_collection in ["ENSO_tel", "test_tel"] and "Map" in metric:
        my_met = metric.replace("Corr", "").replace("Rmse", "").replace("Std", "")
    else:
        my_met = deepcopy(metric)
    reference = plot_param(metric_collection, my_met)["metric_reference"]
    return reference


def read_data1(project, metric_collection):
    lpath = OSpath__join(path_cmip, project + "/" + experiment)
    lname = project + "_" + experiment + "_" + metric_collection + "_v20200430.json"
    filename_js = list(GLOBiglob(OSpath__join(lpath, lname)))[0]
    with open(filename_js) as ff:
        data = json.load(ff)
    ff.close()
    return data["RESULTS"]["model"]


def read_data2(project, metric_collection, model):
    lpath = OSpath__join(path_cmip, project + "/" + experiment)
    lname = project + "_" + experiment + "_" + metric_collection + "_v20200427_allModels_allRuns.json"
    filename_js = list(GLOBiglob(OSpath__join(lpath, lname)))[0]
    with open(filename_js) as ff:
        data = json.load(ff)
    ff.close()
    return data["RESULTS"]["model"][model]


def remove_metrics(list_met, metric_collection):
    """
    Removes some metrics from given list

    Inputs:
    ------
    :param list_met: list of string
        List of metrics.
    :param metric_collection: string
        Name of a metric collection.

    Output:
    ------
    :return list_met_out: list of string
        Input list of metrics minus some metrics depending on given metric collection.
    """
    if metric_collection == "ENSO_perf":
        to_remove = ["BiasSshLatRmse", "BiasSshLonRmse", "BiasSstLatRmse", "BiasTauxLatRmse", "EnsoPrTsRmse",
                     "EnsoSstDiversity_1", "EnsoTauxTsRmse", "NinaSstDur_1", "NinaSstDur_2", "NinaSstLonRmse_1",
                     "NinaSstLonRmse_2", "NinaSstTsRmse_1", "NinaSstTsRmse_2", "NinoSstDiversity_1",
                     "NinoSstDiversity_2", "NinoSstDur_1", "NinoSstDur_2", "NinoSstLonRmse_1", "NinoSstLonRmse_2",
                     "NinoSstTsRmse_1", "NinoSstTsRmse_2", "SeasonalSshLatRmse", "SeasonalSshLonRmse",
                     "SeasonalSstLatRmse", "SeasonalTauxLatRmse"]
    elif metric_collection == "ENSO_proc":
        to_remove = ["BiasSshLonRmse", "EnsodSstOce_1", "EnsoFbSstLhf", "EnsoFbSstLwr", "EnsoFbSstShf", "EnsoFbSstSwr"]
    else:
        to_remove = ["EnsoPrMapCorr", "EnsoPrMapRmse", "EnsoPrMapStd", "EnsoPrMapDjfCorr", "EnsoPrMapDjfStd",
                     "EnsoPrMapJjaCorr", "EnsoPrMapJjaStd", "EnsoSlpMapCorr", "EnsoSlpMapRmse", "EnsoSlpMapStd",
                     "EnsoSlpMapDjfCorr", "EnsoSlpMapDjfRmse", "EnsoSlpMapDjfStd", "EnsoSlpMapJjaCorr",
                     "EnsoSlpMapJjaRmse", "EnsoSlpMapJjaStd", "EnsoSstMapCorr", "EnsoSstMapRmse", "EnsoSstMapStd",
                     "EnsoSstMapDjfCorr", "EnsoSstMapDjfStd", "EnsoSstMapJjaCorr", "EnsoSstMapJjaStd"]
    # remove given metrics
    list_met_out = sorted(list(set(list_met) - set(to_remove)), key=lambda v: v.upper())
    return list_met_out


dict_error = {
    "BiasPrLatRmse": [0, 4], "BiasPrLonRmse": [0, 3], "BiasSstLonRmse": [0, 3], "BiasTauxLonRmse": [0, 23],
    "SeasonalPrLatRmse": [0, 2.4], "SeasonalPrLonRmse": [0, 2.1], "SeasonalSstLonRmse": [0, 0.8],
    "SeasonalTauxLonRmse": [0, 8.8], "EnsoSstLonRmse": [0, 0.6], "EnsoSstTsRmse": [0, 0.4], "EnsoAmpl": [0, 80],
    "EnsoSeasonality": [0, 60], "EnsoSstSkew": [0, 240], "EnsoDuration": [0, 140], "EnsoSstDiversity": [0, 110],
    "EnsoPrMapDjfRmse": [0, 0.4], "EnsoPrMapJjaRmse": [0, 0.4], "EnsoSstMapDjfRmse": [0, 0.4],
    "EnsoSstMapJjaRmse": [0, 0.4], "EnsodSstOce": [0, 70], "EnsoFbSstThf": [0, 120], "EnsoFbSstTaux": [0, 80],
    "EnsoFbTauxSsh": [0, 65], "EnsoFbSshSst": [0, 110]}


def parallelplot(xlabels, ind0=0, fontsize=15, labels=None, title='', yname='', plotobs=False, dot=False, plot_legend=True, legend=[], cname=False, chigh=False, cfram=False):
    ax = plt.subplot(gs[:, ind0])
    xx = [ii for ii, _ in enumerate(xlabels)]
    data1 = dict((met, dict_met_mod1[met]["CMIP"]) for met in xlabels)
    data2 = dict((met, dict_met_mod1[met][my_model]) for met in xlabels)
    if plotobs is True:
        data4 = dict((met, dict_met_mod1[met]["obs"]) for met in xlabels)
    boxproperties = {
        "boxprops": dict(linestyle="-", linewidth=2, color=ref_colors["CMIP"]),
        "capprops": dict(linestyle="-", linewidth=2, color=ref_colors["CMIP"]),
        "flierprops": dict(marker="o", markersize=6.0, markeredgecolor=ref_colors["CMIP"], markerfacecolor=ref_colors["CMIP"], markeredgewidth=0),
        "meanprops":  dict(marker="D", markersize=15.0, markeredgecolor=ref_colors["CMIP"], markerfacecolor=ref_colors["CMIP"], markeredgewidth=0),
        "medianprops": dict(linestyle="-", linewidth=0, color=ref_colors["CMIP"]),
        "whiskerprops": dict(linestyle="-", linewidth=2, color=ref_colors["CMIP"]),
    }
    ax.set_title(title, fontsize=25, y=1.05, loc='left')
    ax.set_ylabel(yname, fontsize=25, labelpad=20)
    for dim, met in enumerate(xlabels):
        ax = plt.subplot(gs[:, ind0+dim])
        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.set_xlim([-1, 1])
        if labels is None:
            tmp_lab = [""]#[met] #[met + "\n(" + dict_met_mod2[met] + ")"]#
        else:
            tmp_lab = [""]#[labels[dim]] #[labels[dim] + "\n(" + dict_met_mod2[met] + ")"]#
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(15)
        for tick in ax.get_xticklabels():
            tick.set_rotation(90)
        if met == 'BiasPrLatRmse':
            mini, maxi = 0, 5
        elif met in ['BiasPrLonRmse', 'BiasSstLonRmse']:
            mini, maxi = 0, 3
        elif met == 'EnsoSstLonRmse':
            mini, maxi = 0.0, 0.6
        elif met == 'EnsoAmpl':
            mini, maxi = 0.4, 1.6
        elif met == 'EnsodSstOce':
            mini, maxi = 0.5, 3.5
        elif met == 'EnsoSeasonality':
            mini, maxi = 0.5, 2.5
        elif met == 'EnsoSstSkew':
            mini, maxi = -0.8, 0.8
        elif met == 'NinoSstDiversity':
            mini, maxi = 0, 70
        elif met == 'EnsoFbSshSst':
            mini, maxi = 0.1, 0.3
        elif met == 'EnsoFbSstTaux':
            mini, maxi = 0, 15
        elif met == 'EnsoFbSstThf':
            mini, maxi = -20, 5
        elif met == 'EnsoFbTauxSsh':
            mini, maxi = 0.1, 0.4
        mini, maxi = dict_error[met]
        ax.set_ylim([mini, maxi])
        #ax.set_yticks([mini, mini + ((maxi-mini) / 2.), maxi])
        #ax.set_yticklabels([str(mini), "", str(maxi)])
        ax.set_yticks([mini, maxi])
        ax.set_yticklabels([str(mini), str(maxi)])
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(15)
        ax.boxplot(data1[met], positions=[0], whis=[5, 95], widths=0.7, labels=tmp_lab, showmeans=True, showfliers=True, zorder=8, **boxproperties)
        if xlabels[dim] in met_o1 or xlabels[dim] + "_1" in met_o1 or xlabels[dim] + "_2" in met_o1:
            cc = "yellowgreen"
        elif xlabels[dim] in met_o2 or xlabels[dim] + "_1" in met_o2 or xlabels[dim] + "_2" in met_o2:
            cc = "plum"
        elif xlabels[dim] in met_o3 or xlabels[dim] + "_1" in met_o3 or xlabels[dim] + "_2" in met_o3:
            cc = "gold"
        else:
            cc = "turquoise"
        if cname is True:
            ax.text(0.1, -0.07, met_names[met], fontsize=15, ha='right', va='top', rotation=45, color=cc, transform=ax.transAxes)
        elif chigh is True:
            boxdict = dict(lw=0, facecolor=cc, pad=3, alpha=1)
            ax.text(0.5, -0.07, met_names[met], fontsize=15, ha='right', va='top', rotation=45, color="k", bbox=boxdict, transform=ax.transAxes)
        else:
            ax.text(0.5, -0.07, met_names[met], fontsize=15, ha='right', va='top', rotation=45, color="k", transform=ax.transAxes)
        if cfram is True:
            x1, x2 = ax.get_xlim()#ax.get_position().x0, ax.get_position().x1
            dx = (x2 - x1) / 100.
            y1, y2 = ax.get_ylim()#ax.get_position().y0, ax.get_position().y1
            dy = (y2 - y1) / 100.
            lix = [[x1-5*dx, x2+5*dx], [x1-5*dx, x2+5*dx]]
            liy = [[y2+5*dy, y2+5*dy], [y1-5*dy, y1-5*dy]]
            lic = [cc] * len (lix)
            lis = ["-"] * len(lix)
            liw = [5] * len(lix)
            for lc, ls, lw, lx, ly in zip(lic, lis, liw, lix, liy):
                line = Line2D(lx, ly, c=lc, lw=lw, ls=ls, zorder=10)
                line.set_clip_on(False)
                ax.add_line(line)
        ax.plot([-0.2] * len(data2[met]), data2[met], ls='None', marker=ref_markers[my_model], mec=ref_colors[my_model], mew=1, mfc=ref_colors[my_model], markersize=12, zorder=11, clip_on=False)
        ax.axhline(y=float(NUMPYma__array(data2[met]).mean()), zorder=6, xmin=0.16, xmax=0.86, color=ref_colors[my_model], lw=4)
        if plotobs is True:
            xxx = 0.2
            if ref_markers["obs"] == ">":
                xxx = -xxx
            ax.plot([xxx], data4[met], ls='None', marker=ref_markers["obs"], mec=ref_colors["obs"], mew=1, mfc=ref_colors["obs"], markersize=16, zorder=10, clip_on=False)
    if plot_legend is True:
        ax = plt.subplot(gs[:, 0])
        for kk, leg in enumerate(legend):
            font = {'color': ref_colors[leg], 'weight': 'bold', 'size': 20}
            ax.text(ref_positi[leg], -0.85, leg, ha="left", va="center", fontdict=font, transform=ax.transAxes)
    return
#---------------------------------------------------#


#---------------------------------------------------#
# Main
#---------------------------------------------------#
print "#---------------------------------------------------#"
print "read json"
print "#---------------------------------------------------#"
# cmip
dict_cmip = dict()
for proj in list_project:
    for mc in metric_collection:
        # open and read json file
        data_json = read_data1(proj, mc)
        # read metrics
        list_models = sorted(data_json.keys(), key=lambda v: v.upper())
        for mod in list_models:
            data_mod = data_json[mod][data_json[mod].keys()[0]]["value"]
            list_metrics = sorted(data_mod.keys(), key=lambda v: v.upper())
            if reduced_set is True:
                list_metrics = remove_metrics(list_metrics, mc)
            for met in list_metrics:
                my_ref = get_reference(mc, met)
                data_met = data_mod[met]["metric"]
                tmp = 1e20 if data_met[my_ref]["value"] is None else data_met[my_ref]["value"]
                try:    dict_cmip[mod]
                except: dict_cmip[mod] = {met: tmp}
                else:   dict_cmip[mod][met] = tmp
                del data_met, my_ref, tmp
            del data_mod, list_metrics
        del data_json, list_models
# models and metrics
tmp_models = sorted([str(mod) for mod in dict_cmip.keys()], key=lambda v: v.upper())
my_metrics = list()
for mod in tmp_models:
    try: dict_cmip[mod].keys()
    except: pass
    else: my_metrics += dict_cmip[mod].keys()
my_metrics = sorted(list(set(my_metrics)), key=lambda v: v.upper())
my_metrics = [met for met in met_order if met in my_metrics]
# shape dict by metric and put models in an array
dict_cmip_by_met = dict()
for met in my_metrics:
    tmp = list()
    for mod in tmp_models:
        try: dict_cmip[mod][met]
        except: pass
        else: tmp += [dict_cmip[mod][met]]
    tmp = NUMPYma__masked_invalid(NUMPYma__array(tmp))
    dict_cmip_by_met[met] = list(NUMPYmasked_where(tmp == 1e20, tmp).compressed())
    del tmp
# given model
dict_model = dict()
for mc in metric_collection:
    # open and read given model
    data_json = read_data2(my_project, mc, my_model)
    # read metrics
    list_members = sorted(data_json.keys(), key=lambda v: v.upper())
    for mem in list_members:
        data_mem = data_json[mem]["value"]
        list_metrics = sorted(data_mem.keys(), key=lambda v: v.upper())
        if reduced_set is True:
            list_metrics = remove_metrics(list_metrics, mc)
        for met in list_metrics:
            my_ref = get_reference(mc, met)
            data_met = data_mem[met]["metric"]
            tmp = 1e20 if data_met[my_ref]["value"] is None else data_met[my_ref]["value"]
            try:    dict_model[mem]
            except: dict_model[mem] = {met: tmp}
            else:   dict_model[mem][met] = tmp
            del data_met, my_ref, tmp
        del data_mem, list_metrics
    del data_json, list_members
# shape dict by metric and put members in an array
dict_model_by_met = dict()
for met in my_metrics:
    tmp = list()
    for mem in dict_model.keys():
        try: dict_model[mem][met]
        except: pass
        else: tmp += [dict_model[mem][met]]
    tmp = NUMPYma__masked_invalid(NUMPYma__array(tmp))
    dict_model_by_met[met] = list(NUMPYmasked_where(tmp == 1e20, tmp).compressed())
    del tmp
del dict_cmip, dict_model, tmp_models


#---------------------------------------------------#
# Plot data
#---------------------------------------------------#
dict_met_mod1 = dict()
for met in my_metrics:
    met2 = met.replace("_1", "").replace("_2", "")
    dict_met_mod1[met2] = {"CMIP": dict_cmip_by_met[met]}
    dict_met_mod1[met2][my_model] = dict_model_by_met[met]
    dict_met_mod1[met2]["obs"] = 0
    tmp = dict_cmip_by_met[met] + dict_model_by_met[met] + [0]
    print met.ljust(20)+str(round(min(tmp)))+" ; "+str(max(tmp))
    del met2, tmp
print "#---------------------------------------------------#"
print "plot data"
print "#---------------------------------------------------#"
if ' ':
    ref_colors = {"CMIP": "forestgreen", "CNRM-CM6-1": "peru", "IPSL-CM6A-LR": "crimson", "obs": "k"}
    ref_positi = {"CMIP": 12., "CNRM-CM6-1": 14.5, "IPSL-CM6A-LR": 14.5, "obs": 10.}
    ref_markers = {"CMIP": "D", "CNRM-CM6-1": ">", "IPSL-CM6A-LR": ">", "obs": "<"}
    #
    # Parallel plot
    #
    # metrics
    nbrl = 4
    nbrc = len(my_metrics)
    fig = plt.figure(0, figsize=(nbrc*0.7, nbrl))
    gs = GridSpec(nbrl, nbrc)
    title = ""#"metric values (% of error)"
    yname = "metric values"#""
    list_met = [met.replace("_1", "").replace("_2", "") for met in my_metrics]
    parallelplot(list_met, ind0=0, fontsize=15, labels=list_met, title=title, yname=yname, plotobs=True, dot=True, plot_legend=True, legend=["obs", "CMIP", my_model], cname=False, chigh=True, cfram=True)
    if isinstance(metric_collection, str) is True:
        name_png = OSpath__join(path_plot, "Figure_08_parallelplot_" + str(metric_collection) + "_" + str(my_model))
    else:
        name_png = OSpath__join(path_plot, "Figure_08_parallelplot_" + str(len(metric_collection)) + "mc_" + str(my_model))
    plt.savefig(name_png, bbox_inches='tight')
    plt.savefig(name_png + ".eps", bbox_inches="tight", format="eps")
    plt.close()


