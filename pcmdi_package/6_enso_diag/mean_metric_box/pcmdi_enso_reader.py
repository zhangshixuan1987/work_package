import os
import glob
from typing import Dict, List, Tuple, Optional

import numpy as np
import xarray as xr

import difflib
import json
import string
from copy import deepcopy

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from numpy import isfinite as NUMPYisfinite
from numpy import linspace as NUMPYlinspace
from numpy import mean as NUMPYmean
from numpy import std as NUMPYstd
from numpy.ma import array as NUMPYma__array
from numpy.ma import masked_invalid as NUMPYma__masked_invalid
from numpy.ma import masked_where as NUMPYmasked_where
from numpy.ma import zeros as NUMPYma__zeros

# ENSO_metrics functions
from EnsoPlotLib import plot_param

class ENSOMetricReader:
    """
    Reader-only: load ENSO metrics from CLIVAR ENSO metrics JSON outputs and
    return numeric arrays + axis labels + reference bookkeeping.
    """

    def __init__(
        self,
        metric_collections,
        list_project,
        dict_json_path,
        list_obs=None,
        met_order=None,
        mod_order=None,
        reduced_set=False,
        sort_y_names=False,
        show_proj_means=False,
        show_ref_row=False,
        show_alt_obs_rows=False,
        debug=False,
    ):
        self.metric_collections = metric_collections
        self.list_project = list_project
        self.dict_json_path = dict_json_path
        self.list_obs = list_obs or []

        self.reduced_set = reduced_set
        self.sort_y_names = sort_y_names
        self.show_proj_means = show_proj_means
        self.show_ref_row = show_ref_row
        self.show_alt_obs_rows = show_alt_obs_rows
        self.debug = debug

        if met_order is None:
            met_order, *_ = self.load_met_order()
        self.met_order = met_order

        if mod_order == "predefined":
            mod_order = self.predefined_mod_order()
        self.mod_order = mod_order
        
    def most_similar_string(self, target, string_list):
        return max(
            string_list, key=lambda s: difflib.SequenceMatcher(None, target, s).ratio()
        )

    def read_metrics(self, return_raw=True):
        tab_norm, tab_raw, x_names, y_names, ref_info = (
            self.json_dict_to_numpy_array_list(
                self.metric_collections,
                self.list_project,
                self.list_obs,
                self.dict_json_path,
                self.reduced_set,
                self.met_order,
                self.mod_order,
                sort_y_names=self.sort_y_names,
                show_proj_means=self.show_proj_means,
                show_ref_row=self.show_ref_row,
                show_alt_obs_rows=self.show_alt_obs_rows,
                debug=self.debug,
            )
        )

        out = {
            "tab_norm": tab_norm,
            "x_names": x_names,
            "y_names": y_names,
            "ref_info": ref_info,
        }
        if return_raw:
            out["tab_raw"] = tab_raw
        return out

    def enso_portrait_plot(self, figure_name="enso_portrait_plot.png"):
        """
        Generates a summary plot for ENSO metrics.

        Returns
        -------
        fig : matplotlib.figure.Figure
        ref_info_dict : dict
        """
        # names for plot
        metric_names_for_plot, met_names = self.load_met_names()

        # met order groups needed by multiportraitplot
        met_order, met_o1, met_o2, met_o3, met_o4 = self.load_met_order()

        # get data (normalized arrays are what portrait plot expects)
        tab_all, _, x_names, y_names, ref_info_dict = self.json_dict_to_numpy_array_list(
            self.metric_collections,
            self.list_project,
            self.list_obs,
            self.dict_json_path,
            self.reduced_set,
            met_order,
            self.mod_order,
            sort_y_names=self.sort_y_names,
            show_proj_means=self.show_proj_means,
            show_ref_row=self.show_ref_row,
            show_alt_obs_rows=self.show_alt_obs_rows,
            debug=self.debug,
        )

        numbering = [ii + ") " for ii in list(string.ascii_lowercase)]
        title = [
            numbering[ii] + metric_names_for_plot[mc]
            for ii, mc in enumerate(self.metric_collections)
        ]

        if "CMIP6" in self.list_project and "CMIP5" in self.list_project:
            text = "* = CMIP6\nmodel"
        else:
            text = None

        levels = list(range(-2, 3))
        fig = self.multiportraitplot(
            tab_all,
            figure_name,
            x_names,
            y_names,
            title=title,
            my_text=text,
            levels=levels,
            highlight=True,
            met_o1=met_o1,
            met_o2=met_o2,
            met_o3=met_o3,
            met_o4=met_o4,
            met_names=met_names,
        )
        del levels, numbering, text, title
        return fig, ref_info_dict

    # ---------------------------------------------------#
    # Support Functions
    # ---------------------------------------------------#
    def my_colorbar(self, mini=-1.0, maxi=1.0, nbins=20):
        """
        Modifies cmo.balance colobar (removes the darkest blue and red)
    
        Inputs:
        ------
        **Optional arguments:**
        :param mini: float
            Minimum value of the colorbar.
        :param maxi: float
            Maximum value of the colorbar.
        :param nbins: integer
            Number of interval in the colorbar.
    
        Outputs:
        -------
        :return newcmp1: object
            Colormap, baseclass for all scalar to RGBA mappings
        :return norm: object
            Normalize, a class which can normalize data into the [0.0, 1.0] interval.
        """
        levels = MaxNLocator(nbins=nbins).tick_values(mini, maxi)
        # cmap = plt.get_cmap("cmo.balance")
        cmap = plt.get_cmap("RdBu_r")
        newcmp1 = cmap(NUMPYlinspace(0.15, 0.85, 256))
        newcmp2 = cmap(NUMPYlinspace(0.0, 1.0, 256))
        newcmp1 = ListedColormap(newcmp1)
        newcmp1.set_over(newcmp2[-30])
        newcmp1.set_under(newcmp2[29])
        newcmp1.set_bad(color="k")  # missing values in black
        norm = BoundaryNorm(levels, ncolors=newcmp1.N)
        return newcmp1, norm
        
    def multiportraitplot(
        self,
        tab,
        name_plot,
        x_names,
        y_names,
        title=[],
        write_metrics=False,
        my_text="",
        levels=None,
        highlight=False,
        nbr_space=2,
        met_o1=None,
        met_o2=None,
        met_o3=None,
        met_o4=None,
        met_names=None,
    ):
        """
        Plot the portraitplot (as in BAMS paper)
        """
        if levels is None:
            levels = [-1.0, -0.5, 0.0, 0.5, 1.0]
    
        # ---- minimal robustness guards (only if needed) ----
        met_o1 = met_o1 or []
        met_o2 = met_o2 or []
        met_o3 = met_o3 or []
        met_o4 = met_o4 or []
        met_names = met_names or {}
        # ---------------------------------------------------
    
        fontdict = {"fontsize": 40, "fontweight": "bold"}
        nbrc = sum([len(tab[ii][0]) for ii in range(len(tab))]) + (len(tab) - 1) * nbr_space
        fig = plt.figure(0, figsize=(0.5 * nbrc, 0.5 * len(tab[0])))
        gs = GridSpec(1, nbrc)
    
        cmap, norm = self.my_colorbar(mini=min(levels), maxi=max(levels))
    
        count = 0
        for kk, tmp in enumerate(tab):
            ax = plt.subplot(gs[0, count : count + len(tmp[0])])
            cs = ax.pcolormesh(tmp, cmap=cmap, norm=norm)
    
            xx1, xx2 = ax.get_xlim()
            dx = 0.5 / (xx2 - xx1)
            yy1, yy2 = ax.get_ylim()
            dy = 0.5 / (yy2 - yy1)
            try:
                ax.set_title(title[kk], fontdict=fontdict, y=1 + dy, loc="center")
            except Exception as e:
                print(f"An error occurred: {e}")
    
            ticks = [ii + 0.5 for ii in range(len(x_names[kk]))]
            ax.set_xticks(ticks)
            ax.set_xticklabels([] * len(ticks))
    
            for ll, txt in enumerate(x_names[kk]):
                label = met_names.get(txt, txt)  # ---- minimal fix for KeyError ----
    
                if highlight is True:
                    if txt in met_o1 or txt + "_1" in met_o1 or txt + "_2" in met_o1:
                        cc = "yellowgreen"
                    elif txt in met_o2 or txt + "_1" in met_o2 or txt + "_2" in met_o2:
                        cc = "plum"
                    elif txt in met_o3 or txt + "_1" in met_o3 or txt + "_2" in met_o3:
                        cc = "gold"
                    else:
                        cc = "turquoise"
    
                    ax.text(
                        ll + 0.5,
                        -0.2,
                        label,
                        fontsize=15,
                        ha="right",
                        va="top",
                        rotation=45,
                        color="k",
                        bbox=dict(lw=0, facecolor=cc, pad=3, alpha=1),
                    )
                else:
                    ax.text(
                        ll + 0.5,
                        -0.2,
                        label,
                        fontsize=20,
                        ha="right",
                        va="top",
                        rotation=45,
                        color="k",
                    )
    
            if highlight is True:
                tmp1 = [met_o1, met_o2, met_o3, met_o4]
    
                nn = 0
                lix = [[0, 0]]
                for tt in tmp1:
                    tmp2 = [
                        txt
                        for ll, txt in enumerate(x_names[kk])
                        if txt in tt or txt + "_1" in tt or txt + "_2" in tt
                    ]
                    nn += len(tmp2)
                    if len(tmp2) > 0:
                        lix += [[nn, nn]]
                    del tmp2
    
                liy = [[0, len(tab[0])]] * len(lix)
                lic, lis = ["k"] * len(lix), ["-"] * len(lix)
                for lc, ls, lx, ly in zip(lic, lis, lix, liy):
                    line = Line2D(lx, ly, c=lc, lw=7, ls=ls, zorder=10)
                    line.set_clip_on(False)
                    ax.add_line(line)
    
                nn = 0
                lic, lix = list(), list()
                for uu, tt in enumerate(tmp1):
                    tmp2 = [
                        txt
                        for ll, txt in enumerate(x_names[kk])
                        if txt in tt or txt + "_1" in tt or txt + "_2" in tt
                    ]
                    if len(tmp2) > 0:
                        if uu == 0:
                            cc = "yellowgreen"
                        elif uu == 1:
                            cc = "plum"
                        elif uu == 2:
                            cc = "gold"
                        else:
                            cc = "turquoise"
                        lic += [cc, cc]
                        if nn > 0:
                            lix += [[nn + 0.2, nn + len(tmp2)], [nn + 0.2, nn + len(tmp2)]]
                        else:
                            lix += [[nn, nn + len(tmp2)], [nn, nn + len(tmp2)]]
                        nn += len(tmp2)
                        del cc
                    del tmp2
    
                liy = [[len(tab[0]), len(tab[0])], [0, 0]] * int(float(len(lix)) / 2)
                lis = ["-"] * len(lix)
                for mm, (lc, ls, lx, ly) in enumerate(zip(lic, lis, lix, liy)):
                    if mm < 2:
                        line = Line2D([lx[0] + 0.05, lx[1]], ly, c=lc, lw=10, ls=ls, zorder=10)
                    elif mm > len(lis) - 3:
                        line = Line2D([lx[0], lx[1] - 0.05], ly, c=lc, lw=10, ls=ls, zorder=10)
                    else:
                        line = Line2D(lx, ly, c=lc, lw=10, ls=ls, zorder=10)
                    line.set_clip_on(False)
                    ax.add_line(line)
    
            ticks = [ii + 0.5 for ii in range(len(y_names))]
            ax.set_yticks(ticks)
            if kk != 0:
                ax.set_yticklabels([""] * len(ticks))
            else:
                ax.text(
                    -5 * dx,
                    -1 * dy,
                    my_text,
                    fontsize=25,
                    ha="right",
                    va="top",
                    transform=ax.transAxes,
                )
                ax.tick_params(axis="y", labelsize=20)
                ax.set_yticklabels(y_names)
            ax.yaxis.set_label_coords(-20 * dx, 0.5)
    
            for ii in range(1, len(tmp)):
                ax.axhline(ii, color="k", linestyle="-", linewidth=1)
            for ii in range(1, len(tmp[0])):
                ax.axvline(ii, color="k", linestyle="-", linewidth=1)
    
            if write_metrics is True:
                for jj in range(len(tmp[0])):
                    for ii in range(len(tmp)):
                        if tmp.mask[ii, jj] is False:
                            plt.text(
                                jj + 0.5,
                                ii + 0.5,
                                str(round(tmp[ii, jj], 1)),
                                fontsize=10,
                                ha="center",
                                va="center",
                            )
    
            if kk == len(tab) - 1:
                x2 = ax.get_position().x1
                y1 = ax.get_position().y0
                y2 = ax.get_position().y1
    
            count += len(tmp[0]) + nbr_space
    
        cax = plt.axes([x2 + 0.03, y1, 0.02, y2 - y1])
        cbar = plt.colorbar(
            cs,
            cax=cax,
            orientation="vertical",
            ticks=levels,
            pad=0.05,
            extend="both",
            aspect=40,
        )
        cbar.ax.set_yticklabels(
            [fr"-2 $\sigma$", "-1", "MMV", "1", fr"2 $\sigma$"],
            fontdict=fontdict,
        )
    
        dict_arrow = dict(facecolor="k", width=8, headwidth=40, headlength=40, shrink=0.0)
        dict_txt = dict(fontsize=40, rotation="vertical", ha="center", weight="bold")
    
        cax.annotate(
            "",
            xy=(3.7, 0.06),
            xycoords="axes fraction",
            xytext=(3.7, 0.45),
            arrowprops=dict_arrow,
        )
        cax.text(5.2, -0.55, "closer to reference", va="top", **dict_txt)
    
        cax.annotate(
            "",
            xy=(3.7, 0.94),
            xycoords="axes fraction",
            xytext=(3.7, 0.55),
            arrowprops=dict_arrow,
        )
        cax.text(5.2, 0.55, "further from reference", va="bottom", **dict_txt)
    
        plt.savefig(name_plot, bbox_inches="tight")
        return fig
        
        
    def find_first_member(self, members, mod=None):
        """
        Finds the first member
    
        Inputs:
        ------
        :param members: list of string
            List of members.
    
        Output:
        ------
        :return mem: string
            First member of the given list.
        """
        if "r1i1p1" in members:
            mem = "r1i1p1"
        elif "r1i1p1f1" in members:
            mem = "r1i1p1f1"
        elif "r1i1p1f2" in members:
            mem = "r1i1p1f2"
        else:
            tmp = deepcopy(members)
            members = list()
            for mem in tmp:
                for ii in range(1, 10):
                    if "r" + str(ii) + "i" in mem:
                        members.append(
                            mem.replace("r" + str(ii) + "i", "r" + str(ii).zfill(2) + "i")
                        )
                    else:
                        members.append(mem)
            del tmp
            mem = sorted(list(set(members)), key=lambda v: v.upper())[0].replace("r0", "r")
        # special case
        if mod == "NorESM2-LM":
            mem = "r2i1p1f1"
        return mem
    
    def json_dict_to_numpy_array_list(
        self,
        metric_collections,
        list_project,
        list_obs,
        dict_json_path,
        reduced_set,
        met_order,
        mod_order,
        sort_y_names=False,
        show_proj_means=False,
        show_ref_row=False,
        show_alt_obs_rows=False,
        debug=False,
    ):
        if debug:
            print("metric_collections:", metric_collections)
            print("list_project:", list_project)
            print("list_obs:", list_obs)

        model_by_proj = dict()
        dict_members = dict()
        for proj in list_project:
            list_models = list()
            dict_members[proj] = dict()
            for mc in metric_collections:
                tmp = ENSOMetricReader.read_data(dict_json_path[proj][mc])
                list_models += list(tmp.keys())
                for mod in list(tmp.keys()):
                    if mod not in dict_members[proj]:
                        dict_members[proj][mod] = list(tmp[mod].keys())
                    else:
                        dict_members[proj][mod] += list(tmp[mod].keys())
                del tmp

            list_models = sorted(list(set(list_models)), key=lambda v: v.upper())
            list_to_remove = [
                "EC-EARTH",
                "FIO-ESM",
                "GFDL-CM2p1",
                "HadGEM2-AO",
                "CIESM",
                "E3SM-1-1-ECA",
                "FGOALS-g3",
                "MCM-UA-1-0",
                "AWI-CM-1-1-MR",
                "AWI-ESM-1-1-LR",
            ]
            for mod in list_to_remove:
                while mod in list_models:
                    list_models.remove(mod)

            for mod in list_models:
                list_members = sorted(
                    list(set(dict_members[proj][mod])), key=lambda v: v.upper()
                )
                mem = self.find_first_member(list_members, mod=mod)
                if proj not in model_by_proj:
                    model_by_proj[proj] = {mod: mem}
                else:
                    if mod not in model_by_proj[proj]:
                        model_by_proj[proj][mod] = mem
                    else:
                        print("this model should not be here")
                del list_members, mem

            del list_models, list_to_remove

        # read json file
        tab_all, tab_all_act, x_names = list(), list(), list()
        different_ref_keys = list()
        ref_info_dict = dict()

        for mc in metric_collections:
            if debug:
                print("mc:", mc)
            dict1 = dict()
            list_models_all = list()
            ref_info_dict[mc] = dict()

            for proj in list_project:
                ref_info_dict[mc][proj] = dict()
                data_json = ENSOMetricReader.read_data(dict_json_path[proj][mc])

                list_models = sorted(
                    list(model_by_proj[proj].keys()), key=lambda v: v.upper()
                )
                list_models_all.extend(list_models)

                for mod in list_models:
                    data_mod = data_json[mod][model_by_proj[proj][mod]]["value"]
                    list_metrics = sorted(list(data_mod.keys()), key=lambda v: v.upper())
                    if reduced_set is True:
                        list_metrics = ENSOMetricReader.remove_metrics(list_metrics, mc)

                    for met in list_metrics:
                        ref = ENSOMetricReader.get_reference(mc, met)
                        ref_key_list = list(data_mod[met]["metric"])
                        ref_key_act = self.most_similar_string(ref, ref_key_list)

                        if ref != ref_key_act:
                            if debug:
                                print(
                                    f"Note: For metrics collection '{mc}' metric '{met}', "
                                    f"reference key in the JSON for the project '{proj}', "
                                    f"'{ref_key_act}', is assumed to be same as the predefined reference, '{ref}'."
                                )
                            different_ref_keys.append([ref, ref_key_act])

                        ref_info_dict[mc][proj][met] = ref_key_act
                        val = (
                            data_mod.get(met, {})
                            .get("metric", {})
                            .get(ref_key_act, {})
                            .get("value", None)
                        )
                        if val is None:
                            val = 1e20

                        if mod not in dict1:
                            dict1[mod] = {met: val}
                        else:
                            dict1[mod][met] = val
                        del ref, ref_key_act, val

                    del data_mod, list_metrics
                del data_json, list_models

            if len(different_ref_keys) > 0:
                # NOTE: this message is a bit misleading (proj is last loop var),
                # but keeping behavior unchanged except avoiding crash.
                unique_different_ref_keys = list(
                    map(list, dict.fromkeys(map(tuple, different_ref_keys)))
                )
                for diff_keys in unique_different_ref_keys:
                    if debug:
                        print(
                            f"Predefined reference: {diff_keys[0]}, reference key in the JSON: {diff_keys[1]}"
                        )

            # models and metrics
            if sort_y_names:
                tmp_models = sorted(
                    [str(mod) for mod in list(dict1.keys())], key=lambda v: v.upper()
                )
            else:
                tmp_models = list_models_all

            my_metrics = list()
            for mod in tmp_models:
                if mod in dict1:
                    my_metrics += list(dict1[mod].keys())

            my_metrics = sorted(list(set(my_metrics)), key=lambda v: v.upper())
            if met_order is not None:
                my_metrics = [met for met in met_order if met in my_metrics]

            if mod_order is not None:
                my_models = [mod for mod in mod_order if mod in tmp_models]
            else:
                my_models = tmp_models

            my_models += sorted(
                list(set(tmp_models) - set(my_models)), key=lambda v: v.upper()
            )
            my_models = list(reversed(my_models))
            del tmp_models

            rows_to_add = list()
            dict_ref_met = dict()

            if show_alt_obs_rows and list_obs is not None:
                rows_to_add += list_obs
                if len(list_obs) > 0 and "obs2obs" in dict_json_path.keys():
                    dict_ref_met = ENSOMetricReader.read_obs(
                        dict_json_path["obs2obs"][mc], list_obs, my_metrics, mc
                    )

            if show_ref_row:
                rows_to_add += ["reference"]

            if show_proj_means:
                rows_to_add += list(reversed(list_project))

            plus = len(rows_to_add)

            tab = NUMPYma__zeros((len(my_models) + plus, len(my_metrics)))
            for ii, mod in enumerate(my_models):
                for jj, met in enumerate(my_metrics):
                    if mod not in dict1 or met not in dict1[mod]:
                        tab[ii + plus, jj] = 1e20
                    else:
                        tab[ii + plus, jj] = dict1[mod][met]

            tab = NUMPYma__masked_invalid(tab)
            tab = NUMPYmasked_where(tab == 1e20, tab)
            tab_act = deepcopy(tab)

            # add values + normalize
            for jj, met in enumerate(my_metrics):
                tmp = tab[plus:, jj].compressed()
                if tmp.size == 0:
                    # nothing to normalize; leave column masked
                    for ii, dd in enumerate(rows_to_add):
                        if dd in list_obs:
                            tab[ii, jj] = dict_ref_met.get(dd, {}).get(met, 1e20)
                            tab_act[ii, jj] = tab[ii, jj]
                        elif dd in list_project:
                            tab[ii, jj] = 1e20
                            tab_act[ii, jj] = 1e20
                        else:
                            tab[ii, jj] = 0
                            tab_act[ii, jj] = 0
                    continue

                mea = float(NUMPYmean(tmp))
                std = float(NUMPYstd(tmp))
                del tmp

                for ii, dd in enumerate(rows_to_add):
                    if dd in list_obs:
                        val = dict_ref_met[dd][met]
                    elif dd in list_project:
                        tmp2 = [
                            tab[i2 + plus, jj]
                            for i2, m2 in enumerate(my_models)
                            if m2 in list(model_by_proj[dd].keys())
                        ]
                        tmp2 = NUMPYma__masked_invalid(NUMPYma__array(tmp2))
                        tmp2 = NUMPYmasked_where(tmp2 == 1e20, tmp2).compressed()
                        val = float(NUMPYmean(tmp2)) if tmp2.size > 0 else 1e20
                        del tmp2
                    else:
                        val = 0
                    tab[ii, jj] = val
                    tab_act[ii, jj] = val
                    del val

                # normalize safely
                if std == 0 or not NUMPYisfinite(std):
                    tab[:, jj] = 0.0
                else:
                    tab[:, jj] = (tab[:, jj] - mea) / std

                del mea, std

            tab = NUMPYma__masked_invalid(tab)
            tab = NUMPYmasked_where(tab > 1e3, tab)
            tab_act = NUMPYma__masked_invalid(tab_act)
            tab_act = NUMPYmasked_where(tab_act > 1e3, tab_act)

            tab_all.append(tab)
            tab_all_act.append(tab_act)

            if reduced_set is True:
                x_names.append([met.replace("_1", "").replace("_2", "") for met in my_metrics])
            else:
                x_names.append(my_metrics)

            if "CMIP6" in list_project and "CMIP5" in list_project:
                my_models = [
                    "* " + mod if mod in list(model_by_proj["CMIP6"].keys()) else mod
                    for mod in my_models
                ]

            if mc == metric_collections[0]:
                y_names = rows_to_add + my_models
                y_names = [
                    "(" + dd + ")" if dd in (list_obs + ["reference"]) else dd
                    for dd in y_names
                ]

            del dict1, dict_ref_met, my_metrics, my_models, plus, tab, tab_act

        return tab_all, tab_all_act, x_names, y_names, ref_info_dict

    @staticmethod
    def load_met_names():
        metric_names_for_plot = {
            "ENSO_perf": "Performance",
            "ENSO_proc": "Processes",
            "ENSO_tel": "Telecon.",
        }
        met_names = {
            "BiasPrLatRmse": "double_ITCZ_bias",
            "BiasPrLonRmse": "eq_PR_bias",
            "BiasSstLonRmse": "eq_SST_bias",
            "BiasTauxLonRmse": "eq_Taux_bias",
            "SeasonalPrLatRmse": "double_ITCZ_sea_cycle",
            "SeasonalPrLonRmse": "eq_PR_sea_cycle",
            "SeasonalSstLonRmse": "eq_SST_sea_cycle",
            "SeasonalTauxLonRmse": "eq_Taux_sea_cycle",
            "EnsoSstLonRmse": "ENSO_pattern",
            "EnsoSstTsRmse": "ENSO_lifecycle",
            "EnsoAmpl": "ENSO_amplitude",
            "EnsoSeasonality": "ENSO_seasonality",
            "EnsoSstSkew": "ENSO_asymmetry",
            "EnsoDuration": "ENSO_duration",
            "EnsoSstDiversity": "ENSO_diversity",
            "EnsoSstDiversity_1": "ENSO_diversity",
            "EnsoSstDiversity_2": "ENSO_diversity",
            "EnsoPrMapDjfRmse": "DJF_PR_teleconnection",
            "EnsoPrMapJjaRmse": "JJA_PR_teleconnection",
            "EnsoSstMapDjfRmse": "DJF_TS_teleconnection",
            "EnsoSstMapJjaRmse": "JJA_TS_teleconnection",
            "EnsoFbSstTaux": "SST-Taux_feedback",
            "EnsoFbTauxSsh": "Taux-SSH_feedback",
            "EnsoFbSshSst": "SSH-SST_feedback",
            "EnsoFbSstThf": "SST-NHF_feedback",
            "EnsodSstOce": "ocean_driven_SST",
            "EnsodSstOce_1": "ocean_driven_SST",
            "EnsodSstOce_2": "ocean_driven_SST",
        }
        return metric_names_for_plot, met_names

    @staticmethod
    def load_met_order():
        met_o1 = [
            "BiasPrLatRmse",
            "BiasPrLonRmse",
            "BiasSshLatRmse",
            "BiasSshLonRmse",
            "BiasSstLatRmse",
            "BiasSstLonRmse",
            "BiasTauxLatRmse",
            "BiasTauxLonRmse",
            "SeasonalPrLatRmse",
            "SeasonalPrLonRmse",
            "SeasonalSshLatRmse",
            "SeasonalSshLonRmse",
            "SeasonalSstLatRmse",
            "SeasonalSstLonRmse",
            "SeasonalTauxLatRmse",
            "SeasonalTauxLonRmse",
        ]
        met_o2 = [
            "EnsoSstLonRmse",
            "EnsoPrTsRmse",
            "EnsoSstTsRmse",
            "EnsoTauxTsRmse",
            "EnsoAmpl",
            "EnsoSeasonality",
            "EnsoSstSkew",
            "EnsoDuration",
            "EnsoSstDiversity",
            "EnsoSstDiversity_1",
            "EnsoSstDiversity_2",
            "NinoSstDiversity",
            "NinoSstDiversity_1",
            "NinoSstDiversity_2",
        ]
        met_o3 = [
            "EnsoPrMapCorr",
            "EnsoPrMapRmse",
            "EnsoPrMapStd",
            "EnsoPrMapDjfCorr",
            "EnsoPrMapDjfRmse",
            "EnsoPrMapDjfStd",
            "EnsoPrMapJjaCorr",
            "EnsoPrMapJjaRmse",
            "EnsoPrMapJjaStd",
            "EnsoSlpMapCorr",
            "EnsoSlpMapRmse",
            "EnsoSlpMapStd",
            "EnsoSlpMapDjfCorr",
            "EnsoSlpMapDjfRmse",
            "EnsoSlpMapDjfStd",
            "EnsoSlpMapJjaCorr",
            "EnsoSlpMapJjaRmse",
            "EnsoSlpMapJjaStd",
            "EnsoSstMapCorr",
            "EnsoSstMapRmse",
            "EnsoSstMapStd",
            "EnsoSstMapDjfCorr",
            "EnsoSstMapDjfRmse",
            "EnsoSstMapDjfStd",
            "EnsoSstMapJjaCorr",
            "EnsoSstMapJjaRmse",
            "EnsoSstMapJjaStd",
        ]
        met_o4 = [
            "EnsoFbSstTaux",
            "EnsoFbTauxSsh",
            "EnsoFbSshSst",
            "EnsoFbSstThf",
            "EnsoFbSstSwr",
            "EnsoFbSstLhf",
            "EnsoFbSstLwr",
            "EnsoFbSstShf",
            "EnsodSstOce",
            "EnsodSstOce_1",
            "EnsodSstOce_2",
        ]
        met_order = met_o1 + met_o2 + met_o3 + met_o4
        return met_order, met_o1, met_o2, met_o3, met_o4

    @staticmethod
    def predefined_mod_order():
        # model order
        mod_order = [
            "ACCESS1-0",
            "ACCESS1-3",
            "ACCESS-CM2",
            "ACCESS-ESM1-5",
            "BCC-CSM1-1",
            "BCC-CSM1-1-M",
            "BCC-CSM2-MR",
            "BCC-ESM1",
            "BNU-ESM",
            "CAMS-CSM1-0",
            "CanCM4",
            "CanESM2",
            "CanESM5",
            "CanESM5-CanOE",
            "CCSM4",
            "CESM1-BGC",
            "CESM1-CAM5",
            "CESM2",
            "CESM2-FV2",
            "CESM1-FASTCHEM",
            "CESM1-WACCM",
            "CESM2-WACCM",
            "CESM2-WACCM-FV2",
            "CMCC-CESM",
            "CMCC-CM",
            "CMCC-CMS",
            "CNRM-CM5",
            "CNRM-CM5-2",
            "CNRM-CM6-1",
            "CNRM-CM6-1-HR",
            "CNRM-ESM2-1",
            "CSIRO-Mk3-6-0",
            "CSIRO-Mk3L-1-2",
            "E3SM-1-0",
            "E3SM-1-1",
            "EC-EARTH",
            "EC-Earth3",
            "EC-Earth3-Veg",
            "FGOALS-f3-L",
            "FGOALS-g2",
            "FGOALS-s2",
            "FIO-ESM",
            "GFDL-CM2p1",
            "GFDL-CM3",
            "GFDL-CM4",
            "GFDL-ESM2G",
            "GFDL-ESM2M",
            "GFDL-ESM4",
            "GISS-E2-1-G",
            "GISS-E2-1-G-CC",
            "GISS-E2-H",
            "GISS-E2-H-CC",
            "GISS-E2-1-H",
            "GISS-E2-R",
            "GISS-E2-R-CC",
            "HadCM3",
            "HadGEM2-AO",
            "HadGEM2-CC",
            "HadGEM2-ES",
            "HadGEM3-GC31-LL",
            "INMCM4",
            "INM-CM4-8",
            "INM-CM5-0",
            "IPSL-CM5A-LR",
            "IPSL-CM5A-MR",
            "IPSL-CM5B-LR",
            "IPSL-CM6A-LR",
            "KACE-1-0-G",
            "MIROC4h",
            "MIROC5",
            "MIROC6",
            "MIROC-ESM",
            "MIROC-ESM-CHEM",
            "MIROC-ES2L",
            "MPI-ESM-LR",
            "MPI-ESM-MR",
            "MPI-ESM-P",
            "MPI-ESM-1-2-HAM",
            "MPI-ESM1-2-HR",
            "MPI-ESM1-2-LR",
            "MRI-CGCM3",
            "MRI-ESM1",
            "MRI-ESM2-0",
            "NESM3",
            "NorESM1-M",
            "NorESM1-ME",
            "NorCPM1",
            "NorESM2-LM",
            "NorESM2-MM",
            "SAM0-UNICON",
            "TaiESM1",
            "UKESM1-0-LL",
        ]
    
        mod_order += [
            "CAS-ESM2-0",
            "CMCC-CM2-HR4",
            "CMCC-CM2-SR5",
            "EC-Earth3-AerChem",
            "EC-Earth3-Veg-LR",
            "FIO-ESM-2-0",
            "HadGEM3-GC31-MM",
            "KIOST-ESM",
        ]
    
        mod_order = sorted(mod_order, key=str.casefold)
    
        return mod_order
        
    @staticmethod
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
        
    @staticmethod
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
            to_remove = [
                "BiasSshLatRmse",
                "BiasSshLonRmse",
                "BiasSstLatRmse",
                "BiasTauxLatRmse",
                "EnsoPrTsRmse",
                "EnsoSstDiversity_1",
                "EnsoTauxTsRmse",
                "NinaSstDur_1",
                "NinaSstDur_2",
                "NinaSstLonRmse_1",
                "NinaSstLonRmse_2",
                "NinaSstTsRmse_1",
                "NinaSstTsRmse_2",
                "NinoSstDiversity_1",
                "NinoSstDiversity_2",
                "NinoSstDur_1",
                "NinoSstDur_2",
                "NinoSstLonRmse_1",
                "NinoSstLonRmse_2",
                "NinoSstTsRmse_1",
                "NinoSstTsRmse_2",
                "SeasonalSshLatRmse",
                "SeasonalSshLonRmse",
                "SeasonalSstLatRmse",
                "SeasonalTauxLatRmse",
            ]
        elif metric_collection == "ENSO_proc":
            to_remove = [
                "BiasSshLonRmse",
                "EnsodSstOce_1",
                "EnsoFbSstLhf",
                "EnsoFbSstLwr",
                "EnsoFbSstShf",
                "EnsoFbSstSwr",
            ]
        else:
            to_remove = [
                "EnsoPrMapCorr",
                "EnsoPrMapRmse",
                "EnsoPrMapStd",
                "EnsoPrMapDjfCorr",
                "EnsoPrMapDjfStd",
                "EnsoPrMapJjaCorr",
                "EnsoPrMapJjaStd",
                "EnsoSlpMapCorr",
                "EnsoSlpMapRmse",
                "EnsoSlpMapStd",
                "EnsoSlpMapDjfCorr",
                "EnsoSlpMapDjfRmse",
                "EnsoSlpMapDjfStd",
                "EnsoSlpMapJjaCorr",
                "EnsoSlpMapJjaRmse",
                "EnsoSlpMapJjaStd",
                "EnsoSstMapCorr",
                "EnsoSstMapRmse",
                "EnsoSstMapStd",
                "EnsoSstMapDjfCorr",
                "EnsoSstMapDjfStd",
                "EnsoSstMapJjaCorr",
                "EnsoSstMapJjaStd",
            ]
        # remove given metrics
        list_met_out = sorted(list(set(list_met) - set(to_remove)), key=lambda v: v.upper())
        return list_met_out

    @staticmethod
    def read_data(filename_json):
        """
        Reads given json file (must have usual PMP's structure)
    
        Input:
        -----
        :param filename_json: string
            Path and name of a json file output of the CLIVAR ENSO metrics package.
    
        Output:
        ------
        :return data: dictionary
            Dictionary output of the CLIVAR ENSO metrics package, first level is models, second is members.
        """
        with open(filename_json) as ff:
            data = json.load(ff)
        ff.close()
        data = data["RESULTS"]["model"]
        return data
        
    @staticmethod
    def read_obs(filename_json, obsvation_names, list_met, metric_collection):
        """
        Reads given json file (must have usual PMP's structure) and read given obs
    
        Input:
        -----
        :param filename_json: string
            Path and name of a json file output of the CLIVAR ENSO metrics package.
        :param obsvation_names: list of string
            Names of wanted additional observations for the portrait plot
        :param list_met: list of string
            List of metrics.
        :param metric_collection: string
            Name of a metric collection.
    
        Output:
        ------
        :return data: list
            Dictionary output of additional observations metric values.
        """
        data_json = read_data(filename_json)
        dict_out = dict()
        for obs in obsvation_names:
            for met in list_met:
                ref = get_reference(metric_collection, met)
                if obs == "20CRv2":
                    if "Ssh" not in met:
                        try:
                            tab = data_json["20CRv2"]["r1i1p1"]["value"][met]["metric"]
                        except KeyError:
                            tab = data_json["20CRv2_20CRv2"]["r1i1p1"]["value"][met][
                                "metric"
                            ]
                elif obs == "NCEP2":
                    if "TauxSsh" in met or "SshSst" in met:
                        tab = data_json["NCEP2_GODAS"]["r1i1p1"]["value"][met]["metric"]
                    elif "Ssh" in met:
                        tab = data_json["GODAS"]["r1i1p1"]["value"][met]["metric"]
                    else:
                        try:
                            tab = data_json["NCEP2"]["r1i1p1"]["value"][met]["metric"]
                        except KeyError:
                            tab = data_json["NCEP2_NCEP2"]["r1i1p1"]["value"][met]["metric"]
                elif obs == "ERA-Interim":
                    if "SstMap" in met:
                        tab = {ref: {"value": 0}}
                    elif "TauxSsh" in met or "SshSst" in met:
                        tab = data_json["ERA-Interim_SODA3.4.2"]["r1i1p1"]["value"][met][
                            "metric"
                        ]
                    elif "Ssh" in met:
                        tab = data_json["SODA3.4.2"]["r1i1p1"]["value"][met]["metric"]
                    else:
                        try:
                            tab = data_json["ERA-Interim"]["r1i1p1"]["value"][met]["metric"]
                        except KeyError:
                            tab = data_json["ERA-Interim_ERA-Interim"]["r1i1p1"]["value"][
                                met
                            ]["metric"]
    
                try:
                    val = tab[ref]["value"]
                except KeyError:
                    val = 1e20
    
                try:
                    dict_out[obs]
                except KeyError:
                    dict_out[obs] = {met: val}
                else:
                    dict_out[obs][met] = val
    
                del ref, val
        return dict_out

class ENSODiagReader:
    """
    Helper to read ENSO diagnostics (perf/proc/tel) across
    historical vs future and ensemble members.

    Usage:
        reader = ENSODiagReader(...)
        da_hist = reader.load("ENSO_perf", "enso_amplitude", period="hist")
        da_fut  = reader.load("ENSO_perf", "enso_amplitude", period="future")
    """

    def __init__(
        self,
        data_dir: str,
        model: str,
        groups: List[str],
        period_list: List[Tuple[int, int]],
        nens: List[int],
        enso_groups: Optional[Dict[str, List[str]]] = None,
        file_suffix_map: Optional[Dict[str, str]] = None,
        members: Optional[List[int]] = None,
        verbose: bool = False,
    ):
        self.data_dir = data_dir
        self.model = model
        self.groups = groups
        self.period_list = period_list
        self.nens = nens

        # Use provided maps if given; otherwise fall back to internal defaults
        self.enso_groups = enso_groups if enso_groups is not None else self.get_enso_var()
        self.file_suffix_map = (
            file_suffix_map if file_suffix_map is not None else self.get_file_suffix()
        )

        # Ensemble member IDs; if not given, infer from directory listing
        self.members = members  # e.g. [51, 91, 101, ...]
        self._cached_member_dirs: Dict[str, List[str]] = {}  # (group_key) -> [paths]

        # Verbosity
        self.verbose = verbose

    # ---------------------------------------------------------
    # Dictionary getters (FLEXIBLE)
    # ---------------------------------------------------------
    def get_enso_var(self) -> Dict[str, List[str]]:
        """Return mapping: ENSO group -> list of variable names."""
        return {
            "ENSO_perf": [
                "pr_lat_rmse", "pr_lon_rmse", "sst_lon_rmse", "taux_lon_rmse",
                "enso_amplitude", "enso_duration", "enso_seasonality",
                "enso_sst_diversity_mode1", "enso_sst_diversity_mode2",
                "enso_sst_lon_rmse", "enso_sst_skewness", "enso_sst_ts_rmse",
                "seasonal_pr_lat_rmse", "seasonal_pr_lon_rmse",
                "seasonal_sst_lon_rmse", "seasonal_taux_lon_rmse",
            ],

            "ENSO_proc": [
                "sst_lon_rmse", "taux_lon_rmse", "enso_amplitude",
                "enso_dsst_oce_mode1", "enso_dsst_oce_mode2",
                "enso_fb_ssh_sst", "enso_fb_sst_taux", "enso_fb_sst_thf",
                "enso_fb_taux_ssh", "enso_seasonality",
                "enso_sst_lon_rmse", "enso_sst_skewness",
            ],

            "ENSO_tel": [
                "enso_amplitude", "enso_pr_map_djf", "enso_pr_map_jja",
                "enso_seasonality", "enso_sst_lon_rmse",
                "enso_sst_map_djf", "enso_sst_map_jja",
            ],
        }

    def get_file_suffix(self) -> Dict[str, str]:
        """Return mapping: logical variable name → file suffix used in filenames."""
        return {
            # ENSO_perf
            "pr_lat_rmse":              "BiasPrLatRmse",
            "pr_lon_rmse":              "BiasPrLonRmse",
            "sst_lon_rmse":             "BiasSstLonRmse",
            "taux_lon_rmse":            "BiasTauxLonRmse",
            "enso_amplitude":           "EnsoAmpl",
            "enso_duration":            "EnsoDuration",
            "enso_seasonality":         "EnsoSeasonality",
            "enso_sst_diversity_mode1": "EnsoSstDiversity_1",
            "enso_sst_diversity_mode2": "EnsoSstDiversity_2",
            "enso_sst_lon_rmse":        "EnsoSstLonRmse",
            "enso_sst_skewness":        "EnsoSstSkew",
            "enso_sst_ts_rmse":         "EnsoSstTsRmse",
            "seasonal_pr_lat_rmse":     "SeasonalPrLatRmse",
            "seasonal_pr_lon_rmse":     "SeasonalPrLonRmse",
            "seasonal_sst_lon_rmse":    "SeasonalSstLonRmse",
            "seasonal_taux_lon_rmse":   "SeasonalTauxLonRmse",

            # ENSO_proc
            "enso_dsst_oce_mode1":      "EnsoDeltaSstOceMode1",
            "enso_dsst_oce_mode2":      "EnsoDeltaSstOceMode2",
            "enso_fb_ssh_sst":          "EnsoFbSshSst",
            "enso_fb_sst_taux":         "EnsoFbSstTaux",
            "enso_fb_sst_thf":          "EnsoFbSstThf",
            "enso_fb_taux_ssh":         "EnsoFbTauxSsh",

            # ENSO_tel
            "enso_pr_map_djf":          "EnsoPrMapDJF",
            "enso_pr_map_jja":          "EnsoPrMapJJA",
            "enso_sst_map_djf":         "EnsoSstMapDJF",
            "enso_sst_map_jja":         "EnsoSstMapJJA",
        }

    # ------------------------------ small helpers ------------------------------

    def available_groups(self) -> List[str]:
        """Return list of available ENSO groups."""
        return list(self.enso_groups.keys())

    def available_vars(self, enso_group: str) -> List[str]:
        """Return list of variable names for a given ENSO group."""
        if enso_group not in self.enso_groups:
            raise ValueError(
                f"Unknown ENSO group '{enso_group}'. "
                f"Must be one of {list(self.enso_groups.keys())}"
            )
        return self.enso_groups[enso_group]

    # ------------------------------ core helpers ------------------------------
    def _unify_longitude_name(
        self,
        ds: xr.Dataset,
        coord_candidates=("longitude", "lon", "LONGITUDE", "LON"),
        unified_name="longitude",
    ) -> xr.Dataset:
        lon_name = None
        for cand in coord_candidates:
            if cand in ds.coords:
                lon_name = cand
                break

        if lon_name is None:
            return ds

        if lon_name != unified_name:
            ds = ds.rename({lon_name: unified_name})

        return ds

    def _get_period_index(self, period: str) -> int:
        if period not in self.groups:
            raise ValueError(f"period must be one of {self.groups}, got {period}")
        return self.groups.index(period)

    def _get_enso_case_name(self, period: str) -> str:
        idx = self._get_period_index(period)
        start, end = self.period_list[idx]
        return f"ENSO_{start}-{end}"

    def _list_member_dirs(self, period: str) -> List[str]:
        """
        Return the list of member directories for a given period, e.g.
        [..., '<DATA_DIR>/hist/v3.LR.historical_0051', ...].
        """
        key = period
        if key in self._cached_member_dirs:
            return self._cached_member_dirs[key]

        base_dir = os.path.join(self.data_dir, period)

        # match v3.LR.historical_****
        pattern = os.path.join(base_dir, f"{self.model}_*")
        dirs = sorted(d for d in glob.glob(pattern) if os.path.isdir(d))

        # If user specified explicit members, filter accordingly
        if self.members is not None:
            keep = []
            for d in dirs:
                mstr = os.path.basename(d).split("_")[-1]
                try:
                    mid = int(mstr)
                except ValueError:
                    continue
                if mid in self.members:
                    keep.append(d)
            dirs = keep

        # Optionally trim to NENS[period_index]
        idx = self._get_period_index(period)
        if len(dirs) > self.nens[idx]:
            dirs = dirs[: self.nens[idx]]

        self._cached_member_dirs[key] = dirs
        return dirs

    def _find_nc_file(
        self,
        member_dir: str,
        enso_group: str,
        suffix: str,
        case_id: str, 
    ) -> str:
        """
        Find the NetCDF file in the enso_group directory whose name ends with
        '_<suffix>.nc'. We keep it flexible w.r.t. the date stamp, etc.
        """
        enso_root = os.path.join(
            member_dir,
            "pcmdi_diags",
            "model_vs_obs",
            "metrics_data",
            "enso_metric",
            enso_group,
        )

        if not os.path.isdir(enso_root):
            raise FileNotFoundError(f"ENSO directory not found: {enso_root}")

        pattern = os.path.join(enso_root, f"*{case_id}*_{suffix}.nc")
        matches = sorted(glob.glob(pattern))
        print(pattern)
        if not matches:
            raise FileNotFoundError(f"No file matching *_{suffix}.nc in {enso_root}")
        if len(matches) > 1:
            # In practice you likely have only one; if multiple, take the last.
            return matches[-1]

        if self.verbose:
            print(f"enso metrics file found: {matches[0]}")

        return matches[0]

    def _choose_default_var(
        self,
        ds: xr.Dataset,
        candidates: List[str],
        member_str: str,
    ) -> str:
        """
        Heuristic to choose a model diagnostic variable when nc_var is not provided.
        Prefer:
          1) vars containing member_str,
          2) 1D longitude vars over 2D maps,
          3) first candidate as a final fallback.
        """
        # prefer vars that contain this member id in the name
        cand_member = [v for v in candidates if member_str in v]

        # prefer 1D longitude vars (no latitude) for amplitude-like metrics
        cand_lon = [
            v for v in cand_member
            if ("longitude" in ds[v].dims and "latitude" not in ds[v].dims)
        ]
        if len(cand_lon) == 1:
            return cand_lon[0]
        if len(cand_lon) > 1:
            return cand_lon[0]

        if len(cand_member) == 1:
            return cand_member[0]
        if len(cand_member) > 1:
            return cand_member[0]

        # if all else fails, just take the first candidate
        return candidates[0]

    def _extract_obs(
        self,
        ds: xr.Dataset,
        base_var: str = "sstStd_lon",
        ref_tag: str = "ERA-Interim",
    ) -> xr.DataArray:
        """
        Load the observational ENSO amplitude zonal std curve
        (e.g., sstStd_lon__ERA-Interim(longitude)) from one file.

        Parameters
        ----------
        base_var : str
            Base pattern for the variable name, e.g., "sstStd_lon".
        ref_tag : str
            Substring that identifies the obs variable name,
            e.g., "ERA-Interim".

        Returns
        -------
        obs_da : DataArray(longitude)
        """
        if base_var is None or ref_tag is None:
            raise ValueError(
                "base_var and ref_tag must be non-None when extracting observations."
            )

        # require BOTH base pattern and ref_tag
        candidates = [
            name for name, var in ds.data_vars.items()
            if (base_var in name and ref_tag in name)
        ]
        if not candidates:
            raise RuntimeError(
                f"No obs vars containing '{base_var}' and '{ref_tag}' "
                f"found in dataset variables."
            )
            
        chosen = candidates[0]          # e.g. sstStd_lon__ERA-Interim
        obs_da = ds[chosen].squeeze()
        # rename to drop the observation name 
        obs_da = obs_da.rename(base_var)
        
        ds.close()
        return obs_da

    # ------------------------------ public API ------------------------------
    def load(
        self,
        enso_group: str,
        var_name: str,
        period: str = "hist",
        nc_var: Optional[str] = None,
        ref_tag: Optional[str] = None,
        case_id: Optional[str] = None, 
    ) -> xr.DataArray:
        """
        Load a given ENSO diagnostic for all members of one period.

        Parameters
        ----------
        enso_group : {"ENSO_perf", "ENSO_proc", "ENSO_tel"}
        var_name   : variable key from self.enso_groups[enso_group], e.g. "enso_amplitude"
        period     : "hist" or "future"
        nc_var     : optional variable name inside the NetCDF file.
                     If None, the first data_var is used.

        Returns
        -------
        da_model, da_obs : xarray.DataArray
            With dims: member + (whatever dims the metric has).
            Coordinates include 'member' (int) and 'member_str'.
        """
        # Sanity checks
        if enso_group not in self.enso_groups:
            raise ValueError(
                f"Unknown ENSO group {enso_group}. "
                f"Must be one of {list(self.enso_groups.keys())}"
            )
        if var_name not in self.enso_groups[enso_group]:
            raise ValueError(
                f"var_name '{var_name}' not in enso_groups['{enso_group}']"
            )
        if var_name not in self.file_suffix_map:
            raise KeyError(
                f"No file suffix mapping for '{var_name}'. "
                f"Add it to the file_suffix_map."
            )

        suffix = self.file_suffix_map[var_name]
        member_dirs = self._list_member_dirs(period)

        da_model = {}
        da_obs = {}
        for i,mdir in enumerate(member_dirs):
            base = os.path.basename(mdir)  # e.g. v3.LR.historical_0051
            mstr = "{:02d}".format(i)
            
            nc_path = self._find_nc_file(mdir, enso_group, suffix, case_id)
            ds = xr.open_dataset(nc_path,decode_times=False)

            data_vars = list(ds.data_vars)
            # drop bounds_* helpers
            candidates = [v for v in data_vars if not v.startswith("bounds_")]

            if nc_var is not None:
                # nc_var is a BASE PATTERN like "sstStd_lon"
                # -> for model, require pattern + member id
                cand_pattern = [v for v in candidates if nc_var in v]
                cand_member  = [v for v in cand_pattern if mstr in v]

                if len(cand_member) >= 1:
                    data_var = cand_member[0]   # e.g. sstStd_lon__v3-LR_0051
                elif len(cand_pattern) >= 1:
                    # fallback: first match on pattern
                    data_var = cand_pattern[0]
                else:
                    # fallback to generic heuristic
                    data_var = self._choose_default_var(ds, candidates, mstr)
            else:
                data_var = self._choose_default_var(ds, candidates, mstr)

            # introduce member coordinate
            mod_da = ds[data_var].squeeze()

            # rename to drop the observation name 
            mod_da = mod_da.rename(nc_var)
        
            da_model[mstr] = mod_da
            
            # extract observation or reference vars 
            do = self._extract_obs(
                ds,
                base_var=nc_var,
                ref_tag=ref_tag
            )
            da_obs[mstr] = do
            
            ds.close()

        if not da_obs or not da_model:
            raise RuntimeError(
                f"No data loaded for {enso_group}/{var_name}/{period}"
            )
            
        return da_model, da_obs

    def load_metric_data(
        self,
        enso_group: str,
        var_name: str,
        nc_var: Optional[str] = None,
        ref_dict: Optional[dict] = None,
        period_list: Optional[list] = None,
        case_id: Optional[str] = None, 
    ) -> xr.Dataset:
        """
        Convenience wrapper: load one or more periods and return a Dataset
        with an extra 'period' dimension.

        Parameters
        ----------
        period_list : sequence of str, optional
            Period tags to pass to `self.load`, e.g. ["hist"] or
            ["hist", "future"]. If None, defaults to ["hist", "future"].

        Returns
        -------
        ds_model, ds_obs : xr.Dataset
            Both with variable:
            - 'metric' : dims (period, member, ...)
        """
        # Decide which period tags to load
        if period_list is None:
            periods = list(self.groups)
        else:
            periods = period_list

        dm_list = {}
        do_list = {}

        for per in periods:
            print(f"processing period: {per}")
            dm, do = self.load(
                enso_group, var_name, period=per, nc_var=nc_var, ref_tag=ref_dict[per], case_id=case_id,  
            )
            dm_list[per] = dm 
            do_list[per] = do 
            
        return dm_list, do_list
        
    def combine_members_to_array(self, member_dict, sample_dim=None):
        """
        Combine member DataArrays into a single DataArray with a new 'member' dim.
    
        Parameters
        ----------
        member_dict : dict
            Mapping like {'00': DataArray, '01': DataArray', ...}
    
        Returns
        -------
        da_model : xarray.DataArray
            dims: ('member', *original_dims)
            coords:
              - 'member'     : numeric or string member index
              - 'member_str' : string label for each member
              - all original coords preserved for non-'member' dims
        """
        if not member_dict:
            raise ValueError("member_dict is empty.")
    
        # Use first DataArray as template
        template = next(iter(member_dict.values()))
        template = template.squeeze()
    
        das = []
    
        for mem_key, da in member_dict.items():
            da = da.squeeze()
    
            # Optional: sanity check dims/sizes match template
            if da.dims != template.dims:
                raise ValueError(
                    f"Member '{mem_key}' has dims {da.dims}, expected {template.dims}"
                )
            for d in template.dims:
                if da.sizes[d] != template.sizes[d]:
                    raise ValueError(
                        f"Member '{mem_key}' has size {da.sizes[d]} on dim '{d}', "
                        f"expected {template.sizes[d]}"
                    )
    
            # Try to use an integer member index if possible, else keep as string
            try:
                mid = int(mem_key)
            except (TypeError, ValueError):
                mid = mem_key
    
            # Add member dimension and coords
            da = da.expand_dims({"member": [mid]})
            da = da.assign_coords(member=("member", [mid]))
            da = da.assign_coords(member_str=("member", [str(mem_key)]))
    
            das.append(da)
    
        # Concatenate along member
        da_model = xr.concat(das, dim="member")
    
        # If member is numeric, sort by member index
        if np.issubdtype(da_model["member"].dtype, np.number):
            order = np.argsort(da_model["member"].values)
            da_model = da_model.isel(member=order)
    
        return da_model

    def pool_members_to_samples(self, member_dict, sample_dim=None):
        """
        Pool member DataArrays along `sample_dim` into a single DataArray with
        a new 'sample' dimension, while keeping all other dims unchanged.
    
        Parameters
        ----------
        member_dict : dict
            Mapping like {'00': DataArray, '01': DataArray, ...}
        sample_dim : str or None, optional
            Name of the dimension to be treated as samples (e.g. 'years').
            If None, and the DataArray is 1D, that single dim is used.
    
        Returns
        -------
        pooled : xarray.DataArray
            dims: ('sample', *other_dims)
            coords:
              - 'sample': integer index
              - 'member': ('sample',) coord indicating which member each sample came from
              - all coords for other_dims copied from the first DataArray
        """
        if not member_dict:
            raise ValueError("member_dict is empty.")
    
        # Use first DataArray as template
        template = next(iter(member_dict.values()))
    
        # Infer sample_dim if not provided
        if sample_dim is None:
            if len(template.dims) == 1:
                sample_dim = template.dims[0]
            else:
                raise ValueError(
                    f"sample_dim is None, but DataArray has multiple dims {template.dims}. "
                    "Please specify sample_dim explicitly."
                )
    
        if sample_dim not in template.dims:
            raise ValueError(
                f"sample_dim='{sample_dim}' not found in dims {template.dims}"
            )
    
        # All dims, with sample_dim excluded
        other_dims = [d for d in template.dims if d != sample_dim]
    
        # Check that non-sample dims are consistent across members
        template_other_sizes = {d: template.sizes[d] for d in other_dims}
        for mem, da in member_dict.items():
            if set(da.dims) != set(template.dims):
                raise ValueError(f"Member {mem} has different dims: {da.dims}")
            for d in other_dims:
                if da.sizes[d] != template_other_sizes[d]:
                    raise ValueError(
                        f"Member {mem} has different size along dim '{d}': "
                        f"{da.sizes[d]} vs {template_other_sizes[d]}"
                    )
    
        vals = []
        members = []
    
        for mem, da in member_dict.items():
            # Put sample_dim first to make concatenation easy
            da_t = da.transpose(sample_dim, *other_dims)
            v = da_t.values  # shape: (n_sample, *other_shapes)
            n = v.shape[0]
            vals.append(v)
            members.extend([mem] * n)
    
        # Concatenate along sample axis (axis=0)
        all_vals = np.concatenate(vals, axis=0)
    
        # Build coords for other dims from template
        other_coords = {
            d: template.coords[d] for d in other_dims if d in template.coords
        }
    
        pooled = xr.DataArray(
            all_vals,
            dims=("sample", *other_dims),
            coords={
                "sample": np.arange(all_vals.shape[0]),
                "member": ("sample", members),
                **other_coords,
            },
            name=template.name,
            attrs=template.attrs,
        )
        return pooled

    def validate_constant_observation(
        self,
        obs_dict,
        ref_group="hist",
        ref_member="00",
        sample_dim=None,
        use_allclose=False,
        rtol=1e-8,
        atol=0.0,
        pool_ensemble=True,
        
    ):
        """
        Validate that observation DataArrays are identical across all periods/members,
        then return a DataArray shaped consistently with the model by pooling along
        the specified sample dimension.
    
        Parameters
        ----------
        self : object
            Object that provides `pool_members_to_samples` as a method.
            If you use this as a standalone function, remove `self` and call the
            function version of `pool_members_to_samples` instead.
        obs_dict : dict
            Nested dict: obs_dict[period][member] -> xarray.DataArray.
        ref_group : str, optional
            Period key used for the reference, e.g. 'hist'.
        ref_member : str, optional
            Member key used for the reference, e.g. '00'.
        sample_dim : str or None, optional
            Dimension to treat as samples when pooling.
            - If None and the DataArray is 1D, the single dim is inferred.
            - If None and the DataArray has multiple dims, `pool_members_to_samples`
              will raise a ValueError asking for an explicit `sample_dim`.
        use_allclose : bool, optional
            If True, use np.allclose instead of np.array_equal.
        rtol, atol : float, optional
            Tolerances passed to np.allclose if use_allclose=True.
    
        Returns
        -------
        out : xarray.DataArray
            Observation DataArray with dims ('sample', *other_dims), created via
            `pool_members_to_samples` so that its shape/dim ordering is consistent
            with the pooled model arrays.
    
        Raises
        ------
        KeyError
            If the reference group/member is missing.
        ValueError
            If any entry differs from the reference.
        """
        # --- get reference ---
        if ref_group not in obs_dict:
            raise KeyError(
                f"Reference group '{ref_group}' not found in obs_dict keys: "
                f"{list(obs_dict.keys())}"
            )
    
        if ref_member not in obs_dict[ref_group]:
            raise KeyError(
                f"Reference member '{ref_member}' not found in obs_dict['{ref_group}'] "
                f"keys: {list(obs_dict[ref_group].keys())}"
            )
    
        ref = obs_dict[ref_group][ref_member]
    
        # --- validate all obs equal to ref (by values) ---
        for period, members in obs_dict.items():
            for member, dx in members.items():
                if use_allclose:
                    same = np.allclose(
                        dx.values, ref.values,
                        rtol=rtol, atol=atol, equal_nan=True
                    )
                else:
                    same = np.array_equal(dx.values, ref.values)
    
                if not same:
                    raise ValueError(
                        f"Observation `obs_dict` differs at period={period}, member={member} "
                        f"(values are not identical to reference {ref_group}/{ref_member})."
                    )
    
        # --- make sure returned obs has pooled/sample shape consistent with model ---
        # Use a fake 1-member dict so we can reuse pool_members_to_samples logic
        if pool_ensemble:
            out = self.pool_members_to_samples({"00": ref}, sample_dim=sample_dim)
        else:
            out = ref 
    
        return out