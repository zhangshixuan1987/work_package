"""Computational cost of alternating FC / FOSI spin-up segments (../computation/Simulation_Cost_Analysis.ipynb)."""
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


@dataclass
class SimulationRun:
    """
    Represents a simulation run configuration.

    Attributes:
        code: Code base used for the run (e.g., "v3.LR").
        exp: Experiment name or tag.
        case_types: List of case type descriptors (e.g., ["DATM", "CPLHIST"]).
        nodes: List of node counts used for different test cases.
        sypd: List of simulated years per day (performance), one per node setting.
    """
    code: str
    exp: str
    case_types: List[str]  # Allow multiple case types
    nodes: List[int]
    sypd: List[Optional[float]]  # May include None if performance not available


class AlternatingSimulationEvaluator:
    """
    Cost of a spin-up that alternates ``len_b`` years of fully coupled (FC, B-case) with
    ``len_g`` years of forced ocean–sea ice (FOSI, G-case), for total lengths ``len_t``.
    """

    def __init__(self, len_b, len_g, len_t):
        self.len_b = len_b
        self.len_g = len_g
        self.len_t = len_t
        plt.rcParams.update({
            'axes.titlesize': 16,
            'axes.labelsize': 14,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 12,
            'legend.title_fontsize': 14
        })
        self.textfontsize = 9
        self.textcolors = ["red", "black"]
        shape = (len(len_b), len(len_g), len(len_t))
        self.total_run_years = np.zeros(shape)
        self.total_cost_months = np.zeros(shape)
        self.cost_save_ratio = np.zeros(shape)

    def calculate_cost(self, sypd, nodes, reference_nodes=104, total_simulation_years=2000):
        """Mean SYPD, wallclock months per simulated year (normalized to ``reference_nodes``), total months."""
        avg_sypd = np.mean([s for s in sypd if s is not None])
        avg_nodes = np.mean(nodes)
        normalized_sypd = avg_sypd * (reference_nodes / avg_nodes)
        days_per_sim_year = 1.0 / normalized_sypd
        months_per_year = days_per_sim_year / 30.44
        cost_months = total_simulation_years * months_per_year
        return avg_sypd, months_per_year, cost_months

    def evaluate_costs_months(self, cost_b, cost_g):
        for i, b in enumerate(self.len_b):
            for j, g in enumerate(self.len_g):
                for k, t in enumerate(self.len_t):
                    if b + g == 0 or t == 0:
                        self.total_run_years[i, j, k] = 0
                        self.total_cost_months[i, j, k] = 0
                        continue

                    cycle_len = b + g
                    n_cycles = t // cycle_len
                    rem = t % cycle_len

                    run_b = n_cycles * b
                    run_g = n_cycles * g

                    if rem > 0:
                        rem_b = min(rem, b)
                        rem_g = rem - rem_b
                        run_b += rem_b
                        run_g += rem_g

                    total_years = run_b + run_g
                    total_months = run_b * cost_b + run_g * cost_g

                    self.total_run_years[i, j, k] = total_years
                    self.total_cost_months[i, j, k] = total_months

        return self.total_run_years, self.total_cost_months

    def evaluate_costs_ratio_to_B(self, cost_b, cost_g):
        for i, b in enumerate(self.len_b):
            for j, g in enumerate(self.len_g):
                for k, t in enumerate(self.len_t):
                    if b + g == 0 or t == 0:
                        self.cost_save_ratio[i, j, k] = 0
                        continue
                    cycle_len = b + g
                    n_cycles = t // cycle_len
                    rem = t % cycle_len
                    alt_cost = n_cycles * (cost_b * b + cost_g * g)
                    if rem > 0:
                        rem_b = min(rem, b)
                        rem_g = rem - rem_b
                        alt_cost += cost_b * rem_b + cost_g * rem_g
                    base_cost = cost_b * t
                    self.cost_save_ratio[i, j, k] = 100 * (base_cost - alt_cost) / base_cost

    def plot_cost_overlay(self, cost_data, x_labels, y_labels, title, colorbar_label, fig_name, overlay_text=None):
        fig, ax = plt.subplots(figsize=(10, 8))

        cmap = plt.get_cmap("YlGn")
        norm = plt.Normalize(vmin=np.nanmin(cost_data), vmax=np.nanmax(cost_data))

        heatmap = ax.imshow(cost_data[::-1, :], cmap=cmap, norm=norm, alpha=0.8)

        ax.set_xticks(np.arange(len(x_labels)))
        ax.set_yticks(np.arange(len(y_labels)))
        ax.set_xticklabels(x_labels)
        ax.set_yticklabels(y_labels[::-1])

        ax.set_xticks(np.arange(len(x_labels) + 1) - 0.5, minor=True)
        ax.set_yticks(np.arange(len(y_labels) + 1) - 0.5, minor=True)
        ax.grid(which="minor", color="black", linestyle='-', linewidth=0.5)
        ax.tick_params(which="minor", bottom=False, left=False)

        overlay_text = overlay_text[::-1, :]
        for i in range(len(y_labels)):
            for j in range(len(x_labels)):
                label = overlay_text[i, j] if overlay_text is not None else f"{cost_data[::-1, :][i, j]:.1f}"
                text_color = self.textcolors[0] if x_labels[j] == 0 or y_labels[::-1][i] == 0 else self.textcolors[1]
                ax.text(j, i, label, ha='center', va='center', fontweight='bold', color=text_color, fontsize=self.textfontsize)

        ax.set_xlabel("FOSI Segment Length (years)")
        ax.set_ylabel("FC Segment Length (years)")
        ax.set_title(title)

        cbar = plt.colorbar(heatmap, ax=ax, label=colorbar_label, fraction=0.035, pad=0.04)
        plt.tight_layout()
        plt.savefig(fig_name)
        plt.close()

    def plot_run_years_fill_cost_overlay(
            self,
            fixed_t_idx=None,
            fig_name="run_years_fill_cost_overlay.pdf",
            fig_title=None,
            cb_title=None,
        ):
        if fixed_t_idx is None:
            fixed_t_idx = len(self.len_t) - 1

        t_val = self.len_t[fixed_t_idx]
        if t_val == 0:
            print("Total simulation length is 0; skipping plot.")
            return

        run_years_map = self.total_run_years[:, :, fixed_t_idx]
        cost_months_map = self.total_cost_months[:, :, fixed_t_idx]

        labels = np.empty_like(run_years_map, dtype=object)
        for i in range(labels.shape[0]):
            for j in range(labels.shape[1]):
                sim_len = run_years_map[i, j]
                cost_per_year = cost_months_map[i, j] / t_val
                total_cost = sim_len * cost_per_year
                if self.len_b[i] + self.len_g[j] == 0:
                    labels[i, j] = "NA"
                else:
                    labels[i, j] = f"{total_cost:.2f}m" if total_cost > 0 else "NA"

        if fig_title is None:
            fig_title = "Overlay: Fill = Cost (months), Label = Cost for B+G Run"

        if cb_title is None:
            cb_title = f"Estimated Cost for {t_val}-yr Ocean-Sea Ice Simulation (months)"

        self.plot_cost_overlay(
            cost_data=cost_months_map,
            x_labels=self.len_g,
            y_labels=self.len_b,
            title=fig_title,
            colorbar_label=cb_title,
            fig_name=fig_name,
            overlay_text=labels
        )

    def plot_run_years_fill_cost_savings_overlay(
            self, cost_b, cost_g,
            fixed_t_idx=None,
            fig_name="cost_savings_overlay.pdf",
            fig_title=None,
            cb_title=None,
        ):
        if fixed_t_idx is None:
            fixed_t_idx = len(self.len_t) - 1

        t_val = self.len_t[fixed_t_idx]
        if t_val == 0:
            print("Total simulation length is 0; skipping plot.")
            return

        run_years_map = self.total_run_years[:, :, fixed_t_idx]
        cost_months_map = self.total_cost_months[:, :, fixed_t_idx]
        labels = np.empty_like(run_years_map, dtype=object)

        base_cost = cost_b * t_val
        cost_data = cost_months_map

        for j in range(labels.shape[1]):
            for i in range(labels.shape[0]):
                alt_cost = cost_months_map[i, j]
                if self.len_b[i] + self.len_g[j] == 0 or base_cost == 0:
                    labels[i, j] = "NA"
                else:
                    savings_pct = 100 * (base_cost - alt_cost) / base_cost
                    labels[i, j] = f"{savings_pct:.1f}%\n{alt_cost:.1f}m"

        if fig_title is None:
            fig_title = "Overlay: Fill = Total Cost (months), Label= %Saving vs FC-only and total wallclock time (months)"

        if cb_title is None:
            cb_title = f"Total Cost of {t_val}-yr Ocean-Sea Ice Simulation (months)"

        self.plot_cost_overlay(
            cost_data=cost_data,
            x_labels=self.len_g,
            y_labels=self.len_b,
            title=fig_title,
            colorbar_label=cb_title,
            fig_name=fig_name,
            overlay_text=labels
        )

    def generate_case_summary(self, sim_runs, reference_nodes=104, total_simulation_years=2000):
        summary_rows = []
        for run in sim_runs:
            for case_type in run.case_types:
                avg_sypd, mpyr, months = self.calculate_cost(
                    sypd=run.sypd,
                    nodes=run.nodes,
                    reference_nodes=reference_nodes,
                    total_simulation_years=total_simulation_years
                )
                summary_rows.append({
                    "code": run.code,
                    "exp": run.exp,
                    "type": case_type,
                    "avg_sypd": round(avg_sypd, 2),
                    "months_per_year": round(mpyr, 6),   # ~1e-3; 2 decimals rounded it to 0.0
                    "total_months": round(months, 2)
                })
        return pd.DataFrame(summary_rows)
