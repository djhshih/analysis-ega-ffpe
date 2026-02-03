#!/usr/bin/env python
import polars as pl
import glob
import os
import sys
import numpy as np
from tqdm import tqdm
from scipy.spatial.distance import jensenshannon
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Union, Tuple

## Local Dependencies
sys.path.append("../common-ffpe-snvf/python")
from common import fdr_cut_pred, ct_filter, natural_sort_variants, snv_filter
from FFPEsig import ffpe_sig_unrepaired, correct_FFPE_profile
from mutation_signatures import annotate_variant_context, build_mutation_count_matrix, plot_sbs96_signature


## Functions
def js_div_violin(
    data: pl.DataFrame,
    x_col: str = "sample_type", 
    y_col: str = "js_divergence", 
    hue_col: Optional[str] = "sample_type", 
    x_label: str = "Sample Type",
    y_label: str = "Jensen-Shannon Divergence",
    title: str = "JS Divergence vs Ground Truth",
    rotate: Optional[int] = None, 
    file: Optional[str] = None,
    # Visual Customization Arguments
    figsize: Tuple[float, float] = (5, 5),
    cut: float = 0,
    palette: str = "muted",
    jitter: float = 0.3,
    point_size: float = 5,
    point_alpha: float = 0.3,
    point_color: str = "black"
) -> None:
    """
    Generates a violin plot with an overlaid strip plot for JS Divergence.

    Parameters:
    -----------
    data : pl.DataFrame
        Input DataFrame containing the data to plot.
    x_col : str, default "sample_type"
        Column name for the x-axis (categorical).
    y_col : str, default "js_divergence"
        Column name for the y-axis (numerical).
    hue_col : str, optional, default "sample_type"
        Column name for color grouping. Set to None to disable hue nesting.
    x_label : str, default "Sample Type"
        Label for the x-axis.
    y_label : str, default "Jensen-Shannon Divergence"
        Label for the y-axis.
    title : str, default "JS Divergence vs Ground Truth"
        Main title of the plot.
    rotate : int, optional
        Degrees to rotate x-axis tick labels.
    file : str, optional
        Path to save the figure. If None, displays the plot.
    figsize : Tuple[float, float], default (5, 5)
        Width and height of the figure in inches.
    cut : float, default 0
        Distance, in units of bandwidth, to extend the density past the extreme 
        datapoints. Set to 0 to limit the violin range to the observed data range.
    palette : str, default "muted"
        Seaborn color palette name.
    jitter : float, default 0.3
        Amount of jitter (noise) to apply to the strip plot points along the categorical axis.
    point_size : float, default 5
        Radius of the strip plot points.
    point_alpha : float, default 0.3
        Transparency level of the strip plot points (0 to 1).
    point_color : str, default "black"
        Color of the strip plot points.
    """
    
    # Set context only for this plot to avoid affecting global settings
    with sns.plotting_context("paper", font_scale=1.5), sns.axes_style("ticks"):
        
        # Initialize figure inside the context
        fig, ax = plt.subplots(figsize=figsize)
        
        # Violin Plot
        sns.violinplot(
            data=data, 
            x=x_col, 
            y=y_col, 
            hue=hue_col, 
            cut=cut, 
            legend=False,
            ax=ax,
            palette=palette
        )
        
        # Strip Plot (Points)
        sns.stripplot(
            data=data, 
            x=x_col, 
            y=y_col, 
            jitter=jitter, 
            size=point_size,
            color=point_color, 
            alpha=point_alpha, 
            ax=ax
        )
        
        # Set labels using the Axes object
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title, pad=20)
        
        # Handle Rotation
        if rotate is not None:
            # FixedLocator is required in newer matplotlib versions before setting labels
            ax.set_xticks(ax.get_xticks()) 
            ax.set_xticklabels(
                ax.get_xticklabels(), 
                rotation=rotate, 
                ha='right', 
                rotation_mode="anchor"
            )

        # Clean up borders (top and right spines)
        sns.despine()

        plt.tight_layout()

        if file:
            plt.savefig(file, dpi=300, bbox_inches="tight")
            print(f"Violin plot saved to {file}")
            plt.close()
        else:
            plt.show()


def js_div_dot(
    data: Union[pl.DataFrame, pd.DataFrame],
    x_col: str = "sample_type", 
    y_col: str = "js_divergence", 
    hue_col: Optional[str] = "sample_type",
    x_label: str = "Sample Type",
    y_label: str = "Jensen-Shannon Divergence",
    title: str = "JS Divergence vs Ground Truth",
    rotate: Optional[int] = None, 
    file: Optional[str] = None,
    # Visual Customization Arguments
    figsize: Tuple[float, float] = (5, 5),
    palette: str = "deep",
    point_size: float = 8,
    point_alpha: float = 0.7,
    point_jitter: float = 0.2,
    linewidth: float = 0,
    edgecolor: str = "black"
) -> None:
    """
    Generates a strip (dot) plot for JS Divergence.
    Accepts Polars or Pandas DataFrames.

    Parameters:
    -----------
    data : pl.DataFrame or pd.DataFrame
        Input DataFrame containing the data to plot.
    x_col : str, default "sample_type"
        Column name for the x-axis.
    y_col : str, default "js_divergence"
        Column name for the y-axis.
    hue_col : str, optional, default "sample_type"
        Column name for color grouping.
    x_label : str, default "Sample Type"
        Label for the x-axis.
    y_label : str, default "Jensen-Shannon Divergence"
        Label for the y-axis.
    title : str, default "JS Divergence vs Ground Truth"
        Main title of the plot.
    rotate : int, optional
        Degrees to rotate x-axis tick labels.
    file : str, optional
        Path to save the figure. If None, displays the plot.
    figsize : Tuple[float, float], default (5, 5)
        Width and height of the figure in inches.
    palette : str, default "deep"
        Seaborn color palette name.
    point_size : float, default 8
        Radius of the points.
    point_alpha : float, default 0.7
        Transparency level of the points (0 to 1).
    point_jitter : float, default 0.2
        Amount of horizontal jitter to separate overlapping points.
    linewidth : float, default 0
        Width of the border around the points. Set to 1 for a defined border.
    edgecolor : str, default "grey"
        Color of the border around the points (if linewidth > 0).
    """

    # Set context only for this plot
    with sns.plotting_context("paper", font_scale=1.5), sns.axes_style("ticks"):
        
        # Initialize figure inside the context
        fig, ax = plt.subplots(figsize=figsize)
        
        # Strip Plot
        sns.stripplot(
            data=data, 
            x=x_col, 
            y=y_col,
            hue=hue_col,
            size=point_size,
            jitter=point_jitter, 
            alpha=point_alpha, 
            linewidth=linewidth,
            edgecolor=edgecolor,
            palette=palette,
            legend=False, # Usually False for strip plots where x=hue
            ax=ax
        )
        
        # Set labels using the Axes object
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title, pad=20)
        
        # Handle Rotation
        if rotate is not None:
            ax.set_xticks(ax.get_xticks())
            ax.set_xticklabels(
                ax.get_xticklabels(), 
                rotation=rotate, 
                ha='right', 
                rotation_mode="anchor"
            )

        # Clean up borders
        sns.despine()

        plt.tight_layout()

        if file:
            plt.savefig(file, dpi=300, bbox_inches="tight")
            print(f"Dot Plot saved to {file}")
            plt.close()
        else:
            plt.show()


def get_ground_truth(orig: pl.DataFrame, truth_set: pl.DataFrame) -> pl.DataFrame:
	ct_snv = orig.pipe(ct_filter)
	nct_snv = orig.pipe(lambda df : ct_filter(df, invert=True))
	true_ct_snv = ct_snv.join(truth_set, on = ["chrom", "pos", "ref", "alt"], how ="semi")
	ground_truth = pl.concat([true_ct_snv, nct_snv]).pipe(natural_sort_variants)
	return ground_truth

def trinuc_context_count(snvs: pl.DataFrame, ref_fasta_path: str) -> np.ndarray: 
	return build_mutation_count_matrix(annotate_variant_context(snvs, ref_fasta_path))["count"].to_numpy()

def filter_mobsnvf(mobsnvf_res: pl.DataFrame, fp_cut: float = 1e-08) -> pl.DataFrame:
	return fdr_cut_pred(mobsnvf_res, score_col="FOBP", fp_cut=fp_cut).filter(pl.col("pred"))

def ffpesig_correct(trinuc_context_count: np.ndarray) -> np.ndarray:
	corrected_profile, corrected_solutions = correct_FFPE_profile(
		V = trinuc_context_count,
		W1 = ffpe_sig_unrepaired
	)
	return corrected_profile

## Setup
ref_fasta_path = "../data/ref/Homo_sapiens_assembly38.fasta"
ffpe_snvf_root = "../ffpe-snvf"
fp_cut = 1e-10

## Convert Wide to long
sample_type_map = {
	'ffpe.js_divergence': "FFPE\nSample", 
	'mobsnvf.js_divergence': "MOBSNVF\nFiltered", 
	'ffpesig.js_divergence': "FFPESig\nCorrected"
}

def process_dataset(dataset: str, variant_set: str, annot_path: str) -> None:

	sample_paths = sorted(glob.glob(f"{ffpe_snvf_root}/{dataset}/{variant_set}/gatk-obmm/*/*.gatk-obmm.tsv"))

	js_divergence = []

	for path in tqdm(sample_paths):
		sample_name = path.split("/")[-2]

		annot = pl.read_csv(annot_path, separator="\t")
		case_id = annot.filter(pl.col("sample_name") == sample_name)[0, "case_id"]

		ffpe = pl.read_csv(path, separator="\t").pipe(snv_filter)

		truth_path = glob.glob(f"{dataset}/{variant_set}/truth_sets/*/{case_id}.tsv")[0]
		truth_set = pl.read_csv(truth_path, separator="\t")
		ground_truth = get_ground_truth(ffpe, truth_set)

		mobsnvf_path = f"{ffpe_snvf_root}/{dataset}/{variant_set}/mobsnvf/{sample_name}/{sample_name}.mobsnvf.snv"
		mobsnvf = pl.read_csv(mobsnvf_path, separator="\t").pipe(lambda df: filter_mobsnvf(df, fp_cut = fp_cut))
		
		spectrum_outdir = f"{dataset}/{variant_set}/plots/mutation-spectrum_fp-cut-{fp_cut:.0e}/{sample_name}"
		os.makedirs(spectrum_outdir, exist_ok=True)
		
		orig_tncc = trinuc_context_count(ffpe, ref_fasta_path)
		truth_tncc = trinuc_context_count(ground_truth, ref_fasta_path)
		mobsnvf_tncc = trinuc_context_count(mobsnvf, ref_fasta_path)
		ffpesig_tncc = ffpesig_correct(trinuc_context_count(ffpe, ref_fasta_path))
		y_lim = orig_tncc.max()

		# Generate Individual Plots
		plot_sbs96_signature(orig_tncc, base_fontsize=6, ylim=y_lim, label="FFPE", file=f"{spectrum_outdir}/{sample_name}_ffpe.pdf")
		plot_sbs96_signature(truth_tncc, base_fontsize=6, ylim=y_lim, label="Ground Truth", file=f"{spectrum_outdir}/{sample_name}_ground-truth.pdf")
		plot_sbs96_signature(mobsnvf_tncc, base_fontsize=6, ylim=y_lim, label = "MOBSNVF", file=f"{spectrum_outdir}/{sample_name}_mobsnvf-filtered.pdf")
		plot_sbs96_signature(ffpesig_tncc, base_fontsize=6, ylim=y_lim, label = "FFPESig", file=f"{spectrum_outdir}/{sample_name}_ffpesig-corrected.pdf")


		# Generate Combined Stacked Plot
		fig, axes = plt.subplots(4, 1, figsize=(12, 8)) 
		plot_sbs96_signature(orig_tncc, base_fontsize=6, ylim=y_lim, label="FFPE", ax=axes[0])
		plot_sbs96_signature(ffpesig_tncc, base_fontsize=6, ylim=y_lim, label="FFPESig", ax=axes[1])
		plot_sbs96_signature(mobsnvf_tncc, base_fontsize=6, ylim=y_lim, label="MOBSNVF", ax=axes[2])
		plot_sbs96_signature(truth_tncc, base_fontsize=6, ylim=y_lim, label="Ground Truth", ax=axes[3])
		plt.tight_layout()
		plt.savefig(f"{spectrum_outdir}/{sample_name}_combined_panel.pdf", bbox_inches="tight", dpi=300)
		plt.close(fig)

		# Jensen Shannon Divergence Calculation
		# Base = 2 keeps range between 0 and 1.
		per_sample_jsdiv = {
			"dataset" : f"ENA_{dataset}",
			"sample_name" : sample_name,
			"ffpe.js_divergence" : jensenshannon(orig_tncc, truth_tncc, base=2)**2,
			"mobsnvf.js_divergence" : jensenshannon(mobsnvf_tncc, truth_tncc, base=2)**2,
			"ffpesig.js_divergence" :  jensenshannon(ffpesig_tncc, truth_tncc, base=2)**2,
		}

		js_divergence.append(per_sample_jsdiv)

	js_outdir = f"{variant_set}/js-divergence"
	os.makedirs(js_outdir, exist_ok=True)
	js_divergence_df = pl.DataFrame(js_divergence)
	js_divergence_df.write_csv(f"{js_outdir}/jensen-shannon-divergence_mobsnvf-fp-cut-{fp_cut}.tsv", separator="\t")

	## Prepare plotting data
	js_divergence_long = (
		js_divergence_df.unpivot(
			index=["dataset", "sample_name"],
			on = ['ffpe.js_divergence', 'mobsnvf.js_divergence', 'ffpesig.js_divergence'],
			variable_name="sample_type",
			value_name="js_divergence"
		)
		.with_columns(pl.col("sample_type").map_elements(lambda x: sample_type_map.get(x, x), return_dtype=str))
	)

	## Plot Jensen-shannon divergence comparison
	js_div_violin(js_divergence_long, file=f"{js_outdir}/jansen-shannon-divergence_mobsnvf-fp-cut-{fp_cut:.0e}_mt_violin-plot.pdf")
	js_div_dot(js_divergence_long, file=f"{js_outdir}/jansen-shannon-divergence_mobsnvf-fp-cut-{fp_cut:.0e}_mt_dot-plot.pdf")

#----------------

process_dataset(
    dataset = "EGAD00001004066",
	variant_set = "somatic_filtered-dp10-blacklist_micr1234-excluded",
    annot_path = "../annot/EGAD00001004066/sample_annotations.tsv"
)

