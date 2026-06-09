#!/usr/bin/env python
from pathlib import Path
from statistics import median
import pysam
import polars as pl
import subprocess

repo_root = Path("..").resolve()


## Functions
def calculate_median_insert_size(bam_path, max_reads=2_000_000, threads=1):
	"""Compute median insert size from proper pairs.
	
	Samples up to max_reads proper pairs from the start of the BAM.
	Raises ValueError if no valid pairs found.
	"""
	sizes = []
	n_zero_tlen = 0
	
	with pysam.AlignmentFile(bam_path, "rb", threads=threads) as bam:
		for read in bam:
			if len(sizes) >= max_reads:
				break
			if (
				read.is_proper_pair          # both mates aligned sensibly together
				and not read.is_secondary    # not an alternate mapping
				and not read.is_supplementary # not a split-read fragment
				and not read.is_duplicate    # not a PCR/optical duplicate
				and not read.is_unmapped     # actually aligned somewhere
			):
				if read.template_length == 0:
					n_zero_tlen += 1
					continue
				sizes.append(abs(read.template_length))
	
	if n_zero_tlen:
		print(
			f"Skipped {n_zero_tlen} proper-pair reads with TLEN=0 in {bam_path}"
		)
	
	if not sizes:
		raise ValueError(f"No valid proper pairs found in {bam_path}")
	
	return int(median(sizes))

def calculate_median_coverage(
	bam_path,
	bed_path,
	threads=2,
	min_mapq=20,
	tmp_outdir="mosdepth",
	skip_if_exists=True,
):
	"""
	Run mosdepth on a BAM restricted to regions in `bed_path` and return the
	median per-region coverage from the total distribution.

	Parameters
	----------
	bam_path : Path
		Input BAM file.
	bed_path : Path
		BED file of regions to evaluate.
	threads : int
		Threads for mosdepth (`--threads`).
	min_mapq : int
		Minimum mapping quality (`--mapq`).
	tmp_outdir : str | Path
		Directory to write mosdepth outputs into.
	skip_if_exists : bool
		If True (default), skip running mosdepth when the expected
		`*.mosdepth.region.dist.txt` already exists. Set False to force re-run.

	Returns
	-------
	int
		Median region coverage (depth at which cumulative proportion >= 0.5).
	"""
	bam_path = Path(bam_path)
	bed_path = Path(bed_path)
	tmp_outdir = Path(tmp_outdir)
	tmp_outdir.mkdir(exist_ok=True, parents=True)

	sample_name = bam_path.name.split(".")[0]
	prefix = tmp_outdir / sample_name
	dist_file = tmp_outdir / f"{sample_name}.mosdepth.region.dist.txt"

	# ---- Run mosdepth (or skip if output already exists) ----
	if skip_if_exists and dist_file.exists():
		print(f"[mosdepth] Found existing {dist_file}, skipping run.")
	else:
		cmd = [
			"mosdepth",
			"--by", str(bed_path),
			"--no-per-base",
			"--threads", str(threads),
			"--mapq", str(min_mapq),
			str(prefix),
			str(bam_path),
		]
		print(f"[mosdepth] Running: {' '.join(cmd)}")
		try:
			subprocess.run(cmd, check=True)
		except FileNotFoundError as e:
			raise RuntimeError(
				"mosdepth executable not found in PATH. Install it or activate the right env."
			) from e
		except subprocess.CalledProcessError as e:
			raise RuntimeError(
				f"mosdepth failed with exit code {e.returncode} for {bam_path}"
			) from e

	if not dist_file.exists():
		raise FileNotFoundError(
			f"Expected mosdepth output not found: {dist_file}"
		)

	# ---- Parse distribution and compute median ----
	mosdepth_dist = pl.read_csv(
		dist_file,
		separator="\t",
		has_header=False,
		new_columns=["chrom", "depth", "proportion"],
	)
	median = (
		mosdepth_dist
		.filter((pl.col("chrom") == "total") & (pl.col("proportion") >= 0.5))
		.sort("depth", descending=True)
		.select("depth")
		.head(1)
		.item()
	)
	return median

## Setup
### Inputs
bed_path = repo_root / "data/regions/sureselect-all-exon-v5_hg38_regions_200bp-pad.bed"
ffperase_launcher = repo_root / "common-ffpe-snvf/templates/ffpe-snvf/ffperase.sh"
ref_path = repo_root / "data/ref/Homo_sapiens_assembly38.fasta"

vcf_paths = sorted([path.resolve() for path in list(Path(repo_root / "vcf/EGAD00001004066/somatic_filtered").glob("*/*FFPE*.vcf.gz")) if "FFPE" in str(path)])

### Outputs
script_outdir = repo_root / "ffpe-snvf/script_ffperase"
script_outdir.mkdir(exist_ok=True, parents=True)


## Create Execution Scripts
for i, vcf_path in enumerate(vcf_paths, start=1):

	sample_name = vcf_path.parent.name
	bam_path = repo_root / "data/EGAD00001004066/bam" / sample_name / f"{sample_name}.bam"

	print(f"{i}. Creating scripts for sample: {sample_name}")

	if not bam_path.exists():
		print(f"\tBAM does not exist: {bam_path}")
		continue

	ins_size = calculate_median_insert_size(bam_path)
	coverage = calculate_median_coverage(bam_path, bed_path)

	exec_script_outpath = script_outdir / f"{sample_name}_ffperase.sh"

	result_outdir = repo_root / f"ffpe-snvf/EGAD00001004066/somatic_filtered/ffperase/{sample_name}"
	result_outdir.mkdir(exist_ok=True, parents=True)

	contents = f"""#!/usr/bin/env bash
set -euo pipefail

bash {ffperase_launcher} \\
	--vcf           {vcf_path} \\
	--bam           {bam_path} \\
	--reference     {ref_path} \\
	--bed           {bed_path} \\
	--coverage      {coverage} \\
	--median-insert {ins_size} \\
	--sample-name   {sample_name} \\
	--step          full \\
	--mutation-type snvs \\
	--outdir {result_outdir}
"""

	with open(exec_script_outpath, "w") as f:
		f.write(contents)


	print(f"Wrote {exec_script_outpath}  (coverage={coverage}, insert={ins_size})")
	
	


