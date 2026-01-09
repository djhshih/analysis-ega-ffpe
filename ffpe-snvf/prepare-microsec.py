#!/usr/bin/env python
import polars as pl
import os
import glob
import pysam

# %% [markdown]
# ### Helper Functions

def return_path_if_exists(path: str, abs=False) -> str:
	if os.path.exists(path):
		return os.path.abspath(path) if abs else path
	else:
		raise FileNotFoundError(f"File not found: {path}")

def read_variants(vcf_path: str, split_multiallelic = True) -> pl.DataFrame:
	variants = (
		pl.read_csv(
			vcf_path, 
			separator="\t", 
			comment_prefix="##", 
			null_values=["."],
			columns=["#CHROM", "POS", "REF", "ALT"], 
			infer_schema_length=100
		)
		.rename(lambda x: x.strip("#").lower())
		.with_columns(pl.col("chrom").str.strip_chars())
	)
 
	if split_multiallelic:
		variants = variants.with_columns(
			pl.col("alt").str.split(",").alias("alt")
		).explode("alt")
	
	return variants

def snv_filter(df: pl.DataFrame) -> pl.DataFrame:
	return df.filter(
		(pl.col("ref").str.len_chars() == 1) & 
		(pl.col("alt").str.len_chars() == 1)
	)

def ct_filter(df: pl.DataFrame) -> pl.DataFrame:
	return df.filter(
		((pl.col("ref") == "C") & (pl.col("alt") == "T")) | 
  		((pl.col("ref") == "G") & (pl.col("alt") == "A"))
	)

# Define a function to be applied to each row-like struct
# Used to retrieve neighbourhood sequence of 20bp in each direction via map_elements which is vectorized
def get_neighborhood_sequence(row_struct: dict, genome: pysam.FastaFile) -> str:
    chrom = row_struct["Chr"]
    pos = row_struct["Pos"] # 1-based VCF position
    ref = row_struct["Ref"]
    alt = row_struct["Alt"]

    # Pysam fetch is 0-based, half-open [start, end)
    # VCF Pos 100 = Index 99
    
    # Left Flank (20bp before the variant start)
    # If Pos=100, we want indices 79 to 99 (length 20)
    left_end_idx = pos - 1
    left_start_idx = max(0, left_end_idx - 20)
    
    # Right Flank (20bp after the REF allele ends)
    # If SNV (Ref len 1) at Pos 100: Ref is index 99. Next base is index 100.
    # If Del (Ref len 5) at Pos 100: Ref is 99..103. Next base is index 104.
    right_start_idx = left_end_idx + len(ref)
    right_end_idx = right_start_idx + 20

    try:
        left_flank = genome.fetch(chrom, left_start_idx, left_end_idx)
        right_flank = genome.fetch(chrom, right_start_idx, right_end_idx)

        full_seq = left_flank + alt + right_flank
        
        # Explicitly removing newline characters (\n) and carriage returns (\r)
        return full_seq.replace("\n", "").replace("\r", "").upper()

    except (ValueError, KeyError):
        # KeyError handles missing chromosomes
        return None

def get_read_length(bam_path: str) -> int:
	bam = pysam.AlignmentFile(bam_path, "rb")
	max_reads = 100000
	
	lengths = []
	for i, read in enumerate(bam):

		if read.is_unmapped:
			continue
		
		if i >= max_reads:
			break
		
		lengths.append(read.query_length)
		
	return pl.Series(lengths).mode()[0]

def natural_sort_variants(df: pl.DataFrame) -> pl.DataFrame:
	return (
		df.with_columns(
			# Create a temporary column 'chr_rank' for sorting
			pl.col("Chr")
			.str.replace("chr", "") # Remove 'chr' prefix
			.str.replace("X", "23") # Handle Sex chromosomes
			.str.replace("Y", "24")
			.str.replace("M", "25") # Handle Mitochondria if present
			.str.replace("MT", "26") # Handle Mitochondria if present
			.cast(pl.Int32, strict=False) # Convert to Integer (strict=False turns unknown contigs to null)
			.fill_null(999) # Put weird contigs at the end
			.alias("chr_rank")
		)
		.sort(["chr_rank", "Pos"]) # Sort by Rank, then Position
		.drop("chr_rank") # Remove the temp column
	)


def annotate_simple_repeats(mut_info: pl.DataFrame, simple_repeats_df: pl.DataFrame) -> pl.DataFrame:
	"""
	Annotates variants with 'Y' or 'N' if they fall within a simple repeat region.
	Assumes both DataFrames contain Chromosome and Position columns.
	"""
	
	# 1. Sort both DataFrames
	# join_asof requires the 'on' column (Pos/start) to be sorted within the 'by' groups (Chr).
	# Sorting by ["Chr", "Pos"] achieves this and keeps the data human-readable.
	variants_sorted = natural_sort_variants(mut_info) #.with_columns(pl.col("Pos").set_sorted())
	repeats_sorted = simple_repeats_df.sort(["chrom", "start"]) #.with_columns(pl.col("start").set_sorted())

	# 2. Perform the ASOF Join
	# We find the repeat segment that started most recently before (or exactly at) the variant position.
	joined = variants_sorted.join_asof(
		repeats_sorted,
		left_on="Pos",
		right_on="start",
		by_left="Chr",
		by_right="chrom",
		strategy="backward" # Look for the closest start position <= variant position
	)

	# 3. Check overlap
	# We found the closest start. Now, is the variant actually INSIDE that segment?
	# i.e., is Variant Pos <= Repeat End?
	final_df = joined.with_columns(
		pl.when(pl.col("Pos") <= pl.col("end"))
		.then(pl.lit("Y"))
		.otherwise(pl.lit("N"))
		.alias("SimpleRepeat_TRF")
	).drop(["start", "end"]) # Clean up BED columns

	return final_df


def make_mut_info(vcf_path: str, genome: pysam.FastaFile, sample_name: str, simple_repeats_bed: pl.DataFrame | None = None, ct_only: bool = False, snv_only: bool = False) -> pl.DataFrame:
	
	# Read Variants
	mut_info = read_variants(vcf_path)

	# Keep only SNVs if requested
	if snv_only:
		mut_info = snv_filter(mut_info)
  
	# Keep only C>T and G>A mutations if requested
	if ct_only:
		mut_info = ct_filter(mut_info)

	# Rename columns to match MicroSEC expected input
	mut_info = mut_info.rename({
		"chrom": "Chr",
		"pos": "Pos",
		"ref": "Ref",
		"alt": "Alt"
	})
 
	# Add Sample Name
	mut_info = mut_info.with_columns(
		pl.lit(sample_name).alias("Sample")
	)

	# MicroSEC requires variants to be at least 200 nucleotides away from the terminal ends of the chromosome/contig where they are located
	chrom_lengths = pl.DataFrame({
		"Chr": genome.references,
		"chr_len_ref": genome.lengths
	})

	# 2. Join and Filter
	mut_info = (
		mut_info.join(chrom_lengths, on="Chr", how="inner")
		.filter(
			(pl.col("Pos") > 200) & 
			(pl.col("Pos") < (pl.col("chr_len_ref") - 200))
		)
		.drop("chr_len_ref")
	)

	# Create a struct of the columns needed for the neighborhood sequence calculation
	# and then apply the function.
	mut_info = mut_info.with_columns(
		pl.struct(["Chr", "Pos", "Ref", "Alt"])
		.map_elements(lambda x: get_neighborhood_sequence(x, genome), return_dtype=pl.Utf8)
		.alias("Neighborhood_sequence")
	)
	
	# Create Mutation Type column as this is a requirement for MicroSEC input
	snv_mask = mut_info["Ref"].str.len_chars() == mut_info["Alt"].str.len_chars()
	deletion_mask = (mut_info["Ref"].str.len_chars() > 1) & (mut_info["Alt"].str.len_chars() == 1)
	insertion_mask = (mut_info["Ref"].str.len_chars() == 1) & (mut_info["Alt"].str.len_chars() > 1)
	
	mut_info = (
		mut_info.with_columns([
			pl.when(snv_mask)
			.then(mut_info["Ref"].str.len_chars().cast(pl.Utf8) + "-snv")
			.when(deletion_mask)
			.then((mut_info["Ref"].str.len_chars() - 1).cast(pl.Utf8) + "-del")
			.when(insertion_mask)
			.then((mut_info["Alt"].str.len_chars() - 1).cast(pl.Utf8) + "-ins")
			.otherwise(pl.lit("-"))
			.alias("Mut_type")
		])
	)


	# Add the simpleRepeat_TRF column by checking if the variant is in a simple repeat region
	# simple repeats bed obtained from ucsc table browser is used
	if simple_repeats_bed is not None:
		mut_info = annotate_simple_repeats(mut_info, simple_repeats_bed)
	else:
		# If no simple repeats bed provided, fill with "-"
		mut_info = mut_info.with_columns(
			pl.lit("-").alias("SimpleRepeat_TRF")
		)


	# Reoder columns accodign to MicroSEC documentation
	mut_info = mut_info.select([
		"Sample",
		"Mut_type",
		"Chr",
		"Pos",
		"Ref",
		"Alt",
		"SimpleRepeat_TRF",
		"Neighborhood_sequence"
	])

	return mut_info

def make_sample_info(
		sample_name: str,
		mut_info_path: str,
		bam_path: str,
		ref_path: str,
		simple_repeat_path: str,
		genome_build: str = "hg38",
		adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
		adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
	) -> pl.DataFrame:

	read_length = get_read_length(bam_path)

	info = {
		"sample_name": sample_name,
		"mutation_file": mut_info_path,
		"bam_file": bam_path,
		"read_length": read_length,
		"adapter_1": adapter_1,
		"adapter_2": adapter_2, # Can be None/NA
		"genome_build": genome_build, # e.g. "hg38" or "Human"
		"panel_name": "TOP",
		"ref_genome": ref_path, # Path to fasta if needed for CRAM
		"simple_repeat_bed": simple_repeat_path # Optional
	}

	return pl.DataFrame(info)

# %% [markdown]
# ### Main Execution
# #### Prepare Annotations

ffpe_vcf_paths = [os.path.abspath(path) for path in glob.glob("../vcf/EGAD00001004066/somatic_vcf-dp10/*/*.vcf.gz") if "Frozen" not in path]

sample_paths = pl.DataFrame({
    "sample_name" : [path.split("/")[-2] for path in ffpe_vcf_paths],
    "vcf_path" : ffpe_vcf_paths,
    "bam_path" : [return_path_if_exists(f"../data/EGAD00001004066/bam/{path.split("/")[-2]}/{path.split("/")[-2]}.bam", abs=True) for path in ffpe_vcf_paths],
})

# %% [markdown]
# #### Resource Loading and Directory Setup

outdir_root = "."
## MicroSEC batch script dir
batch_script_dir = f"{outdir_root}/script_microsec"
os.makedirs(batch_script_dir, exist_ok=True)

simple_repeat_path = return_path_if_exists("../data/ref/ucsc_simple-repeat_hg38.bed", abs=True)
ref_path = return_path_if_exists("../data/ref/Homo_sapiens_assembly38.fasta", abs=True)
microsec_script_path = return_path_if_exists(f"../common-ffpe-snvf/R/microsec.R", abs=True)

## Wrapper function to prepare inputs and execution script for a dataset
def process_dataset(
		dataset,
		variant_set,
		sample_paths: pl.DataFrame, 
		vcf_path_col: str = "vcf_path", 
		bam_path_col: str = "bam_path",
		sample_name_col: str = "sample_name",
		simple_repeat_path: str = simple_repeat_path,
		microsec_script_path: str = microsec_script_path,
		ref_path: str = ref_path,
		batch_script_dir: str = batch_script_dir,
		ct_only: bool = False,
		snv_only: bool = False
	):

	microsec_root = f"{outdir_root}/{dataset}/{variant_set}/microsec"

	for i in range(sample_paths.height):

		## sample specific variable setup
		sample_name = sample_paths[i, sample_name_col]
		vcf_path = sample_paths[i, vcf_path_col]
		bam_path = sample_paths[i, bam_path_col]

		print(f"{i+1}. Creating MicroSEC input for {sample_name} from VCF: {vcf_path}")

		###--------------------------------------------

		## Read reference genome
		genome = pysam.FastaFile(ref_path)

		# Read in the simple repeats bed file to annotate variants later
		simple_repeats_bed = pl.read_csv(
			simple_repeat_path, 
			separator="\t", 
			has_header=False, 
			new_columns=["chrom", "start", "end", "type", "score"],
		).select(["chrom", "start", "end"])
		
		###--------------------------------------------

		## Prepare mut_info for the sample
		print("\tPreparing mutation info")
		mut_info = make_mut_info(
			vcf_path=vcf_path, 
			genome=genome, 
			sample_name=sample_name, 
			simple_repeats_bed=simple_repeats_bed,
			ct_only=ct_only, 
			snv_only=snv_only
		)
	
		## Write mut_info
		output_dir = f"{microsec_root}/{sample_name}/inputs"
		os.makedirs(output_dir, exist_ok=True)
		mut_info_path = os.path.abspath(f"{output_dir}/{sample_name}.microsec.mut-info.tsv")
		mut_info.write_csv(mut_info_path, separator="\t")

		###----------------------------------------------
		
		## Prepare sample_info for the sample
		print("\tPreparing sample info")
		sample_info = make_sample_info(sample_name, mut_info_path, bam_path, ref_path, simple_repeat_path)

		## Write sample_info
		sample_info_path = f"{microsec_root}/{sample_name}/inputs/{sample_name}.microsec.sample_info.tsv"
		pl.DataFrame(sample_info).write_csv(sample_info_path, separator="\t")

		###-----------------------------------------------

		# Write batch execution scripts
		print("\tPreparing execution script")
		script_path = f"{batch_script_dir}/{sample_name}.microsec.sh"
		microsec_outdir = os.path.abspath(f"{microsec_root}/{sample_name}")
		
		with open(script_path, "w") as file:
			file.writelines("#!/usr/bin/env bash\n")
			file.writelines(f"Rscript {microsec_script_path} --sample_info '{sample_info_path}' --output_dir '{microsec_outdir}'\n")



process_dataset("EGAD00001004066", "somatic_vcf-dp10", sample_paths)


