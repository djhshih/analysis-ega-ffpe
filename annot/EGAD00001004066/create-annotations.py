import polars as pl

sample_metadata = pl.read_csv("samples.csv")
sample_metadata

### Create sample annotation file

sample_metadata_refined = (
	sample_metadata
 	.select([
      	'accession_id',
		'alias',
		'title',
		'biological_sex',
		'phenotype',
	])
	.rename({
    	"accession_id" : "sample_accession_id", 
     	"alias" : "sample_alias"
    })
	.with_columns(
		pl.col("title").str.split(" ").alias("title_split"),
		pl.col("title").str.replace_all(" ", "-")
	)
	.with_columns(
		pl.col("title_split").list.get(0).alias("preservation"),
        pl.col("title_split").list.get(1).alias("tissue_type"),
        pl.col("title_split").list.get(2).alias("sample_type"),
		(pl.col("title") + pl.lit("_") + pl.col("sample_alias")).alias("sample_name")
	)
	.drop("title_split")
	.select(['sample_name','sample_alias','sample_accession_id','preservation','tissue_type','sample_type','phenotype','biological_sex'])
)

sample_metadata_refined.write_csv("sample_annotations.tsv", separator="\t")

### Create fastq annotation file

fastq_metadata = pl.read_csv("sample_file.csv")

fastq_metadata = (
	fastq_metadata
	.rename({
		"file_name" : "fastq_file_name",
		"file_accession_id" : "fastq_accession_id",
	})
)

fastq_annotation = (
	fastq_metadata
 	.join(sample_metadata_refined, on=["sample_alias","sample_accession_id"], how="inner")
	.select(['title','fastq_file_name','fastq_accession_id','sample_alias','sample_accession_id','preservation','tissue_type','sample_type','phenotype','biological_sex'])
)

fastq_annotation.write_csv("fastq_annotations.tsv", separator="\t")

### Create an annotation file for both samples and fastqs

sample_fastq_annotation = (
	sample_metadata_refined
 	.join(fastq_metadata, on=["sample_accession_id","sample_alias"], how="inner")
	.group_by(['title','sample_alias','sample_accession_id','preservation','tissue_type','sample_type','phenotype','biological_sex'])
	.agg(
		pl.col("fastq_file_name").implode().list.join(";"),
		pl.col("fastq_accession_id").implode().list.join(";"),
	)
)

sample_fastq_annotation.write_csv("sample_fastq_annotations.tsv", separator="\t")

### Generate a count of samples per tissue type

sample_count = (
	sample_metadata_refined
	.group_by(["tissue_type"])
	.agg([
		((pl.col("preservation") == "FFPE") & (pl.col("sample_type") == "Tumoral")).sum().alias("n_tumor_ffpe"),
		((pl.col("preservation") == "Frozen") & (pl.col("sample_type") == "Normal")).sum().alias("n_tumor_frozen"),
		((pl.col("preservation") == "FFPE") & (pl.col("sample_type") == "Normal")).sum().alias("n_normal_ffpe"),
		((pl.col("preservation") == "Frozen") & (pl.col("sample_type") == "Tumoral")).sum().alias("n_normal_frozen"),
	])
)

sample_count.write_csv("sample_count.tsv", separator="\t")
