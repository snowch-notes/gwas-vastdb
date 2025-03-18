**V. Hands-on Demo with Sample GWAS Data**

This chapter provides a practical, step-by-step guide to working with sample GWAS data in VAST DB, allowing you to apply the concepts learned in previous chapters.

**5.1. Obtaining and Preparing Sample GWAS Data**

* **Links to Publicly Available Datasets:**
    * **1000 Genomes Project:**
        * Link: [https://www.internationalgenome.org/](https://www.internationalgenome.org/)
        * Description: This project provides a comprehensive catalog of human genetic variation, including SNP data, across diverse populations.
        * For this demo, we recommend downloading a subset of chromosome 22 VCF files to keep the data volume manageable.
        * **Download Instructions:**
            * Access the FTP server: `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/`
            * Download the following files:
                * `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
                * `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi`
            * Use `wget`, an FTP client, or your web browser to download the files. Example using `wget`:
                ```bash
                wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
                wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
                ```
            * **Subsetting the VCF:**
                * To create a smaller subset, use `bcftools`. Example:
                    ```bash
                    bcftools view -r 22:16000000-17000000 ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -Oz -o chr22_subset.vcf.gz
                    bcftools index chr22_subset.vcf.gz
                    ```
    * **Simulated GWAS Data:**
        * For a simpler and quicker start, consider using simulated GWAS data. Tools like `PLINK` or `Hail` can generate simulated genotype and phenotype data.
        * This is especially helpful if you don't have immediate access to real GWAS datasets.
* **Steps for Data Conversion and Formatting:**
    * **VCF to Apache Arrow:**
        * Use Python libraries like `pysam` and `pyarrow` to convert VCF files to Apache Arrow format, which is compatible with VAST DB.
        * Example Python code snippet (illustrative):

```python
        import pysam
        import pyarrow as pa
        import pyarrow.parquet as pq
        import pandas as pd

        def vcf_to_arrow(vcf_file, output_parquet):
            vcf = pysam.VariantFile(vcf_file)
            data = []
            for record in vcf:
                row = [record.chrom, record.pos]
                for sample in record.samples:
                    genotype = record.samples[sample]['GT']
                    if genotype is None:
                        row.append(None) #handle missing data.
                    else:
                        alleles = [str(record.alleles[i]) for i in genotype if i is not None]
                        row.append("|".join(alleles))

                data.append(row)

            columns = ["Chromosome", "Position"] + list(vcf.header.samples)
            table = pa.Table.from_pandas(pd.DataFrame(data, columns=columns))
            pq.write_table(table, output_parquet)

        # Example usage
        vcf_to_arrow("chr22_subset.vcf.gz", "chr22_subset.parquet")
```

    * **Phenotype Data:**
        * Prepare a CSV or TSV file with individual IDs and corresponding phenotype values.
        * Use `pandas` and `pyarrow` to convert the phenotype data to Apache Arrow format.
    * **Metadata:**
        * Create CSV or TSV files for SNP annotations and individual demographics.
        * Convert these files to Apache Arrow format using `pandas` and `pyarrow`.
    * Ensure that Individual IDs are consistent between the phenotype data and the genotype data.

**5.2. Setting Up VAST DB for the Demo**

* **Instructions for Creating Schemas and Tables:**
    * Provide example `vastdb` SDK code to create schemas and tables for genotype, phenotype, and metadata data.
    * Example:

```python
        import vastdb
        import pyarrow as pa

        # Connect to VAST DB
        session = vastdb.connect(
            endpoint='http://vip-pool.v123-xy.VastENG.lab',
            access=AWS_ACCESS_KEY_ID,
            secret=AWS_SECRET_ACCESS_KEY
        )

        with session.transaction() as tx:
            bucket = tx.bucket("my_gwas_bucket")
            schema = bucket.create_schema("gwas_demo")

            # Genotype Table Schema
            genotype_schema = pa.schema([
                ('Chromosome', pa.string()),
                ('Position', pa.int64()),
                # Add columns for each individual based on the VCF file
            ])
            genotype_table = schema.create_table("Genotype_chr22", genotype_schema)

            # Phenotype Table Schema
            phenotype_schema = pa.schema([
                ('Individual_ID', pa.string()),
                ('Disease_Status', pa.bool_()),
            ])
            phenotype_table = schema.create_table("Phenotype", phenotype_schema)

            # Metadata Table Schema
            metadata_schema = pa.schema([
                ('SNP_ID', pa.string()),
                ('Gene', pa.string()),
            ])
            metadata_table = schema.create_table("Metadata", metadata_schema)
```

    * Emphasize partitioning the genotype table by chromosome or genomic region.
* **Example Data Loading Scripts Using the `vastdb` SDK:**
    * Provide Python scripts to load the prepared Arrow data into the VAST DB tables.
    * Example:

```python
        import vastdb
        import pyarrow.parquet as pq

        # Connect to VAST DB
        session = vastdb.connect(
            endpoint='http://vip-pool.v123-xy.VastENG.lab',
            access=AWS_ACCESS_KEY_ID,
            secret=AWS_SECRET_ACCESS_KEY
        )

        with session.transaction() as tx:
            bucket = tx.bucket("my_gwas_bucket")
            schema = bucket.schema("gwas_demo")
            genotype_table = schema.table("Genotype_chr22")
            phenotype_table = schema.table("Phenotype")
            metadata_table = schema.table("Metadata")

            # Load genotype data
            genotype_data = pq.read_table("chr22_subset.parquet")
            genotype_table.insert(genotype_data)

            # Load phenotype data
            phenotype_data = pq.read_table("phenotype.parquet")
            phenotype_table.insert(phenotype_data)

            #Load metadata.
            metadata_data = pq.read_table("metadata.parquet")
            metadata_table.insert(metadata_data)
```

You are absolutely right. My apologies! Let's complete the content for section 5.3 and 5.4.

**5.3. Running Example Queries and Analyses (Continued)**

* **Demonstrating Basic and Advanced Queries from Chapter IV:**
    * Provide example queries that can be run on the loaded data.
    * Example:

```python
        # Select individuals with a specific genotype
        reader = genotype_table.select(predicate=(_.Position == 16050000))
        result = reader.read_all()
        print(result)

        # Select individuals with a specific phenotype
        reader = phenotype_table.select(predicate=(_.Disease_Status == True))
        result = reader.read_all()
        print(result)

        # Join genotype and phenotype tables
        # (Requires fetching data and using pandas or pyarrow for the join)
        genotype_reader = genotype_table.select(columns=['Chromosome','Position'])
        phenotype_reader = phenotype_table.select(columns=['Individual_ID','Disease_Status'])

        genotype_result = genotype_reader.read_all().to_pandas()
        phenotype_result = phenotype_reader.read_all().to_pandas()

        # Assuming Individual_ID is available in genotype_result
        joined_result = genotype_result.merge(phenotype_result, on='Individual_ID', how='inner')
        print(joined_result)
```

* **Performing Simple Statistical Analyses:**
    * Demonstrate how to calculate allele frequencies or perform simple association tests using Python libraries.
    * Example:

```python
        # Calculate allele frequencies (using pandas after fetching data)
        genotype_reader = genotype_table.select(columns=['Chromosome', 'Position'])
        genotype_results = genotype_reader.read_all().to_pandas()
        # count the number of times each unique value occurs in the position column.
        allele_counts = genotype_results['Position'].value_counts()
        print(allele_counts)
```

* **Visualizing Results:**
    * Show how to create simple visualizations using `matplotlib` or `seaborn`.
    * Example:

```python
        import matplotlib.pyplot as plt

        # Create a histogram of allele frequencies
        allele_counts.plot(kind='bar')
        plt.title("Allele Frequency Histogram")
        plt.xlabel("Position")
        plt.ylabel("Frequency")
        plt.show()
```

**5.4. Performing Comparative Timing Tests**

* **Comparing Query Times Between `vastdb` and Other Methods:**
    * Provide Python scripts that measure the execution time of queries in `vastdb` and compare them to other methods (e.g., using `pandas` on the same data loaded into memory).
    * Example:

```python
        import time
        import pandas as pd

        # Time a query in vastdb
        start_time_vastdb = time.time()
        reader = genotype_table.select(predicate=(_.Position == 16050000))
        result_vastdb = reader.read_all().to_pandas()
        vastdb_time = time.time() - start_time_vastdb

        # Time the same query using pandas (if data is loaded into memory)
        start_time_pandas = time.time()
        # load the parquet file into a pandas dataframe.
        pandas_df = pd.read_parquet('chr22_subset.parquet')
        pandas_result = pandas_df[pandas_df['Position'] == 16050000]
        pandas_time = time.time() - start_time_pandas

        print(f"vastdb time: {vastdb_time}")
        print(f"pandas time: {pandas_time}")
```

* **Providing the Scripts Used to Perform the Timing Tests:**
    * Make the Python scripts used for timing tests available for download.
    * Explain how to run the scripts and interpret the results.
    * Ensure that the tests are repeatable.
    * Ensure that the tests are fair, and that the data is loaded into memory before the pandas test.
