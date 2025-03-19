**5. Database Design for GWAS Data**

**5.1. Relational vs. Columnar Databases: A Brief Overview**

* **Relational Databases (Row-Oriented):**
    * Data is stored in rows, where each row represents a record and each column represents an attribute.
    * Optimized for transactional workloads (OLTP), where individual records are frequently inserted, updated, or deleted.
    * Examples: PostgreSQL, MySQL, Oracle.
    * For GWAS, row-oriented storage leads to inefficient reads when analyzing specific SNPs, because large amounts of irrelevant data is read.
* **Columnar Databases (Column-Oriented):**
    * Data is stored in columns, where each column represents an attribute and values of the same attribute are stored contiguously.
    * Optimized for analytical workloads (OLAP), where large datasets are queried for analysis.
    * Examples: VAST DB, Apache Parquet, Apache Arrow.
    * For GWAS, columnar storage allows for efficient retrieval of specific SNPs, enabling faster analysis.
* **Key Differences:**
    * I/O patterns: Row-oriented reads entire records, column-oriented reads specific columns.
    * Compression: Column-oriented allows for better compression of homogeneous data.
    * Query performance: Column-oriented excels at aggregations and analytical queries.

**5.2. Why VAST DB is Well-Suited for GWAS**

* **Columnar Architecture:**
    * VAST DB's columnar storage aligns perfectly with the structure of GWAS data, where SNPs are analyzed as individual columns.
    * This architecture allows for efficient retrieval of specific SNPs, significantly speeding up analytical queries.
* **Predicate Pushdown:**
    * VAST DB's ability to push down predicates (filters) to the storage layer minimizes the amount of data transferred and processed.
    * This optimization is crucial for GWAS, where queries often involve filtering on specific SNPs or phenotype values.
* **Compression:**
    * VAST DB's compression algorithms effectively reduce the storage footprint of GWAS data, which is often highly repetitive.
    * Since genotype data has a limited number of values, great compression ratios are possible.
* **Scalability:**
    * VAST DB is designed to handle petabyte-scale datasets, which is essential for large-scale GWAS studies.
* **Arrow Data Types:**
    * VAST DB's support of Arrow data types, allows for very efficient handling of data.

**5.3. Transposed approach (Recommended)**

* **Genotype Tables: Transposed and Partitioned by Chromosome or Genomic Region:**
    * Explain the transposed design: each row represents a genotype call for an individual at a specific position.
    * Emphasize the benefits of partitioning by chromosome or genomic region:
        * Reduces table size.
        * Improves query performance.
        * Enables parallel processing.
        * Minimizes data scanned per query.
* **Phenotype Tables: Separate Tables for Different Traits:**
    * Storing phenotype data in separate tables allows for flexibility and efficient querying of specific traits.
    * This approach also makes it easier to manage and update phenotype data.
* **Metadata Tables: Storing Annotations and Demographics:**
    * Storing metadata in separate tables ensures that it is well-organized and easily accessible.
    * This approach also allows for efficient querying of metadata, such as SNP annotations or individual demographics.

**5.4. Single-Table (Wide-Column) Approach: Considerations and Trade-offs**

* **Advantages:**
    * Simpler schema design.
    * Potential for atomic updates of related data.
* **Disadvantages:**
    * Inefficient analytical queries due to row-oriented access patterns.
    * Reduced effectiveness of predicate pushdown.
    * Increased data size due to repeated row keys.
    * Complicated data type handling.
    * Partitioning complexities.
* **Why it's less optimal for GWAS:**
    * Columnar databases like VAST DB are optimized for reading entire columns, not individual cells.
    * Retrieving all values for a single SNP would require scanning the entire table.

**5.5. Data Types and Schema Design**

* **Choosing Appropriate Data Types:**
    * `STRING` for genotype data (e.g., "0/1", "1/1").
    * `INT64` for individual IDs and chromosome positions.
    * `FLOAT` or `DOUBLE` for continuous phenotype data.
    * `BOOL` for binary phenotype data.
    * `STRING` for categorical phenotype data and metadata.
* **Designing Efficient Table Schemas:**
    * Defining clear and consistent column names.
    * Using appropriate data types for each column.
    * Creating indexes or partitions to optimize query performance.

**5.6. Example Schema and Data Representation**

* **Genotype Table Example (Transposed):**
    * `Chromosome` (STRING), `Position` (INT64), `Individual_ID` (STRING), `Genotype` (STRING).
    * Explain how VCF genotype data is converted to this schema.
* **Phenotype Table Example:**
    * `Individual_ID` (STRING), `Disease_X_Status` (BOOL), `Height` (FLOAT).
* **Metadata Table Example:**
    * `SNP_ID` (STRING), `Gene` (STRING).
    * `Individual_ID` (STRING), `Age` (INT), `Sex` (STRING).
* **Data Representation:**
    * Illustrative examples of data within the defined schema.
    * How missing data is handled.
    * Example of the genotype string representing the alleles.
Absolutely. Let's add some estimated data volumes for a large genomics organization to the "Data Representation" section of Chapter 5, considering our transposed design.

**5.6. Example Schema and Data Representation**

* **Genotype Table Example (Transposed):**
    * `Chromosome` (STRING), `Position` (INT64), `Individual_ID` (STRING), `Genotype` (STRING).
    * Explain how VCF genotype data is converted to this schema.
* **Phenotype Table Example:**
    * `Individual_ID` (STRING), `Disease_X_Status` (BOOL), `Height` (FLOAT).
* **Metadata Table Example:**
    * `SNP_ID` (STRING), `Gene` (STRING).
    * `Individual_ID` (STRING), `Age` (INT), `Sex` (STRING).
* **Data Representation:**
    * Illustrative examples of data within the defined schema.
    * How missing data is handled.
    * Example of the genotype string representing the alleles.
* **Data Volume Estimates for a Large Genomics Organization:**
  * Imagine a large genomics organization with 1 million individuals genotyped at approximately 800,000 SNPs.
  * **Genotype Table:**
      * Rows: 1,000,000 individuals x 800,000 SNPs = 800 billion rows.
      * Data Size (estimated):
          * Assuming each row (Chromosome, Position, Individual_ID, Genotype) occupies approximately 50 bytes (this is a rough estimate and depends on compression).
          * Total: 800,000,000,000 rows x 50 bytes/row = 40,000,000,000,000 bytes = 40 TB (raw).
          * After compression, this number could be significantly smaller.
  * **Phenotype Table:**
      * Rows: 1,000,000 rows (one per individual).
      * Data Size (estimated):
          * Assuming 10 phenotype columns, and 100 bytes per row.
          * Total: 1,000,000 rows x 100 bytes/row = 100,000,000 bytes = 100 MB.
  * **Metadata Table:**
      * SNP Metadata: 800,000 SNPs.
      * Individual Metadata: 1,000,000 individuals.
      * Data Size (estimated):
          * Depending on the complexity of the metadata, this could range from a few GB to tens of GB.
  * **Key Considerations:**
      * These are rough estimates and actual data volumes will vary depending on factors such as the number of SNPs, the complexity of phenotypes, and the amount of metadata.
      * Compression is crucial for managing the large data volumes generated by GWAS.
      * Partitioning strategies can significantly impact query performance and data management.
      * The genotype table, due to the transposed design, is the table that takes up the vast majority of storage space.
