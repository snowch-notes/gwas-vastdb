**IV. Querying and Analysis in VAST DB (Using the vastdb Python SDK)**

This chapter will guide you through querying and analyzing GWAS data stored in VAST DB using the `vastdb` Python SDK. We'll cover basic and advanced queries, performance optimization, and snapshot strategies, all within the context of the SDK's capabilities.

**4.1. Basic Queries: Selecting Individuals and SNPs**

* **Connecting to VAST DB and Accessing Tables:**
    * First, you'll need to establish a connection to your VAST DB cluster using the `vastdb.connect()` function.
    * Example:

```python
import pyarrow as pa
import vastdb
from vastdb import _ #used for predicate pushdown

AWS_ACCESS_KEY_ID = "YOUR_AWS_ACCESS_KEY_ID"
AWS_SECRET_ACCESS_KEY = "YOUR_AWS_SECRET_ACCESS_KEY"

session = vastdb.connect(
    endpoint='http://vip-pool.v123-xy.VastENG.lab',
    access=AWS_ACCESS_KEY_ID,
    secret=AWS_SECRET_ACCESS_KEY
)

with session.transaction() as tx:
    bucket = tx.bucket("my_gwas_bucket")
    schema = bucket.schema("gwas_schema")
    genotype_table = schema.table("Genotype_Chromosome_1")
    phenotype_table = schema.table("Phenotype_Disease_X")
```

* **Using `WHERE` Clauses to Filter Data (Predicate Pushdown):**
    * The `vastdb` SDK allows you to apply filters (predicates) directly to your queries.
    * Example: "Retrieve all individuals with `SNP_1 = 2`."

```python
    reader = genotype_table.select(predicate=(_.SNP_1 == 2))
    result = reader.read_all()
    print(result)
```

    * Example: "Retrieve all individuals with `Disease_X_Status = TRUE`."

```python
    reader = phenotype_table.select(predicate=(_.Disease_X_Status == True))
    result = reader.read_all()
    print(result)
```

    * Example: "Retrieve all individuals with `SNP_1 = 1` and `SNP_2 = 0`."

```python
    reader = genotype_table.select(predicate=(_.SNP_1 == 1) & (_.SNP_2 == 0))
    result = reader.read_all()
    print(result)
```

    * Example: "Retrieve all individuals with `Height > 180`."

```python
    reader = phenotype_table.select(predicate=(_.Height > 180))
    result = reader.read_all()
    print(result)
```

    * Example: "Retrieve all individuals with `SNP_5 within the isin list [1,2]`"

```python
    reader = genotype_table.select(predicate=(_.SNP_5.isin([1,2])))
    result = reader.read_all()
    print(result)
```

* **Retrieving Specific SNP Values or Phenotype Data (Projections):**
    * You can select specific columns using the `columns` parameter in `table.select()`.
    * Example: "Retrieve the `SNP_3` value for individual `1001`."

```python
    reader = genotype_table.select(
        columns=['SNP_3'],
        predicate=(_.Individual_ID == 1001)
    )
    result = reader.read_all()
    print(result)
```

    * Example: "Retrieve the `Height` and `Weight` for all individuals."

```python
    reader = phenotype_table.select(columns=['Height', 'Weight'])
    result = reader.read_all()
    print(result)
```

**4.2. Advanced Queries: Statistical Analysis**

* **Aggregations and Calculations Across SNPs:**
    * While the `vastdb` SDK primarily focuses on data retrieval, you can perform aggregations using `pyarrow` after fetching the data.
    * Example: "Calculate the average `SNP_1` value."

```python
    reader = genotype_table.select(columns=['SNP_1'])
    result = reader.read_all()
    snp_1_array = result['SNP_1'].to_numpy()
    average_snp_1 = snp_1_array.mean()
    print(average_snp_1)
```

* **Joining Genotype and Phenotype Tables:**
    * Joining tables requires fetching the data and then using libraries like `pandas` or `pyarrow` for the join operation.
    * Example: "Retrieve the `SNP_1` value and `Disease_X_Status` for all individuals."

```python
    genotype_reader = genotype_table.select(columns=['Individual_ID', 'SNP_1'])
    phenotype_reader = phenotype_table.select(columns=['Individual_ID', 'Disease_X_Status'])

    genotype_result = genotype_reader.read_all().to_pandas()
    phenotype_result = phenotype_reader.read_all().to_pandas()

    joined_result = genotype_result.merge(phenotype_result, on='Individual_ID')
    print(joined_result)
```

* **Performing basic GWAS calculations within VAST DB:**
    * Calculations such as allelic frequencies and linear regressions, are normally performed using statistical packages after the data has been retrieved.

**4.3. Performance Optimization: Predicate Pushdown and Partitioning**

* **How VAST DB Optimizes Query Execution:**
    * The `vastdb` SDK leverages VAST DB's predicate pushdown capabilities, as demonstrated in the query examples.
    * Partitioning strategies, implemented during schema design, further enhance performance.
    * Using the columnar selection of the SDK, only the needed columns are read.
* **Best Practices for Query Optimization:**
    * Use predicates to filter data as early as possible.
    * Select only the necessary columns.
    * Ensure your tables are properly partitioned.

**4.4. Snapshot Strategies**

* **Why and When to Take Snapshots:**
    * While the `vastdb` SDK doesn't directly manage snapshots, VAST DB itself supports snapshots.
    * Snapshots are essential for data backup and recovery, especially after major data ingestion or processing steps.
    * Consult the VAST DB documentation for snapshot management procedures.
* **Considerations for Large Datasets:**
    * Because of snapshot limitations, replication of data to another vast cluster should be considered.
    * For extremely large datasets, prioritize snapshots of critical subsets.
    * Implement a robust snapshot retention policy.
