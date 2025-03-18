**II. GWAS Data: Characteristics and Representation**

**2.1. Genotype Data: SNPs and Numerical Encoding**

* **Representing SNPs as Numerical Values (0, 1, 2):**
    * Humans are diploid, meaning they inherit two copies of each chromosome, one from each parent. Therefore, at any given SNP location, an individual will have two alleles.
    * To represent these allele combinations numerically:
        * **0:** Represents a homozygous genotype for the reference allele (e.g., both alleles are "A," denoted as AA).
        * **1:** Represents a heterozygous genotype (e.g., one allele is "A," and the other is "C," denoted as AC).
        * **2:** Represents a homozygous genotype for the alternate allele (e.g., both alleles are "C," denoted as CC).
    * Other encoding schemes exist, such as using "A," "C," "G," and "T" directly, or using binary encoding. The 0,1,2 encoding however is very common in GWAS.
* **The Genotype Matrix: Individuals as Rows, SNPs as Columns:**
    * The genotype data is typically organized into a matrix where:
        * Rows represent individual participants in the GWAS.
        * Columns represent the SNPs being analyzed.
    * Each cell in the matrix contains the numerical genotype value (0, 1, or 2) for a specific individual at a specific SNP.
    * This matrix structure is ideal for efficient statistical analysis, as it allows for vectorized operations and calculations across large datasets.
    * Missing data is very common in GWAS. Missing data is often represented as a -1, or as a null value.
    * Phased data refers to data where the alleles inherited from each parent are known, and unphased data refers to data where the allelic combination is known, but not the parental origin.
* **Data Compression and Storage Considerations:**
    * Genotype data often contains repetitive patterns, as many individuals share the same genotypes at certain SNPs.
    * Columnar storage is highly advantageous because it stores values for each SNP (column) contiguously, enabling efficient compression algorithms to reduce storage space.
    * Since the data is mostly composed of 0, 1, and 2, very good compression ratios can be achieved.

**2.2. Phenotype Data: Traits and Measurements**

* **Types of Phenotype Data: Binary, Continuous, Categorical:**
    * **Binary:**
        * Represents two distinct categories, such as disease presence (1) or absence (0).
        * Example: A study investigating the genetic basis of type 2 diabetes might use a binary phenotype indicating whether individuals have the disease.
    * **Continuous:**
        * Represents measurements on a continuous scale.
        * Examples: Height (measured in centimeters), blood pressure (measured in mmHg), or levels of a specific biomarker in the blood.
    * **Categorical:**
        * Represents distinct categories with no inherent order.
        * Examples: Blood type (A, B, AB, O), or different subtypes of a disease.
* **Representing Phenotype Data in a Structured Format:**
    * Phenotype data is typically stored in tables, with each row representing an individual and each column representing a specific trait or measurement.
    * Consistent data types (e.g., integers, floating-point numbers, strings) are essential for accurate analysis.
    * Missing phenotype data is also very common. It is important to have a standard way of representing missing phenotype data.
    * Derived phenotype data is data that is calculated from other phenotype data. For example, BMI is derived from height and weight.
* **Linking Phenotype Data to Genotype Data:**
    * A unique individual identifier (ID) is crucial for linking phenotype data to genotype data.
    * This ID ensures that the correct phenotype information is associated with the corresponding genotype information for each individual.
    * Data integrity is very important, because if the phenotype and genotype data are not correctly linked, the results of the GWAS will be invalid.

**2.3. Metadata: Annotations and Demographics**

* **SNP Annotations: Location, Gene Association, Function:**
    * SNP annotations provide essential contextual information about each SNP.
    * Chromosome and genomic position: Where the SNP is located within the genome.
    * Gene associations: Which genes are located near the SNP.
    * Functional annotations: The potential functional impact of the SNP (e.g., whether it is located in a coding region, regulatory region, or intergenic region).
    * These annotations are crucial for interpreting GWAS results and understanding the biological mechanisms underlying genetic associations.
* **Individual Demographics: Age, Sex, Ethnicity:**
    * Demographic data plays a critical role in controlling for confounding factors that can influence GWAS results.
    * Age, sex, and ethnicity can all affect the prevalence of certain traits or diseases.
    * This data is collected through questionnaires, medical records, or other sources.
    * Because demographic data can be sensitive, it is important to handle it with care, and to ensure that it is protected.
* **The Importance of Metadata for Interpretation:**
    * Metadata enhances the value of GWAS data by providing context and enabling more accurate interpretation of results.
    * Comprehensive and well-documented metadata is essential for reproducibility and data sharing.
    * It is important to track changes to metadata, and to have a way of handling different versions of metadata.
