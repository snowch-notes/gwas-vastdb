**II. GWAS Data: Characteristics and Representation**

**2.1. Genotype Data: SNPs and Numerical Encoding**

* **Representing SNPs as Numerical Values (0, 1, 2):**
    * Explanation of how diploid organisms (like humans) have two alleles at each SNP location.
    * How these allele combinations are encoded:
        * 0: Homozygous for the reference allele (e.g., AA).
        * 1: Heterozygous (e.g., AC).
        * 2: Homozygous for the alternate allele (e.g., CC).
    * Discussion of other possible encoding schemes, and when they might be used.
* **The Genotype Matrix: Individuals as Rows, SNPs as Columns:**
    * Visual representation of the genotype matrix.
    * Explanation of how this matrix structure facilitates analysis.
    * Explanation of how missing data is represented.
    * Explanation of how data is sometimes phased, and sometimes unphased.
* **Data Compression and Storage Considerations:**
    * The repetitive nature of genotype data and its implications for compression.
    * The advantages of columnar storage for genotype data.

**2.2. Phenotype Data: Traits and Measurements**

* **Types of Phenotype Data: Binary, Continuous, Categorical:**
    * Binary: Disease status (affected/unaffected), presence/absence of a trait.
    * Continuous: Height, weight, blood pressure, biomarker levels.
    * Categorical: Disease subtypes, population groups, survey responses.
    * Examples of each type, and how they relate to GWAS studies.
* **Representing Phenotype Data in a Structured Format:**
    * How phenotype data is typically organized in tables.
    * The importance of consistent data types and formats.
    * How to handle missing phenotype data.
    * How to handle derived phenotype data.
* **Linking Phenotype Data to Genotype Data:**
    * The role of individual identifiers (IDs) in connecting genotype and phenotype data.
    * Ensuring data integrity and consistency across datasets.

**2.3. Metadata: Annotations and Demographics**

* **SNP Annotations: Location, Gene Association, Function:**
    * Chromosome and genomic position of SNPs.
    * Gene associations (e.g., which genes are located near SNPs).
    * Functional annotations (e.g., whether a SNP is in a coding region, regulatory region, or intergenic region).
    * The importance of SNP annotations for interpreting GWAS results.
* **Individual Demographics: Age, Sex, Ethnicity:**
    * The role of demographic data in controlling for confounding factors.
    * How demographic data is collected and stored.
    * How to handle sensitive demographic data, while still maintaining data integrity.
* **The Importance of Metadata for Interpretation:**
    * How metadata enhances the value of GWAS data.
    * The need for comprehensive and well-documented metadata.
    * How to handle different versions of metadata, and how to track changes to metadata.
