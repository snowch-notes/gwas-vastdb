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
 
**2.4. Exploring Real GWAS Data**

* **Introduction:**
    * This section aims to illustrate the concepts discussed in Chapter II by examining a simplified example of GWAS data.
    * It's important to note that this is a small subset for educational purposes and does not represent the complexity of a full-scale GWAS.

* **Data Source:**
    * For this example, we will use a small subset of data derived from the 1000 Genomes Project.
    * The 1000 Genomes Project provides a comprehensive catalog of human genetic variation.
    * You can access the data here: \[Provide specific link to a small, manageable subset if available, or reference the main 1000 Genomes Project site: [https://www.internationalgenome.org/](https://www.internationalgenome.org/) ]
    * For this example, we will focus on a small region of Chromosome 22.

* **Data Overview:**

    * **Genotype Data:**
        * Here's a simplified representation of what genotype data might look like (transposed format):

        |   Chromosome   |   Position   |   Individual\_ID   |   Genotype   |
        | :-----------: | :----------: | :---------------: | :---------: |
        |       22      |   16050000   |      Sample\_1      |      A\|G      |
        |       22      |   16050000   |      Sample\_2      |      G\|G      |
        |       22      |   16050000   |      Sample\_3      |      A\|A      |
        |       22      |   16050100   |      Sample\_1      |      C\|T      |
        |       22      |   16050100   |      Sample\_2      |      C\|C      |
        |       22      |   16050100   |      Sample\_3      |      T\|T      |

        * Explanation:
            * Chromosome and Position indicate the location of the SNP.
            * Individual\_ID represents a participant in the study.
            * Genotype shows the alleles for that individual at that SNP (e.g., A\|G means one allele is A, and the other is G).
    * **Phenotype Data:**
        * Here's an example of phenotype data for a binary trait (e.g., Disease X):

        |   Individual\_ID   |   Disease\_X\_Status   |
        | :---------------: | :-------------------: |
        |      Sample\_1      |          True         |
        |      Sample\_2      |          False        |
        |      Sample\_3      |          True         |

        * Explanation:
            * Individual\_ID links this data to the genotype data.
            * Disease\_X\_Status indicates whether the individual has the disease (True) or not (False).
    * **Metadata:**
        * Example of SNP metadata:

        |   Chromosome   |   Position   |   Gene   |   Annotation   |
        | :-----------: | :----------: | :-----: | :-----------: |
        |       22      |   16050000   |   BRCA2  |    Intronic    |
        |       22      |   16050100   |   BRCA2  |     Exonic     |

        * Explanation:
            * Chromosome and Position identify the SNP.
            * Gene indicates the closest known gene.
            * Annotation provides functional information about the SNP's location.

* **Illustrative Examples:**

    * **Numerical Encoding:**
        * In some cases, genotypes are numerically encoded. For example:
            * A\|A might be represented as 0
            * A\|G might be represented as 1
            * G\|G might be represented as 2
    * **Allele Frequencies:**
        * Let's take the SNP at position 16050000.
        * Allele A appears in Sample\_1 and Sample\_3 (3 times).
        * Allele G appears in Sample\_1 and Sample\_2 (2 times).
        * If these were all the samples, the frequency of allele A would be 3/5, and the frequency of allele G would be 2/5.
    * **Linking Genotype and Phenotype:**
        * To analyze if the SNP at 16050000 is associated with Disease X, we would link the genotype and phenotype data using the Individual\_ID.
        * For example, Sample\_1 has genotype A\|G and Disease\_X\_Status = True.
    * **Metadata Enrichment:**
        * The metadata tells us that the SNP at 16050000 is located within the BRCA2 gene. This information can be used to interpret the results of the GWAS.

* **Visualizations (Optional):**

    * A simple bar chart could be used to visualize the allele frequencies calculated above.

* **Transition to Database Design:**

    * The structure of this example data (chromosome, position, individual ID, genotype, phenotype, metadata) demonstrates the need for a database design that can efficiently store and query this information.
    * The next chapter will explore how to design a database to handle this type of data, focusing on the use of VAST DB.
