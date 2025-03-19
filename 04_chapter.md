**4. GWAS Data File Types and Structures**

**4.1. VCF (Variant Call Format)**

* **Overview:**
    * The Variant Call Format (VCF) is a standard text file format used to store genetic variant information.
    * It's designed to represent SNPs, indels, and other genetic variations found in sequencing data.
    * VCF is highly versatile, accommodating a wide range of variant data, including genotype calls, quality scores, and annotations.
* **Structure:**
    * A VCF file is divided into two primary sections: the header and the data body.
    * **Header Section:**
        * The header contains metadata, such as information about the reference genome, sample IDs, and annotations.
        * Lines in the header begin with "##".
        * It defines the structure of the data section, including contig definitions (chromosomes), INFO field descriptions, and FORMAT field descriptions.
    * **Data Section:**
        * The data section contains the actual variant information, with each row representing a single variant.
        * Key columns include:
            * **CHROM:** Chromosome or contig name.
            * **POS:** Position of the variant on the chromosome.
            * **ID:** Variant identifier (e.g., dbSNP rsID).
            * **REF:** Reference allele(s).
            * **ALT:** Alternate allele(s).
            * **QUAL:** Quality score of the variant call.
            * **FILTER:** Filter status (e.g., PASS, FAIL).
            * **INFO:** Additional annotations (e.g., allele frequencies, gene annotations).
            * **FORMAT:** Describes the format of the sample genotype data.
            * **Sample Columns:** Genotype calls and other sample-specific data.
* **Compressed VCF and Index Files:**
    * Due to their large size, VCF files are often compressed using gzip (VCF.gz).
    * Index files (VCF.gz.tbi) enable fast random access to compressed VCF files, crucial for efficient data retrieval.
    * Tools like `bcftools` and `tabix` are used to manage these files.
* **Example VCF Data Snippet:**

```
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1 Sample2
    22      16050000        rs123   A       G       100     PASS    AF=0.5  GT      0/1     1/1
    22      16050100        rs456   C       T       120     PASS    AF=0.25 GT      0/0     1/0
```

**4.2. Understanding VCF in the Context of GWAS**

* **VCF Usage in GWAS:**
    * VCF files are the primary format for storing and exchanging genotype data in GWAS.
    * They provide a comprehensive record of genetic variations across a population, enabling researchers to identify associations between variants and traits.
* **Genotype Data in VCF:**
    * The core data for GWAS is found in the FORMAT and sample columns.
    * **FORMAT column:**
        * Defines the fields present in the sample columns, such as GT (genotype), AD (allele depths), and GQ (genotype quality).
    * **Sample columns:**
        * Contain the actual genotype calls for each individual at each variant position.
        * Genotype calls are typically represented as "0/0" (homozygous reference), "0/1" (heterozygous), or "1/1" (homozygous alternate).
* **Extracting Genotype Information:**
    * Tools like `pysam` (Python) and `bcftools` (command-line) are used to parse VCF files and extract genotype data.
    * This extracted data is then used for statistical analyses and association testing.
    * The genotype data is often transposed from the wide format of the vcf file, to the long format needed for database storage, as shown in the example code in the demo chapter.

**4.3. Other Relevant File Considerations**

* **Phenotype and Metadata Files:**
    * CSV and TSV files are commonly used to store phenotype data (traits, measurements) and metadata (SNP annotations, demographic information).
    * These formats are simple and widely compatible with data analysis tools.
* **Other Genotype Formats:**
    * While VCF is the primary format in this demo, it's important to acknowledge other formats like PLINK (BED, BIM, FAM).
    * PLINK files are used by the PLINK software.
    * This tutorial and demo will concentrate on VCF files, due to their wide use, and the way that they are easily transposed.
* **Summary Statistics Files:**
    * Summary statistics files contain the results of GWAS analyses (SNP IDs, p-values, effect sizes).
    * These files are used for meta-analyses and fine-mapping.
