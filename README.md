# all code
One of the genetic risk factors of multiple sclerosis (MS) is well known as HLA class II alleles, but the genetic association with human T cell receptor (TCR) repertoire remain poorly understood. Here we investigated monozygotic twin discordant for MS with multi- validation cohorts to dissect how HLA genotype and disease status shape TCRβ composition at different aspects -- including clonal diversity, similarity, V-gene uage frequencies, CDR3 sequence composition, biochemical properties, and potential epitope recognition. We found that HLA background exerts minimal influence on overall diversity, but strongly constrains Vβ gene usage in an allele-specific manner. MS–healthy twin pairs shared the most TCRβ sequences, though shared sequences represented <5% of the total repertoire, indicating that the TCR repertoire from identical stem cells is still hypervariable despite identical genetic and environmental backgrounds. MS-associated TCRβ sequence patterns were identified in HLA-DRB1*15-positive individuals and validated in independent cohorts.
Within the TRBV19*01/TRBJ02-01*01 subset, we observed a robust and reproducible CDR3β signature characterized by aspartic acid (D) at position 6 and arginine (R) at position 7, respectively. With biochemical properties and DBSCAN clustering we also defined DRSS as a distinct motif pattern enriched in HD. GLIPH analysis independently confirmed preferential DRSS motif occurrence in healthy donors.
To explore potential antigen specification, AlphaFold structural modeling was applied to viral epitopes. A single EBV-derived peptide showed consistently higher predicted binding to DRSS-bearing TCRs, whereas no viral epitopes preferentially bound MS-associated sequences. Together, these findings identify complementary MS- and HD-associated TCRβ signatures that likely reflect distinct immunological histories and HLA-driven selection pressures.
## Basic Information
### MS twins & Validation cohorts -- propensity score matching
<p align="center">
  <img src="plots/basic_info.PNG" width=500 title="MS twins & validation cohorts">
</p>

## Clonality, generation probability, v gene usage frequency, and identical TCR sequences vs. identical HLA genes
1. TCRβ repertoire diversity and clonality were quantified using the Simpson diversity index.
2. TCRβ sequence generation probabilities (Pgene) were estimated using the OLGA (Optimized Likelihood estimate of immunoGlobulin Amino-acid sequences) algorithm implemented via the tcrdist3 framework.
3. Vβ gene usage frequencies were greately influenced by HLA genes. Distinct HLA loci exhibit characteristic V-gene usage biases.
4. Repertoire sharing was assessed by counting the number of exactly matching productive TCRβ rearrangements between each sample pair. 

## HLA-DRB1*15-MS associated TCR sequences
1. Productive TCRβ sequences observed in at least two MS-affected individuals were tested for association using Fisher’s exact test on a 2 × 2 contingency table summarizing presence/absence across MS and healthy twins. Sequences with P < 0.05 were retained as candidates.
2. Candidate sequences were then evaluated in an independent cohort of 623 HLA-typed individuals (Emerson et al., 2016). Analyses were restricted to HLA-DRB1*15–positive subjects, and the same 2 × 2 Fisher test was applied, using a stricter threshold of P < 0.001.

## Amino acid distribution
## Atchley facorts matrix, clustering and quantification
## AlphaFold modeling prediction
