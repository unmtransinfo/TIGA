# TIGA: Target Illumination GWAS Analytics

Genome-wide association studies (GWAS) generate large volumes of associations between
genomic variations and phenotypic traits. However, assessing strength, specificity
and relevance of associations, and semantically valid aggregation of associations,
for applications such as drug target prioritization, is challenging. This project
addresses this challenge.

* [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/)
* [Experimental Factor Ontology (EFO)](https://www.ebi.ac.uk/efo/)

The __NHGRI-EBI GWAS Catalog__ is itself an expertly designed and curated aggregation of
GWAS results and metadata, and the primary data source for this project. Our effort
builds upon the __GWAS Catalog__ with more specific applications and use cases,
focused on protein-coding genes and well defined traits semantically related to disease
states relevant to discovery of drugs and druggable targets.

See related [poster](https://www.slideshare.net/jeremyjyang/gwas-explorer-drug-target-illumination-by-evidence-aggregation-and-multivariate-mu-scoring) presented at Open Targets - IDG meeting, Wellcome Genome Campus, Hinxton UK, Nov 15, 2019.

## Dependencies

* Python 3.6+; package [BioClients](https://github.com/jeremyjyang/BioClients) 
* R 3.6+; packages readr, data.table, shiny, DT, shinyBS, plotly, igraph
* [dqshiny](https://github.com/daqana/dqshiny) dev version (late 2019) via `remotes::install_github("daqana/dqshiny")` to resolve update\_autocomplete\_input bug.

## GWAS Catalog features

* __GWAS Catalog__ studies each have a `study_accession`.
Also are associated with a publication (PubMedID), but not uniquely.
* `OR_or_BETA`: Reported odds ratio or beta-coefficient associated with strongest
SNP risk allele. Note that if an OR &lt;1 is reported this is inverted, along with
the reported allele, so that all ORs included in the Catalog are &gt;1. Appropriate
unit and increase/decrease are included for beta coefficients.
* `MAPPED_GENE`: Gene(s) mapped to the strongest SNP. If the SNP is located
within a gene, that gene is listed. If the SNP is intergenic, the upstream
and downstream genes are listed, separated by a hyphen. May be chromosomal
location or range (e.g. "LOC102723594 - LOC285043").
* Documentation:
  * <https://www.ebi.ac.uk/gwas/docs/fileheaders>
  * <https://www.ebi.ac.uk/gwas/docs/methods>
  * <https://www.ebi.ac.uk/gwas/docs/methods/curation>
* Reference: Welter D, MacArthur J, Morales J, Burdett T, Hall P, Junkins H,
Klemm A, Flicek P, Manolio T, Hindorff L, and Parkinson H. The NHGRI
GWAS Catalog, a curated resource of SNP-trait associations. Nucleic
Acids Research, 2014, Vol. 42 (Database issue): D1001-D1006.

### Issues

* Beta coefficients require units and thus are not comparable between
non-convertible units (e.g. mg vs mm). Nor are beta
coefficients comparable with OR, so it is questionable that these values
are combined in one field `OR_or_BETA`.

## Features of TIGA

* Protein-coding gene to disease association focus.
* Evidence assessment based on confirmatory statistics.
* __iCite__ annotations from __iCite API__, via PMIDs from __GWAS Catalog__.
* Visualization of associations for a given disease by scatter plot of
effect size versus __meanRankScore__, inverse multivariate mean rank
of benchmark-validated variables.

## TIGA Workflow

See [WORKFLOW.md](doc/WORKFLOW.md) for details describing
how to update the TIGA dataset from sources.

## TIGA Application

* Currently at <http://unmtid-shinyapps.net/tiga/>
