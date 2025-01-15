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

* ["TIGA: Target illumination GWAS analytics", Jeremy Yang, Dhouha Grissa,
Christophe Lambert, Cristian Bologa, Stephen Mathias, Anna Waller, David Wild,
Lars Juhl Jensen, Tudor Oprea, Bioinformatics, btab427,
https://doi.org/10.1093/bioinformatics/btab427, published 04 June 2021.](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btab427/6292081)
* [Poster presented at February 2021 Annual IDG Meeting](https://zenodo.org/record/4594201)

## Dependencies

* R 4.2+; readr, data.table, igraph, muStat, RMySQL (Webapp: shiny, DT, shinyBS, [shinysky](https://github.com/AnalytixWare/ShinySky), plotly)
* Python 3.9+; pandas, [BioClients](https://github.com/jeremyjyang/BioClients)
* Java 8+; Jena, [IU\_IDSL\_JENA](https://github.com/IUIDSL/iu_idsl_jena)

## GWAS Catalog features

* __GWAS Catalog__ studies each have a `study_accession`.
Also are associated with a publication (PubMedID), but not uniquely.
* `OR_or_BETA`: Reported odds ratio or beta-coefficient associated with strongest
SNP risk allele. Note that if an OR &lt;1 is reported this is inverted, along with
the reported allele, so that all ORs included in the Catalog are &gt;1. Appropriate
unit and increase/decrease are included for beta coefficients.
* `MAPPED_GENE`: Gene(s) mapped to the strongest SNP. If the SNP is located
within a gene, that gene is listed. If the SNP is intergenic, upstream
and downstream genes are listed. May be chromosomal
location or range (e.g. "LOC102723594 - LOC285043").
* Documentation: [methods](https://www.ebi.ac.uk/gwas/docs/methods); [curation](https://www.ebi.ac.uk/gwas/docs/methods/curation); [fileheaders](https://www.ebi.ac.uk/gwas/docs/fileheaders)
* Reference: Buniello, A. et al. (2019) The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics 2019. Nucleic Acids Res., 47, D1005–D1012.

### Issues

* Beta coefficients require units and thus are not comparable between
non-convertible units (e.g. mg vs mm). Nor are beta
coefficients comparable with OR, so it is questionable that these values
are combined in one field `OR_or_BETA`.  Current TIGA workaround is to use
simple count of beta values supporting gene-trait association.
* GWAS Catalog developers have devoted major effort to more precisely mapping traits
to EFO, and EFO has increasingly aligned with MONDO. This represents a major
improvement with regard to semantic precision and scientific rigor. However,
this also means that results from the Catalog and TIGA have changed from
release to release, which can be confusing, and presents a challenge for
aggregating studies by trait.

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

* <https://unmtid-shinyapps.net/tiga/>

## TIGA Downloads

Latest release and archives of full dataset and utility files.

* <https://unmtid-dbs.net/download/TIGA/latest/>

## Docker

TIGA may be deployed via Docker container, built with [Dockerfile.shiny](Dockerfile.shiny) and [Go_DockerBuild.sh](sh/Go_DockerBuild.sh), from [rocker/shiny](https://hub.docker.com/r/rocker/shiny).

