# TIGA Workflow

Steps for updating the TIGA dataset from sources.

## Dependencies

* R 4.2+; readr, data.table, igraph, muStat, RMySQL (Webapp: shiny, DT, shinyBS, shinysky, plotly)
* Python 3.9+; pandas, venv, [BioClients](https://github.com/jeremyjyang/BioClients)
* Java 8+; Jena, [IU\_IDSL\_JENA](https://github.com/IUIDSL/iu_idsl_jena)

## Steps

1. Download latest files from the [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/downloads). See [FTP site](ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases) for latest and all releases. Required files:
    * gwas-catalog-studies\_ontology-annotated.tsv
    * gwas-catalog-associations\_ontology-annotated.tsv
1. Download from [Experimental Factor Ontology (EFO)](https://www.ebi.ac.uk/efo/):
    * efo.owl
1. Edit LATEST\_RELEASE\_GWC.txt and LATEST\_RELEASE\_EFO.txt accordingly.
1. Create `venv` virtual environment for Python.
  1. `mkdir venv`
  1. `cd venv`
  1. `venv -i ../venv\_requirements.txt`
1. RUN [Go\_TIGA\_Workflow.sh](sh/Go_TIGA_Workflow.sh). Commands can also be run
manually as described here.
  1. Clean studies:
    * [gwascat\_gwas.R](R/gwascat_gwas.R)
  1. Clean, separate OR\_or\_beta into oddsratio, beta columns:
    * [gwascat\_assn.R](R/gwascat_assn.R)
  1. Convert EFO OWL to TSV:
    * `java -jar iu_idsl_jena-0.0.1-SNAPSHOT-jar-with-dependencies.jar`
  1. From EFO TSV create GraphML:
    * [efo\_graph.R](R/efo_graph.R)
  1. Clean traits:
    * [gwascat\_trait.R](R/gwascat_trait.R)
  1. MAPPED GENES: Separate mapped into up-/down-stream.
    * [snp2gene\_mapped.pl](perl/snp2gene_mapped.pl)
  1. Get iCite RCRs for studies via PMIDs:
    * `python3 -m BioClients.icite.Client get_stats`
  1. Get Ensembl annotations for mapped genes via EnsemblIds:
    * `python3 -m BioClients.ensembl.Client get_info`
  1. Get IDG TCRD gene annotations:
    * `python3 -m BioClients.idg.tcrd.Client listTargets`
  1. Run commands in [Go\_gwascat\_DbCreate.sh](sh/Go_gwascat_DbCreate.sh) building MySql db. Writes file `gwas_counts.tsv`.
  1. Pre-process and filter. Studies, genes and traits may be removed due to insufficient evidence, with reasons recorded.
    * [tiga\_gt\_prepfilter.R](R/tiga_gt_prepfilter.R)
  1.  Provenance for gene-trait pairs (STUDY\_ACCESSION, PUBMEDID).
    * [tiga\_gt\_provenance.R](R/tiga_gt_provenance.R)
  1. Generate variables, statistics, evidence features for gene-trait pairs.
    * [tiga\_gt\_variables.R](R/tiga_gt_variables.R)
  1. Score and rank gene-trait pairs based on selected variables.
    * [tiga\_gt\_stats.R](R/tiga_gt_stats.R)
1. TIGA web app requires files:
    1. gwascat\_gwas.tsv
    1. filtered\_genes.tsv
    1. filtered\_studies.tsv
    1. filtered\_traits.tsv
    1. gt\_provenance.tsv.gz
    1. gt\_stats.tsv.gz
    1. efo\_graph.graphml.gz
    1. gwascat\_release.txt
    1. efo\_release.txt
    1. tcrd\_info.tsv
1. TIGA download files should be copied to the [TIGA Download Directory](https://unmtid-shinyapps.net/download/TIGA/) for automated access.

## Notes

* Split comma separated fields, convert to UTF-8 characters.
* Gene-trait association variables:
  * `N_study`: studies supporting gene-trait association
  * `N_snp`: SNPs involved with gene-trait association
  * `N_snpw`(\*): SNPs involved with gene-trait association weighted by genomic distance
  * `RCRAS`(\*): RCR Aggregated Score
  * `pValue`(\*): max SNP pValues
  * `OR`: median(OR), where OR = odds ratio
  * `N_beta`: count of supporting beta values
  * `geneNtrait`: total traits associated with gene
  * `traitNgene`: total genes associated with trait
* Gene-trait scores and ranks:
  * `meanRank`: meanRank based on variables selected(\*) by benchmark validation.
  * `meanRankScore`: `100 - Percentile(meanRank)`
* MySql database currently not required for TIGA app.
