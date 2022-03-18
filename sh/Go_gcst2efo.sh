#!/bin/bash
###

function efoid2efoterm {
	efo_id=$1
	efo_name=$(cat /home/data/GWASCatalog/releases/2022/01/26/gwas-efo-trait-mappings.tsv |awk -F $'\t' '{print $2 "\t" $3}' |grep $efo_id |awk -F $'\t' '{print $1}' |sort -u)
	echo "$efo_name"
}
function efoid2diseasenames {
	efo_id=$1
	disease_names=$(cat /home/data/GWASCatalog/releases/2022/01/26/gwas-efo-trait-mappings.tsv |awk -F $'\t' '{print $1 "\t" $3}' |grep $efo_id |awk -F $'\t' '{print $1}' |sort -u |tr " " "_")
	echo "$disease_names"
}
# Studies for EFO_0002508 (Parkinson's disease) in GWASCatalog release 2021-09-15.
gcsts="\
GCST000002
GCST000005
GCST000261
GCST000490
GCST000528
GCST000530
GCST000567
GCST000772
GCST000855
GCST000874
GCST000950
GCST000959
GCST001126
GCST001189
GCST001205
GCST001430
GCST001445
GCST001546
GCST001617
GCST002077
GCST002353
GCST002455
GCST002544
GCST002778
GCST003059
GCST003060
GCST003600
GCST003651
GCST003652
GCST003653
GCST003922
GCST003984
GCST004902
GCST005548
GCST006153
GCST007198
GCST007780
GCST008910
GCST009324
GCST009325
GCST009373
GCST009374
GCST009512
GCST009513
GCST010049
GCST010765
GCST010990
GCST010991
GCST011037
GCST011038
GCST011039
GCST011677
GCST90000014
"
#
for gcst in $gcsts ; do
	efo_uris=$(cat /home/data/GWASCatalog/releases/2022/01/26/gwas-catalog-studies_ontology-annotated.tsv |awk -F $'\t' '{print $15 "\t" $14}' |grep $gcst |awk -F $'\t' '{print $2}' |sort -u)
	if [ "${efo_uris}" ]; then
		efo_uris=$(echo ${efo_uris} |sed 's/,/ /g')
		for efo_uri in $efo_uris ; do
			efo_id=$(basename $efo_uri)
			printf "${gcst} -> ${efo_id} EFO term: \"$(efoid2efoterm ${efo_id})\"\n"
			diseasenames=$(efoid2diseasenames ${efo_id})
			for diseasename in $diseasenames ; do
				diseasename=$(echo $diseasename |tr "_" " ")
				printf "${gcst} -> ${efo_id} disease name: \"${diseasename}\"\n"
			done
		done
	else
		printf "${gcst} -> NO EFOs\n"
	fi
done
