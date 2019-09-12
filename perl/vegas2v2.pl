#!/usr/bin/perl
use Cwd;
use File::Spec::Functions qw(rel2abs); 
use File::Basename; 
use Scalar::Util qw(looks_like_number);
use File::Copy;
use List::Util 'shuffle';
use List::Util qw(sum);

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# VEGAS
$version = "2.01.17"; # VEGAS2 version 2 and Month and Year of last major modification

###--- BEGIN FRONT MATTER ---###

print "\n##############################################################\n";
print "#                                                            #\n";
print "# VEGAS2: Versatile Gene-based Association Study 2 $version #\n";
print "#                                                            #\n";
print "##############################################################\n\n";

print "This version was developed by Dr. Aniket Mishra and Asso. Prof. Stuart Macgregor.\n";
print 'If you have any issues/queries or suggestions please send an email to aniket.mishra@qimrbeghofer.edu.au';

my $USAGE = <<'END_USAGE';
Usage: perl vegas2v2 -G [-snpandp] <snpandp.file> [-custom] <plinkbinaryfile> [-glist] <genelocationfile> [-upper] 0 [-lower] 0 [-chr] 1..23 [-out] <output.name> [-genelist] <genelist.file> 
          [-max] 10000000 [-bestsnp/-top 1..100]
       perl vegas2v2 -P [-geneandp] <geneandp.file> [-geneandpath] <genepathway.annot> [-maxsample] 1000000
Options:
    -G/-P	     : The first argument should be -G or -P to perform gene-based or pathway-based analysis respectively.
    -snpandp     : Immediately after -G argument and space provide -snpandp option and the name of text file containing snp id and p value.
    -custom      : Provide the plink binary format genotype files to compute pairwise LD between variants
    -glist       : Provide human genome annotated gene location file. You can download it from NCBI, UCSC databases. We also provide hg19 gene annotation file on VEGAS2 webpage.
    -upper       : Provide the 5' UTR end gene-boundary.
    -lower       : Provide the 3' UTR end gene-boundary.
    -chr         : Provide the chromosome number to perform analysis on, by default it work on all autosomes.
    -out         : Provide outfile name.
    -genelist    : Provide a file with list of gene Symbols to perform test on.
    -max         : Provide the maximum number of simulation to perform. Must be above 1e6.
    -bestsnp     : To perform best SNP test.
    -top         : To perform top percentage test. Imediately after space provide a number between 1 to 100 to perfom top percentage test.
    -geneandp    : Provide the text file containing gene id and p-value respectively.
    -geneandpath : Provide gene-pathway annotation file.
    -maxsample   : Provide the maximum number of resamples to perform. Must be above 1e5.
END_USAGE



$dir = dirname(rel2abs($0));
# By defualt VEGAS performs -top 100 test
$percentage = "1";
$outfile = "gene-basedoutput.out";

# Options
for($i = 1; $i <= $#ARGV; ++$i){
	$field = $ARGV[$i];
	if(length($field)>1){
		if($field eq "-snpandp"){ # read input for gene-based analysis
			$snpandp = "$ARGV[$i+1]";
		}
#For pathway analysis
		if($field eq "-geneandp"){
			$geneandp = "$ARGV[$i+1]";
		}
		if($field eq "-geneandpath"){
			$geneandpath = "$ARGV[$i+1]";
		}
		if($field eq "-maxsample"){ # Specify max number of sims (default is 1000000)
			$maxsample = "$ARGV[$i+1]";
			if($maxsample < 1000000){
				print "$USAGE\n";
				die "Error: maximum number of -maxsample must be greater than 100000\n";
			}
		}
#For pathway analysis finished
#Options for gene-based analysis
		if($field eq "-top"){ # Top X test
			$topten = "$field";
			$perc = $ARGV[$i+1];
			if($perc > 100 || $perc < 0){
				print "$USAGE\n";
				die "Error: Percentage must be between 0 and 100\n";
			}
			$percentage = $perc/100;
		}
		if($field eq "-topsnp"){ #Do TopSNP test
			$topsnp = 1;
		}
		if($field eq "-chr"){ # Do test for specific chromosome
			$dochr = "$ARGV[$i+1]";
			if ($dochr > 24 || $dochr < 1){
				print "$USAGE\n";
				die "Error: -chr must be a number between 1 and 22\n";
			}
		}
		if($field eq "-out"){ # Specify outfile
			$outfile = "$ARGV[$i+1]";
		}
		if($field eq "-genelist"){ # Do test for list of genes
			$list = "$ARGV[$i+1]";
			if(!-e "$list"){
				print "$USAGE\n";
				die "Error: $list not found\n";
			}
			open(LIST, "$list");
			@genelist = <LIST>;
			close LIST;
		}
		if($field eq "-max"){ # Specify max number of sims (default is 1000000)
			$max = "$ARGV[$i+1]";
			if($max < 1000000){
				print "$USAGE\n";
				die "Error: maximum number of sims must be greater than 1000000\n";
			}
		}
		if($field eq "-keeptimestamp"){ # Do not delete time_stamp folder
			$keeptimestamp = 1;
		}
		if($field eq "-upper"){ # Upper boundary (use with -custom)
			$upper = "$ARGV[$i+1]";
		}
		if($field eq "-lower"){ # Lower boundary (use with -custom)
			$lower = "$ARGV[$i+1]";
		}
		if($field eq "-custom"){ # Use custom individual genotypes - do not confuse with previous -custom (replaced with -pop)
			$custom = "$ARGV[$i+1]";
		}
		if($field eq "-glist"){ #gene location file hg19/hg18 etc
			$glist = "$ARGV[$i+1]";
		}
		# For now just keep it here
		if($field eq "-prune"){ #gene location file hg19/hg18 etc
			$prunersquare = $ARGV[$i+1];
			if($prunersquare > 0.9999 || $prunersquare < 0.001){
				die "Error: Prune r2 parameter must be more than 0.001 and less than 0.9999\n";
			}
		}
	}
}

# Make working directory
my ($Second, $Minute, $Hour, $Day, $Month, $Year, $WeekDay, $DayOfYear, $IsDST) = localtime(time) ;
$time_stamp = $DayOfYear.$Hour.$Minute.$Second.$outfile.$dochr;
mkdir "$time_stamp", 0777 or die "Error: cannot make working directory. Check that vegas has correct permission settings: $!";
print "\n\nCreating working directory: $time_stamp\n";


# If -G check all required arguments are available
if ( $ARGV[0] eq "-G" ) {
	print "Performing gene-based analysis\n";
	print "Options in effect: \n";
	
	# Check if snpandp file exist and copy to working directory
	if( !-e "$snpandp" ) {
		print "$USAGE\n";
		die "Error: snpandp file not found\n";
	}
#	copy( "$snpandp", "$time_stamp/snpandp" ) or die "Copy failed: $!";
	system("sort -k1,1 $snpandp >$time_stamp/snpandp");

	print "\t-snpandp $snpandp\n";
	
	# Check if glist file exist and copy to working directory
	if( !-e "$glist" ) {
		print "$USAGE\n";
		die "Error: glist file not found. Please provide gene annotation file\n";
	}
	copy( "$glist", "$time_stamp/tempglist" ) or die "Copy failed: $!";
	print "\t-glist $glist\n";
	
	# 3 Check if genotype file exist
	if( defined( $custom ) ) {
		if( !-e "$custom.bed" ) {
			print "$USAGE\n";
			die "Error: $custom.bed not found\n";
		}
		if( !-e "$custom.bim" ) {
			print "$USAGE\n";
			die "Error: $custom.bim not found\n";
		}
		if( !-e "$custom.fam" ) {
			print "$USAGE\n";
			die "Error: $custom.fam not found\n";
		}
#		copy("$custom.bed","$time_stamp/custom.bed") or die "Copy failed: $!";
#		copy("$custom.bim","$time_stamp/custom.bim") or die "Copy failed: $!";
#		copy("$custom.fam","$time_stamp/custom.fam") or die "Copy failed: $!";
		}
	else {
		print "$USAGE\n";
		die "Error: please provide reference genotype files to perform gene-based test\n";
	}
	print "\t-custom $custom\n";
	if( defined( $lower ) ) {
	print "\t-lower $lower\n";
	}
	else {
	$lower = "0";
	print "\t-lower $lower\n";
	}
	if( defined( $upper ) ) {
	print "\t-lower $upper\n";
	}
	else {
	$upper = "0";
	print "\t-upper $upper\n";
	}
	if( defined($prunersquare) ) {
		print "\t-prune $prunersquare\n";
	}
	if( defined($topsnp) ) {
		print "\t-topsnp\n";
	}
	else {
		$perc = $percentage * 100;
		print "\t-top $perc\n";
	}
	if( defined($dochr) ) {
		print "\t-chr $dochr\n";
	}
	if( defined($list) ) {
		print "\t-genelist $list\n";
		if(!-e "$list"){
			print "$USAGE\n";
			die "Error: genelist file not found\n";
		}
		copy("$list","$time_stamp/customgenelist") or die "Copy failed: $!";
	}
	if( defined($max)) {
		print "\t-max $max\n";
	}
	if( defined($keeptimestamp) ) {
		print "\t-keeptimestamp\n";
	}
	if( defined($outfile) ) {
		print "\t-out $outfile";
	}
	
	# Other validations
	if(defined($topten) and defined($topsnp)){
		print "$USAGE\n";
		die "Error: Cannont do -topten and -topsnp test concurrently\n";
	}
	if (defined($list) and defined($dochr)){
		print "$USAGE\n";
		die "Error: Cannot do -genelist and -chr test concurrently\n";
	}
}
elsif( $ARGV[0] eq "-P" ) {
	print "Performing pathway-based analysis\n";
	print "Options in effect: \n";
	# Check if geneandp file exist and copy to working directory
	if( !-e "$geneandp" ) {
		print "$USAGE\n";
		die "Error: -geneandp file not found\n";
	}
#	copy( "$geneandp", "$time_stamp/geneandp" ) or die "Copy failed: $!";
	system("sort -k1,1 $geneandp >$time_stamp/geneandp"); # first column should be a gene id
	print "\t-geneandp $geneandp\n";
	# Check if glist file exist and copy to working directory
	if( !-e "$glist" ) {
		print "$USAGE\n";
		die "Error: -glist file not found. Please provide gene annotation file\n";
	}
#	copy( "$glist", "$time_stamp/tempglist" ) or die "Copy failed: $!";
	system("sort -k4,4 $glist >$time_stamp/tempglist");
	print "\t-glist $glist\n";
	# Check if geneandpath file exist and copy to working directory
	if( !-e "$geneandpath" ) {
		print "$USAGE\n";
		die "Error: -geneandpath file not found\n";
	}
#	copy( "$geneandpath", "$time_stamp/geneandpath" ) or die "Copy failed: $!";
	system("sort -k1,1 $geneandpath >$time_stamp/geneandpath"); # Make sure the first column is gene and second is pathway 
	print "\t-geneandpath $geneandpath\n";
	if( defined($outfile) ) {
		print "\t-out $outfile";
	}
	if( defined($maxsample)) {
		print "\t-maxsample $maxsample\n";
	}
}
else
{
	print "$USAGE\n";
	die "Error: First argument must be either -G or -P to perform gene-based test or pathway-based test respectively\n";
}
###--- END FRONT MATTER ---###

###--- BEGIN FILE VALIDATION and ANALYSIS ---###

#Chnging to working directory
print "\nChanging working directory...";
chdir "$time_stamp" or die "Error: cannot change to working directory: $!";
print "done\n";

# Check that R and plink exists
&checkplinkandr;
# Check that mvtnorm and corpcor exists
&checkpackages;

#Validating input files and running analysis
if( $ARGV[0] eq "-G" ) {
	&testingsnpandp;
	&testingglist;


	#Now mapping input variants to genotype file 
	system("awk '{print $1}' snpandp >snpIPlist");
	$input_variant_count = `wc -l <snpIPlist`;
	print "Number of variants provided in input file - $input_variant_count\n";
	system("plink --bfile $custom --extract snpIPlist --make-bed --out custom --noweb --silent > /dev/null");
	$input_variant_merged_with_genofile_count = `wc -l <custom.bim`;
	print "Number of variants available after merging with genotype file - $input_variant_merged_with_genofile_count";
#	system("plink --bfile $custom --make-set NEWglist --write-set --out Originalglistsetfile --noweb --silent > /dev/null");
	
	# Now make allgene files in geneset directory
	print "\nWritting allgene files...";
	mkdir "geneset", 0777 or die "Error: cannot make working directory. Check that vegas has correct permission settings: $!";

	for my $chromosome (1..24) {
		open(ALLGENES, ">geneset/allgene$chromosome");

		foreach ( @hglist ) {
	    	chomp($_);
	    	@each_gene = split(/\s/,$_);
			$given_chromosome = @each_gene[0];
			$given_symbol = @each_gene[3];
	
			if($given_chromosome == $chromosome) {
			print ALLGENES "$given_symbol\n";
			}
		}
		close ALLGENES;
	}

# Identify genes with more than 0 variants in it
	
	print "done.\n";
	
	
	# Do test
# Make one more file for genes with p-value less than 2e-6 reporting all clumps at r2 0.20

	if(defined($topten)){
		system("echo \"Chr Gene nSNPs nSims Start Stop Test Pvalue Top-$percentage-pvalue Best-SNP SNP-pvalue\" > gene-basedoutput.out");
	}
	if(defined($topsnp)){
		system("echo \"Chr Gene nSNPs nSims Start Stop Test Pvalue Top-SNP-pvalue Best-SNP SNP-pvalue\" > gene-basedoutput.out");
	}
	unless(defined($topsnp) || defined($topten)){ # Define vanilla test
		system("echo \"Chr Gene nSNPs nSims Start Stop Test Pvalue Best-SNP SNP-pvalue\" > gene-basedoutput.out");
	}
	
	if(defined($dochr)){
		&docustomchr;
	}
	if(defined($list)){
		&docustomlist;
	}
	if(!defined($list) and !defined($dochr)){
		&customtest;
	}
	if(defined($clustertest)){
		&clustertest;
	}
	if(defined($outfile)){
		system("cp gene-basedoutput.out ../$outfile.out");
	}
	else{
		system("cp gene-basedoutput.out ../gene-basedoutput.out");
	}	
}
elsif( $ARGV[0] eq "-P" ) { 
	&testinggeneandp;
	&testingglist;
	print "Merging geneandp, geneandpath and glist files...";
	system("join -1 1 -2 1 geneandp geneandpath >genePandpath"); 
	system("join -1 1 -2 4 genePandpath tempglist >genePpathandlocTEMP");
	system("awk '{print \$0,\$4 * 1e9 + \$5, \$4 * 1e9 + \$6}' genePpathandlocTEMP >genePpathandloc");
	print "done\n";
	
# Report statistics for compititive pathway analysis
	print "Following are the statistics for compititive pathway analysis:\n";
	system("awk '{print \$3}' genePpathandloc | sort -u >pathwaylist");
	$pathway_count = `wc -l <pathwaylist`;
	print "\tNumber of pathways to test - $pathway_count\n";
	system("awk '{print \$1,\$2}' genePpathandloc | sort -u >gandpforresampling");
	$gene_count_forsampling = `wc -l <gandpforresampling`;
	print "\tNumber of genes mapped to pathways - $gene_count_forsampling\n";
	
	&makepforsampling;
	
	system("echo \"Pathway nGenesMapped nGenesUsed nSamples ObservedP empiricalP Genes\" > pathway-basedoutput.out");
	
	&defaultpathwaytest;
	
	if(defined($outfile)){
		system("cp pathway-basedoutput.out ../$outfile.out");
	}
	else{
		system("cp pathway-basedoutput.out ../pathway-basedoutput.out");
	}
	
#	if( defined( pathanalysisforomicsdata ) ) {	 # Develop method for pathway analysis of omics data analysis
#	}
#	else { # be defualt run pathway analysis for GWAS data
#	system("echo \"Pathway nGenesMapped nGenesUsed nSamples ObservedP empiricalP Genes\" > pathway-basedoutput.out");
#	}
	
#	&testinggeneandpath; # Make these two sub routines
}

# Clean up files

chdir "../" or die "Error: Cannot change directory: $!";
unless(defined($keeptimestamp)){
	system("rm -R $time_stamp");
}
#rmdir $time_stamp or warn "cannot remove $time_stamp\n";

print "\nThank you for using VEGAS2 version 2!\n";

###--- END GENEBASED TEST ---###

###--- START SUBROUTINES ---###

# Now make a file for random sampling
sub makepforsampling
{

	$vanillatest=1;

	system("R --vanilla --slave <<EOF
	options(warn=-1)
	rm(list = ls(all = TRUE))

	Gene_and_P<-read.table('geneandp',header=F)
	X2_for_Sampling<-qchisq(Gene_and_P[,2],1,lower.tail=F)
	write.table(X2_for_Sampling,'x2forsampling',quote=F,col.names=F,row.names=F)

EOF"); # EOF should not have tab or spaces 

}

sub defaultpathwaytest
{
	print "\nRunning pathway analysis...\n";
	print "Pathway nGenesMapped nGenesUsed nSamples ObservedP empiricalP Genelist\n";
# Open pathwaylist
	open(ALLPATHWAYS, "pathwaylist") or die "Error: Cannot find pathwaylist\n";
	@allpathways = <ALLPATHWAYS>;
	close(ALLPATHWAYS);

	foreach $pathway (@allpathways) {
		chomp($pathway);
		system("grep -w \"$pathway\" genePpathandloc > genesinPath.position");
	
		# Test only those pathways with length more than 5 genes in it
		$pathway_length = `wc -l genesinPath.position`;
		if($pathway_length < 5){next;}
		system("sort -gk7,7 genesinPath.position >genesinPath.positionsorted");

		# Drpping all but one gene from cluster of neiboughring genes
		system("

		R --vanilla --slave <<EOF
		options(warn=-1)
# Now try doing it with data.table

		PATHWAY_DATA <- read.table('genesinPath.positionsorted',header=F,colClasses=c('character','numeric','character',rep('numeric',5)))
	
		Pathway_Name    <- PATHWAY_DATA[1,3]							## VARIABLE 1
		Original_Length <- length(PATHWAY_DATA[,1]) 				## VARIABLE 2

		tooclose <- rep(0,Original_Length)
		for (w in seq(1,(Original_Length-1))){
			if((PATHWAY_DATA[w+1,7] - PATHWAY_DATA[w,8]) < 500e3){
				tooclose[w] <- 1
			}
		}
	
		PATHWAY_DATA2 <- data.frame(PATHWAY_DATA,tooclose)
		
		PATHWAY_DATA3 <- PATHWAY_DATA2[PATHWAY_DATA2[,9]==0,]
	
		Neighbour_Filter_Length <- length(PATHWAY_DATA3[,1])				## VARIABLE 3
		
		if(Neighbour_Filter_Length > 1){
			ObservedX2	 <- sum(qchisq(PATHWAY_DATA3[,2],1,lower.tail=F))
			ObservedP    <- pchisq(ObservedX2, Neighbour_Filter_Length, lower.tail=F) ## VARIABLE 5
			
			Gene_List    <- paste(PATHWAY_DATA3[,1],collapse='_')		## VARIABLE 7
			Temp_Pathway <- data.frame(Pathway_Name,Original_Length,Neighbour_Filter_Length,ObservedX2,ObservedP,Gene_List)
			write.table(Temp_Pathway,'temppathway',row.names=F,col.names=F,quote=F)
		}
	
EOF");

		# Now compute empirical p-value
		open(WORKINGPATHWAY, "temppathway") or die "Error: Cannot find temppathway\n";
		@workingpathway = <WORKINGPATHWAY>;
		close(WORKINGPATHWAY);
		
		open(X2FORSAMPLING, "x2forsampling") or die "Error: Cannot find P_for_Sampling\n";
		@x2forsampling = <X2FORSAMPLING>;
		close(X2FORSAMPLING);
				
		foreach ( @workingpathway ) {
		    chomp($_);
		    @each_pathway_data = split(/\s/,$_);
			
			$pathway_name = @each_pathway_data[0];
			$pathway_original_length = @each_pathway_data[1];
			$nogenes_to_resample = @each_pathway_data[2];
			$observedX2 = @each_pathway_data[3];
			$pathway_observedP = @each_pathway_data[4];
			$pathway_listofgenes = @each_pathway_data[5];
			
			print "$pathway_name $pathway_original_length $nogenes_to_resample ";
						
			die "Too few elements (".scalar(@x2forsampling).") to select $nogenes_to_resample from\n"
			        unless $nogenes_to_resample < @x2forsampling;
			
			# Shuffled list of indexes into @deck
			# Using empP formula 1+r/1+n
			$Number_of_reps_above_observed = 1; 
			$number_of_shuffle = 1000;
							
			for my $shuffleCOUNT (1 .. $number_of_shuffle) {
				@shuffled_indexes = shuffle(0..$#x2forsampling);
				# Get just N of them.
				@pick_indexes = @shuffled_indexes[ 0 .. $nogenes_to_resample - 1 ];  
				# Pick cards from array
				@picks = @x2forsampling[ @pick_indexes ];
				
				if(sum(@picks) >= $observedX2 ){
					$Number_of_reps_above_observed = $Number_of_reps_above_observed + 1;
				}
			}
			$empiricalP = $Number_of_reps_above_observed / (1 + $number_of_shuffle);			
			
			if($empiricalP < 0.01){
				$number_of_shuffle = 10000;
				for my $shuffleCOUNT (1 .. $number_of_shuffle) {
					@shuffled_indexes = shuffle(0..$#x2forsampling);
					# Get just N of them.
					@pick_indexes = @shuffled_indexes[ 0 .. $nogenes_to_resample - 1 ];  
					# Pick cards from array
					@picks = @x2forsampling[ @pick_indexes ];
				
					if(sum(@picks) >= $observedX2 ){
						$Number_of_reps_above_observed = $Number_of_reps_above_observed + 1;
					}
				}
				$empiricalP = $Number_of_reps_above_observed / (1 + $number_of_shuffle);
			}
			
			if($empiricalP < 0.001){
				$number_of_shuffle = 50000;		
				for my $shuffleCOUNT (1 .. $number_of_shuffle) {
					@shuffled_indexes = shuffle(0..$#x2forsampling);
					# Get just N of them.
					@pick_indexes = @shuffled_indexes[ 0 .. $nogenes_to_resample - 1 ];  
					# Pick cards from array
					@picks = @x2forsampling[ @pick_indexes ];
				
					if(sum(@picks) >= $observedX2 ){
						$Number_of_reps_above_observed = $Number_of_reps_above_observed + 1;
					}
				}
				$empiricalP = $Number_of_reps_above_observed / (1 + $number_of_shuffle);
			}
			
			if($empiricalP < 0.0005){
				$number_of_shuffle = 100000;				
				for my $shuffleCOUNT (1 .. $number_of_shuffle) {
					@shuffled_indexes = shuffle(0..$#x2forsampling);
					# Get just N of them.
					@pick_indexes = @shuffled_indexes[ 0 .. $nogenes_to_resample - 1 ];  
					# Pick cards from array
					@picks = @x2forsampling[ @pick_indexes ];
				
					if(sum(@picks) >= $observedX2 ){
						$Number_of_reps_above_observed = $Number_of_reps_above_observed + 1;
					}
				}
				$empiricalP = $Number_of_reps_above_observed / (1 + $number_of_shuffle);
			}
			
			if($empiricalP < 0.0001){
				$number_of_shuffle = 500000;				
				for my $shuffleCOUNT (1 .. $number_of_shuffle) {
					@shuffled_indexes = shuffle(0..$#x2forsampling);
					# Get just N of them.
					@pick_indexes = @shuffled_indexes[ 0 .. $nogenes_to_resample - 1 ];  
					# Pick cards from array
					@picks = @x2forsampling[ @pick_indexes ];
				
					if(sum(@picks) >= $observedX2 ){
						$Number_of_reps_above_observed = $Number_of_reps_above_observed + 1;
					}
				}
				$empiricalP = $Number_of_reps_above_observed / (1 + $number_of_shuffle);
			}
			
			if($empiricalP < 0.00002){
				$number_of_shuffle = 1000000;
				if(defined($maxsample)){
					$number_of_shuffle = $maxsample;
				}
							
				for my $shuffleCOUNT (1 .. $number_of_shuffle) {
					@shuffled_indexes = shuffle(0..$#x2forsampling);
					# Get just N of them.
					@pick_indexes = @shuffled_indexes[ 0 .. $nogenes_to_resample - 1 ];  
					# Pick cards from array
					@picks = @x2forsampling[ @pick_indexes ];
				
					if(sum(@picks) >= $observedX2 ){
						$Number_of_reps_above_observed = $Number_of_reps_above_observed + 1;
					}
				}
				$empiricalP = $Number_of_reps_above_observed / (1 + $number_of_shuffle);
			}

		# Now print results #Majority of variables defined at start of foreach loop
		open(RESULTSPATHWAY, ">>pathway-basedoutput.out");
		print RESULTSPATHWAY "$pathway_name $pathway_original_length $nogenes_to_resample $number_of_shuffle $pathway_observedP $empiricalP $pathway_listofgenes\n";
		print "$number_of_shuffle $pathway_observedP $empiricalP $pathway_listofgenes\n";

		} # Foreach @workingpathway
	} # foreach $pathway
	print "\ndone\n";
} # sub defaultpathwaytest


sub customtest{ # VEGAS using custom set of SNPs - default
	print "\nReading custom genotypes...done\n\n";
	unless(defined($topsnp) || defined($topten)){ # Define vanilla test
		$vanillatest=1;
	}
	
	for ($chr = 1; $chr <= 22; $chr++) {
	#for ($chr=22){
#		system("plink --bfile ../$custom --chr $chr --make-bed --out custom$chr --noweb --silent > /dev/null");
		
		open(ALLGENES, "geneset/allgene$chr") or die "Error: Cannot find geneset/allgene$chr\n";
		@allgenes = <ALLGENES>;
		close(ALLGENES);
		
		foreach $gene (@allgenes){

			if(-e "tempgene.snp"){system("rm tempgene.snp");}
			if(-e "tempgene.pvalue"){system("rm tempgene.pvalue");}
			if(-e "gene.position"){system("rm gene.position");}
			
			chomp($gene);
			system("grep -w \"$gene\" tempglist|head -1 > gene.position");
			
			@glistpos = grep(/\b$gene\b/, @hglist);	
			chomp(@glistpos);
			@glistpos = split(/ /,@glistpos[0]);
		
			if(!defined($lower)){
				$start = @glistpos[1]-50000;
			}
			else{
				$start = @glistpos[1]-$lower;
			}
			if(!defined($upper)){
				$stop = @glistpos[2]+50000;
			}
			else{
				$stop = @glistpos[2]+$upper;
			}
	
#			system("plink --bfile custom$chr --chr $chr --from-bp $start --to-bp $stop --write-snplist --noweb --silent > /dev/null"); # It would be super if you could do these three lines with one PLINK command
			system("plink --bfile $custom --chr $chr --from-bp $start --to-bp $stop --write-snplist --noweb --silent > /dev/null");

			if(-e "plink.snplist"){
				system("mv plink.snplist tempgene.snp");
			}
			
			open(SNPLIST,"tempgene.snp");
			@snplist = <SNPLIST>;
			close(SNPLIST);
			
			if(scalar(@snplist == 0)){next;}
			
#			foreach $snp(@snplist){
#				chomp($snp);
#				system("grep -w \"$snp\" snpandp >> tempgene.pvalue");
#			}

			system("sort -k1,1 tempgene.snp >tempgene.snp2");
			system("join -1 1 -2 1 tempgene.snp2 snpandp >tempgene.pvalue");
			
			# Send to R to do make into tempgene format, then send to &hapsmims? - files prepared: tempgene.pvalue plink.ld gene.position
			
			&maketempgener;
			&hapsims;
#			system("rm custom$chr.bed custom$chr.bim custom$chr.fam");
		}
	}
}

sub docustomchr { # VEGAS using custom set of SNPs - chromosome
	print "\nReading custom genotypes...done\n\n";
	unless(defined($topsnp) || defined($topten)){ # Define vanilla test
		$vanillatest=1;
	}
	

	#for ($chr = 1; $chr <= 22; $chr++) {
	for ($chr=$dochr){
#		system("plink --bfile ../$custom --chr $chr --make-bed --out custom$chr --noweb --silent > /dev/null");
		
		open(ALLGENES, "geneset/allgene$chr");
		@allgenes = <ALLGENES>;
		close(ALLGENES);
		
		foreach $gene (@allgenes){

			if(-e "tempgene.snp"){system("rm tempgene.snp");}
			if(-e "tempgene.pvalue"){system("rm tempgene.pvalue");}
			if(-e "gene.position"){system("rm gene.position");}
			
			chomp($gene);
			system("grep -w \"$gene\" tempglist|head -1 > gene.position");
			
			@glistpos = grep(/\b$gene\b/, @hglist);	
			chomp(@glistpos);
			@glistpos = split(/ /,@glistpos[0]);
		
			if(!defined($lower)){
				$start = @glistpos[1]-50000;
			}
			else{
				$start = @glistpos[1]-$lower;
			}
			if(!defined($upper)){
				$stop = @glistpos[2]+50000;
			}
			else{
				$stop = @glistpos[2]+$upper;
			}
	
#			system("plink --bfile custom$chr --chr $chr --from-bp $start --to-bp $stop --write-snplist --noweb --silent > /dev/null"); # It would be super if you could do these three lines with one PLINK command
			system("plink --bfile $custom --chr $chr --from-bp $start --to-bp $stop --write-snplist --noweb --silent > /dev/null"); 
			if(-e "plink.snplist"){system("mv plink.snplist tempgene.snp");}
			
			open(SNPLIST,"tempgene.snp");
			@snplist = <SNPLIST>;
			close(SNPLIST);
			
			if(scalar(@snplist == 0)){next;}
			
#			foreach $snp(@snplist){
#				chomp($snp);
#				system("grep -w \"$snp\" snpandp >> tempgene.pvalue");
#			}
			system("sort -k1,1 tempgene.snp >tempgene.snp2");
			system("join -1 1 -2 1 tempgene.snp2 snpandp >tempgene.pvalue");
			
			# Send to R to do make into tempgene format, then send to &hapsmims? - files prepared: tempgene.pvalue plink.ld gene.position
			
			&maketempgener;
			&hapsims;
		}
	}
}

sub docustomlist { # VEGAS using custom set of SNPs - genelist
	print "\nReading custom genotypes...done\n\n";
	unless(defined($topsnp) || defined($topten)){ # Define vanilla test
		$vanillatest=1;
	}
	
	open(HG18, "tempglist");
	@hg18 = <HG18>;
	close HG18;
	
	open(GENELIST, "customgenelist");
	@genelist = <GENELIST>;
	close GENELIST;
	
	# Read in genelist and match to chromosomes
	foreach $genelist (@genelist){
		chomp($genelist);
		$genelist =~ s/^\s+//;
		$genelist =~ s/\s+$//;
		foreach $hg18 (@hg18){
			chomp($hg18);
			$hg18 =~ s/^\s+//;
			$hg18 =~ s/\s+$//;
			@hglist = split(/ /,$hg18);
			if($genelist eq @hglist[3]){
				push(@dochroms,"@hglist[0]\n");
				push(@dogenes,"@hglist[3]\n");
			
				open(WRITECHR, ">>allgene@hglist[0]");
				print WRITECHR "$genelist\n";
				close WRITECHR;			
			}
		}
	}
	
	# Extract unique chromosomes
	%seen = ();
	foreach $item (@dochroms){
		push(@uniqchr, $item) unless $seen{$item}++;
	}


	foreach $chr(@uniqchr) {
	#for ($chr=22){
		chomp($chr);
		system("plink --bfile $custom --chr $chr --make-bed --out custom$chr --noweb --silent > /dev/null");
		
		open(ALLGENES, "allgene$chr");
		@allgenes = <ALLGENES>;
		close(ALLGENES);
		
		foreach $gene (@allgenes){

			if(-e "tempgene.snp"){system("rm tempgene.snp");}
			if(-e "tempgene.pvalue"){system("rm tempgene.pvalue");}
			if(-e "gene.position"){system("rm gene.position");}
			
			chomp($gene);
			system("grep -w \"$gene\" tempglist|head -1 > gene.position");
			
			@glistpos = grep(/\b$gene\b/, @hg18);	
			chomp(@glistpos);
			@glistpos = split(/ /,@glistpos[0]);
		
			if(!defined($lower)){
				$start = @glistpos[1]-50000;
			}
			else{
				$start = @glistpos[1]-$lower;
			}
			if(!defined($upper)){
				$stop = @glistpos[2]+50000;
			}
			else{
				$stop = @glistpos[2]+$upper;
			}
	
			system("plink --bfile custom$chr --chr $chr --from-bp $start --to-bp $stop --write-snplist --noweb --silent > /dev/null");
#			system("plink --bfile custom --chr $chr --from-bp $start --to-bp $stop --write-snplist --noweb --silent > /dev/null");
			if(-e "plink.snplist"){system("mv plink.snplist tempgene.snp");}

# Check if extra variants are provided add the same to docustomchr and customtest
			
			open(SNPLIST,"tempgene.snp");
			@snplist = <SNPLIST>;
			close(SNPLIST);
			
			if(scalar(@snplist == 0)){next;}
			
			system("sort -k1,1 tempgene.snp >tempgene.snp2");
			system("join -1 1 -2 1 tempgene.snp2 snpandp >tempgene.pvalue");
	
# Send to R to do make into tempgene format, then send to &hapsmims? - files prepared: tempgene.pvalue plink.ld gene.position
			
			&maketempgener;
# Check the line count after removing duplicate SNPs
			$line_count = `wc -l <tempgene.pvalue`;
			if ( $line_count > 1 ) {
				&hapsims;
			}
			else{
			next;
			}
		}
	}
}

sub maketempgener{
unless(-z "tempgene.pvalue"){

	if(defined($prunersquare)){
	$line_count = `wc -l <tempgene.pvalue`;
		
		if ( $line_count > 1 ) {
#		system("awk '{print \$1}' tempgene.pvalue >prune.input");
#		system("awk 'BEGIN{print \"SNP P\"}1' tempgene.pvalue >clump.input");
		system("echo 'SNP P' >clump.header");
		system("cat clump.header tempgene.pvalue >clump.input");
#		system("plink --bfile custom$chr --extract prune.input --indep-pairwise 100 5 $prunersquare --noweb --silent");
# Clump would work better than the prune sice it will choose the variant based on the p-values insted of randomly 
#		system("plink --bfile custom$chr --clump clump.input --clump-p1 1 --clump-p2 1 --clump-r2 0.99 --clump-kb 1000 --noweb --silent");
#		system("plink --bfile ../$custom --clump clump.input --clump-p1 1 --clump-p2 1 --clump-r2 0.99 --clump-kb 1000 --noweb --silent > /dev/null");
		system("plink --bfile custom --clump clump.input --clump-p1 1 --clump-p2 1 --clump-r2 0.99 --clump-kb 1000 --noweb --silent > /dev/null");
		
#		system("sort -k1,1 plink.prune.in >plink.prune.in2");
#		system("sort -k1,1 tempgene.pvalue >tempgene.pvalue2");
#		system("join -1 1 -2 1 plink.prune.in2 tempgene.pvalue2 >tempgene.pvalue");
		system("awk 'NR>1{print \$3,\$5}' plink.clumped >tempgene.pvalue");
		}
	}

	system("

R --vanilla --slave <<EOF
options(warn=-1)
pvals <- read.table('tempgene.pvalue',header=F)
positions <- read.table('gene.position',header=F)
tempgene <- data.frame(pvals[,1],rep(positions[,4],length(pvals[,1])),rep(positions[,1],length(pvals[,1])),rep(positions[,2],length(pvals[,1])),rep(positions[,3],length(pvals[,1])),qchisq(pvals[,2],1,lower.tail=F))
tempgene2<- subset(tempgene,!duplicated(tempgene[,1]))
write.table(tempgene2,'tempgene',row.names=F,col.names=F,quote=F)

EOF");

system("awk '{print \$1;}' tempgene > tempgene.snp");
}

# if provided summary concatanate tempgene files

}
	
sub hapsims {
	if (-e "needmoresims"){system("rm needmoresims");}
	if (-e "correctedp.txt"){system("rm correctedp.txt");}
	if (-e "plink.ld"){system("rm plink.ld");}
	if (-e "do100000sims"){system("rm do100000sims");}
	if (-e "alwayswritep"){system("rm alwayswritep");}
	
	&plink;
	
	# 1000 sims
	if(defined($vanillatest)){
		&mvsimsr(1000);
	}
	if(defined($topten)){
		&mvsimstoptenr(1000);
	}
	if(defined($topsnp)){
		&mvsimstopsnpr(1000);
	}
	
	
	if (-e "correctedp.txt"){ # If 1000 sims is enough
		open(CORRECTEDP, "correctedp.txt");
		$correctedp = <CORRECTEDP>;
		open(PRINTCORRECTEDP, ">>gene-basedoutput.out");
		print PRINTCORRECTEDP "$correctedp";
		print "$correctedp";
	}
	
	if (-e "do100000sims"){ # Do 100000 sims
		if(defined($vanillatest)){
			&mvsimsr(100000);
		}
		if(defined($topten)){
			&mvsimstoptenr(100000);
		}
		if(defined($topsnp)){
			&mvsimstopsnpr(100000);
		}
		
		if(-e "correctedp.txt"){ # If p>0.1 (not likely) write output
			open(CORRECTEDP, "correctedp.txt");
			$correctedp = <CORRECTEDP>;
			open(PRINTCORRECTEDP, ">>gene-basedoutput.out");
			print PRINTCORRECTEDP "$correctedp";
			print "$correctedp";
		}
		if(-e "do100000sims"){ # If 0.001 < p < 0.1, write output
			open(CORRECTEDP, "do100000sims");
			$correctedp = <CORRECTEDP>;
			open(PRINTCORRECTEDP, ">>gene-basedoutput.out");
			print PRINTCORRECTEDP "$correctedp";
			print "$correctedp";
		}
		if(-e "needmoresims"){ #Do max sims, then write output
			if(defined($max)){
				if(defined($vanillatest)){
					&mvsimsr($max);
				}
				if(defined($topten)){
					&mvsimstoptenr($max);
				}
				if(defined($topsnp)){
					&mvsimstopsnpr($max);
				}
			}
			else{
				if(defined($vanillatest)){
					&mvsimsr(1000000);
				}
				if(defined($topten)){
					&mvsimstoptenr(1000000);
				}
				if(defined($topsnp)){
					&mvsimstopsnpr(1000000);
				}
			}
			open(CORRECTEDP, "alwayswritep");
## Run clumping for genes with p-value less than 2e-6 reporting all clumps at r2 0.20
			$correctedp = <CORRECTEDP>;
			open(PRINTCORRECTEDP, ">>gene-basedoutput.out");
			print PRINTCORRECTEDP "$correctedp";
			print "$correctedp";
			
		}
	}
	
	if (-e "needmoresims"){ # Do max sims
		if(defined($max)){
			if(defined($vanillatest)){
				&mvsimsr($max);
			}
			if(defined($topten)){
				&mvsimstoptenr($max);
			}
			if(defined($topsnp)){
				&mvsimstopsnpr($max);
			}
		}
		else{
			if(defined($vanillatest)){
				&mvsimsr(1000000);
			}
			if(defined($topten)){
				&mvsimstoptenr(1000000);
			}
			if(defined($topsnp)){
				&mvsimstopsnpr(1000000);
			}
		}
		open(CORRECTEDP, "alwayswritep");
		$correctedp = <CORRECTEDP>;
		open(PRINTCORRECTEDP, ">>gene-basedoutput.out");
		print PRINTCORRECTEDP "$correctedp";
		print "$correctedp";
	}
}

sub mvsimsr{

	if (-e "needmoresims"){system("rm needmoresims");}
	if (-e "correctedp.txt"){system("rm correctedp.txt");}
	if (-e "do100000sims"){system("rm do100000sims");}
	if (-e "alwayswritep"){system("rm alwayswritep");}	

system("

R --vanilla --slave  <<EOF

options(warn=-1)

rm(list = ls(all = TRUE))

#read in ld matrix and generate lots of null chi distns based on correlated snps

#library(mvtnorm)
numsnps <- $numsnps
matrix(scan('plink.ld',quiet=T),nc=numsnps) -> co
#co <- sqrt(co)
#head(co)

#check that co is positive definite. Make diagonals 1.0001.

library(corpcor)

reps <- $_[0]


if(is.positive.definite(co)==F){
	co <- make.positive.definite(co)
}
if(is.positive.definite(co)==F){
	matrix(scan('plink.ld',quiet=T),nc=numsnps) -> co
	for(i in 1:numsnps){
		co[i,i] <- 1.0001
	}
}
if(is.positive.definite(co)==F){
	for(i in 1:numsnps){
		co[i,i] <- 1.001
	}
}
if(is.positive.definite(co)==F){
	for(i in 1:numsnps){
		co[i,i] <- 1.01
	}
}


library(mvtnorm)
rmvnorm(reps,mean=rep(0,numsnps),sigma=co) -> rd

rowSums(rd ^2) -> summedsq


#write.table(data.frame(summedsq,summedsqtopten),file='geneXsummedsq.txt',row.names=F,quote=F,col.names=F)

read.table('tempgene',header=F) -> gene
ass <- gene[,6]
#ass
sum(ass,na.rm=T) -> sumass
sum(!is.na(ass)) -> len # get number of snps with non-missing test stat
length(ass) -> asslen

#scan('geneXsummedsq.txt') -> simvals # read in simulation replicates

#head(simvals)
(1 + sum(ifelse(sumass > summedsq, 0, 1)))/(1 + length(summedsq)) -> empp
#empp

snpp <- pchisq(max(ass),1,lower.tail=F)
snpname <- gene[which(ass==max(ass)),1][1]

if(empp>0.1){
write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,snpname,snpp),'correctedp.txt',row.names=F,col.names=F)
}

if(empp<=0.1 & empp>0.001){
write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,snpname,snpp),'do100000sims',row.names=F,col.names=F)
}

if(empp<=0.001){
write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,snpname,snpp),'needmoresims',row.names=F,col.names=F)
}

write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,snpname,snpp),'alwayswritep',row.names=F,col.names=F)

EOF
");
}

sub mvsimstoptenr{


if (-e "needmoresims"){system("rm needmoresims");}
if (-e "correctedp.txt"){system("rm correctedp.txt");}
if (-e "do100000sims"){system("rm do100000sims");}
if (-e "alwayswritep"){system("rm alwayswritep");}

if(-e "tempempp"){system("rm tempempp")}
if(-e "temprd2"){system("rm temprd2")}
if(-e "temprd-sorted"){system("rm temprd-sorted")}

if($_[0] == 1000){

system("R --vanilla --slave <<EOF

options(warn=-1)
rm(list = ls(all = TRUE))

#read in ld matrix and generate lots of null chi distns based on correlated snps

#library(mvtnorm)
$numsnps -> numsnps
matrix(scan('plink.ld',quiet=T),nc=numsnps) -> co
#co <- sqrt(co)
#head(co)

#check that co is positive definite. Make diagonals 1.0001.

library(corpcor)

reps <- $_[0]


if(is.positive.definite(co)==F){
	co <- make.positive.definite(co)
}
if(is.positive.definite(co)==F){
	matrix(scan('plink.ld',quiet=T),nc=numsnps) -> co
	for(i in 1:numsnps){
		co[i,i] <- 1.0001
	}
}
if(is.positive.definite(co)==F){
	for(i in 1:numsnps){
		co[i,i] <- 1.001
	}
}
if(is.positive.definite(co)==F){
	for(i in 1:numsnps){
		co[i,i] <- 1.01
	}
}


library(mvtnorm)
rmvnorm(reps,mean=rep(0,numsnps),sigma=co) -> rd

rowSums(rd ^2) -> summedsq

percentage <- $percentage

#Old bug in all previous version highlited by Julian Hecker et al
#topten <- function(x){
#	sum(sort(x,decreasing=T)[1:(1+length(diag(co))*percentage)]^2)
#}

# Update in version 3
topten <- function(x){
	sum(sort(x^2,decreasing=T)[1:(1+length(diag(co))*percentage)])
}


apply(rd,1,topten) -> summedsqtopten

read.table('tempgene',header=F) -> gene
ass <- gene[,6]
#ass
sum(ass,na.rm=T) -> sumass
sum(!is.na(ass)) -> len # get number of snps with non-missing test stat
length(ass) -> asslen

asstopten <- sort(ass,decreasing=T)[1:(1+length(diag(co))*percentage)]
sumasstopten <- sum(asstopten)
(1 + sum(ifelse(sumasstopten > summedsqtopten, 0, 1)))/(1 + length(summedsqtopten)) -> empptopten
#empptopten

(1 + sum(ifelse(sumass > summedsq, 0, 1)))/(1 + length(summedsq)) -> empp
#empp

snpp <- pchisq(max(ass),1,lower.tail=F)
snpname <- gene[which(ass==max(ass)),1][1]

if(empptopten >= 0.1){
write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'correctedp.txt',row.names=F,col.names=F)
}

#if(empp < 0.1 & empp > 0.001){
#write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'do100000sims',row.names=F,col.names=F)
#}

if(empptopten < 0.1 & empptopten > 0.001){
write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'do100000sims',row.names=F,col.names=F)
}

if(empptopten <= 0.001){
write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'needmoresims',row.names=F,col.names=F)
}

write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'alwayswritep',row.names=F,col.names=F)

EOF");
}


if($_[0] > 1000){ # Do 1e5 simulations and topten sort etc, save to disk, then repeat until reach 1e6 (or $max)
	my $reps = int($_[0]/100000);
	for($i = 1; $i <= $reps; $i++){
		system("R --vanilla --slave <<EOF
numsnps <- $numsnps
reps <- 100000
percentage <- $percentage
options(warn=-1)

matrix(scan('plink.ld',quiet=T),nc=numsnps) -> co

library(corpcor)

if(is.positive.definite(co)==F){
	co <- make.positive.definite(co)
}
if(is.positive.definite(co)==F){
	matrix(scan('plink.ld',quiet=T),nc=numsnps) -> co
	for(i in 1:numsnps){
		co[i,i] <- 1.0001
	}
}
if(is.positive.definite(co)==F){
	for(i in 1:numsnps){
		co[i,i] <- 1.001
	}
}
if(is.positive.definite(co)==F){
	for(i in 1:numsnps){
		co[i,i] <- 1.01
	}
}

library(mvtnorm)
rmvnorm(reps,mean=rep(0,numsnps),sigma=co) -> rd
rowSums(rd ^2) -> summedsq

#Old
#topten <- function(x){
#	sum(sort(x,decreasing=T)[1:(1+length(diag(co))*percentage)]^2)
#}
#Update in version 3
topten <- function(x){
	sum(sort(x^2,decreasing=T)[1:(1+length(diag(co))*percentage)])
}

apply(rd,1,topten) -> summedsqtopten

write.table(summedsq,'vanilla',append=T,col.names=F,row.names=F,quote=F)
write.table(summedsqtopten,'topten',append=T,col.names=F,row.names=F,quote=F)
EOF")
	}
	system("R --vanilla --slave <<EOF
read.table('tempgene',header=F) -> gene
summedsqtopten <- scan('topten')
summedsq <- scan('vanilla')
percentage <- $percentage
reps <- $_[0]

ass <- gene[,6]
#ass
sum(ass,na.rm=T) -> sumass
sum(!is.na(ass)) -> len # get number of snps with non-missing test stat
length(ass) -> asslen

asstopten <- sort(ass,decreasing=T)[1:(1+length(gene[,1])*percentage)]
sumasstopten <- sum(asstopten)
(1 + sum(ifelse(sumasstopten > summedsqtopten, 0, 1)))/(1 + length(summedsqtopten)) -> empptopten
#empptopten

(1 + sum(ifelse(sumass > summedsq, 0, 1)))/(1 + length(summedsq)) -> empp
#empp

snpp <- pchisq(max(ass),1,lower.tail=F)
snpname <- gene[which(ass==max(ass)),1][1]

if(empptopten >= 0.1){
write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'correctedp.txt',row.names=F,col.names=F)
}

#if(empp < 0.1 & empp > 0){
#write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'do100000sims',row.names=F,col.names=F)
#}

if(empptopten < 0.1 & empptopten > 0.001){
write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'do100000sims',row.names=F,col.names=F)
}

if(empptopten < 0.001){
write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'needmoresims',row.names=F,col.names=F)
}

write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'alwayswritep',row.names=F,col.names=F)

EOF");

system("rm vanilla topten");

}
}

sub mvsimstopsnpr{
	if (-e "needmoresims"){system("rm needmoresims");}
	if (-e "correctedp.txt"){system("rm correctedp.txt");}
	if (-e "do100000sims"){system("rm do100000sims");}
	if (-e "alwayswritep"){system("rm alwayswritep");}
	if(-e "tempempp"){system("rm tempempp")}
	if(-e "temprd2"){system("rm temprd2")}
	if(-e "temprd-sorted"){system("rm temprd-sorted")}

	system("R --vanilla --slave <<EOF
	options(warn=-1)
	rm(list = ls(all = TRUE))

#read in ld matrix and generate lots of null chi distns based on correlated snps
	$numsnps -> numsnps
	matrix(scan('plink.ld',quiet=T),nc=numsnps) -> co

#check that co is positive definite. Make diagonals 1.0001.
	library(corpcor)
	reps <- $_[0]
	if(is.positive.definite(co)==F){
		co <- make.positive.definite(co)
	}
	if(is.positive.definite(co)==F){
		matrix(scan('plink.ld',quiet=T),nc=numsnps) -> co
		for(i in 1:numsnps){
			co[i,i] <- 1.0001
		}
	}
	if(is.positive.definite(co)==F){
		for(i in 1:numsnps){
			co[i,i] <- 1.001
		}
	}
	if(is.positive.definite(co)==F){
		for(i in 1:numsnps){
			co[i,i] <- 1.01
		}
	}

	library(mvtnorm)
	rmvnorm(reps,mean=rep(0,numsnps),sigma=co) -> rd
	rowSums(rd ^2) -> summedsq
	percentage <- $percentage
	topsnps <- function(x){
		max(x^2) -> topten
	}
	apply(rd,1,topsnps) -> summedsqtopten

	read.table('tempgene',header=F) -> gene
	ass <- gene[,6]
	sum(ass,na.rm=T) -> sumass
	sum(!is.na(ass)) -> len # get number of snps with non-missing test stat
	length(ass) -> asslen
	sumasstopten <- max(ass)
	(1 + sum(ifelse(sumasstopten > summedsqtopten, 0, 1)))/(1 + length(summedsqtopten)) -> empptopten
	(1 + sum(ifelse(sumass > summedsq, 0, 1)))/(1 + length(summedsq)) -> empp
	snpp <- pchisq(max(ass),1,lower.tail=F)
	snpname <- gene[which(ass==max(ass)),1][1]
	if(empptopten >= 0.1){
		write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'correctedp.txt',row.names=F,col.names=F)
	}
#	if(empp < 0.1 & empp > 0.001){
#		write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'do100000sims',row.names=F,col.names=F)
#	}
	if(empptopten < 0.1 & empptopten > 0.001){
		write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'do100000sims',row.names=F,col.names=F)
	}
	if(empptopten <= 0.001){
		write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'needmoresims',row.names=F,col.names=F)
	}
	write.table(data.frame('$chr','$gene',asslen,reps,unique(gene[4]),unique(gene[5]),sumass,empp,empptopten,snpname,snpp),'alwayswritep',row.names=F,col.names=F)

EOF");

}


sub sort_by_number{ $b <=> $a }

sub plink{ # Use plink to generate ld matrix, then check for errors
	
#	system("plink --bfile custom$chr --extract tempgene.snp --matrix --r --noweb --silent > /dev/null");
#	system("plink --bfile ../$custom --extract tempgene.snp --matrix --r --noweb --silent > /dev/null");
	system("plink --bfile custom --extract tempgene.snp --matrix --r --noweb --silent > /dev/null");

	open(PLINKLD, "plink.ld");
	@plinkld = <PLINKLD>;
	foreach $plinkld(@plinkld){
		$plinkld =~ s/nan/0/g; #replace nan's with 0 in ld matrix
	}
	if(scalar(@plinkld == 1)){
		@plinkld = 1;
	}
	
	open(PRINTPLINKLD, ">plink.ld");
	print PRINTPLINKLD "@plinkld";
	$numsnps = scalar(@plinkld);
	close PRINTPLINKLD;
	close PLINKLD;
	
}

sub checkplinkandr {
	print "\nChecking if plink and R exist...";
	system("echo -n \$PATH > binarypath");

	open(BINPATH, "binarypath");
	$binpath = <BINPATH>;
	@binpath = split(/:/, $binpath);

	$Rtimes=0;
	$plinktimes=0;

	foreach $i(@binpath){
		#print "$i\n";
		if(-e "$i/R"){
			$Rtimes=$Rtimes+1;
		}
		if(-e "$i/r"){
                        $Rtimes=$Rtimes+1;
                }
		if(-e "$i/plink"){
			$plinktimes=$plinktimes+1;
		}
	}
	if($Rtimes==0){
		chdir "../" or die "cannot change directory: $!";
		system("rm -R $time_stamp");
		die "Error: Cannot find R - make sure it can be accessed through \$PATH \n$!\n";
	}
	if($plinktimes==0){
		chdir "../" or die "cannot change directory: $!";
		system("rm -R $time_stamp");
		die "Error: Cannot find plink - make sure it can be accessed through \$PATH \n$!\n";
	}

	close PATH;
	system("rm binarypath");
	print "done\n";
}

sub checkpackages {
	print "\nChecking if mvtnorm and corpcor R packages exist...";
	system("
R --vanilla --slave <<EOF
options(warn=-1)
write.table(.packages(all.available=T),'R-packagelist',quote=F,col.names=F,row.names=F)
EOF");

	system("grep -w \"mvtnorm\" R-packagelist >> testpackagelist");
	system("grep -w \"corpcor\" R-packagelist >> testpackagelist");
	open(TESTPACKAGELIST,"testpackagelist");
	my @testpackagelist = <TESTPACKAGELIST>;
	if(scalar(@testpackagelist ne 2)){
		die "Error: Missing R library: VEGAS requires corpcor and mvtnorm to run\n";
	}
	close TESTPACKAGELIST;
	print "done\n";
}

sub testinggeneandp
{

# Validations

# Check pvalues in file are valid - in decimal or standard scientific notation (0.0034, 1e-6 etc...)
# Check if p-values are between 0 and 1
print "\nValidating input file...";

open(PVALS, "geneandp"); # Location of GWAS results file
@pvals = <PVALS>;

foreach ( @pvals ) {
    chomp($_);
    @each_line = split(/\s/,$_);
    if ( looks_like_number(@each_line[1] ) ) {
        $checked_linenumber = 0;
        if ( @each_line[1] < 0 || @each_line[1] > 1 ) {
            print WRITELOG "Error: Invalid p-values in $user_filename detected.\n P-values should be in between 0 to 1";
			print WRITELOG "Please check the column 2 (p-value) for @each_line\n";
			print "\nPlease check column 2 (p-value) for @each_line\n";
            die "Error: Invalid p-values in $user_filename detected.\n";
        }
    }
    else {
        print WRITELOG "Error: Invalid p-values in $user_filename detected.\n";
		print WRITELOG "Please check the column 2 (p-value) for @each_line\n";
		print "\nPlease check column 2 (p-value) for @each_line\n";
        die "Error: Invalid p-values in $user_filename detected.\n";
    }
}
close PVALS;
print "done.\n";

}

sub testingsnpandp
{

# Validations

# Check pvalues in file are valid - in decimal or standard scientific notation (0.0034, 1e-6 etc...)
# Check if p-values are between 0 and 1
print "\nValidating input file...";

open(PVALS, "snpandp"); # Location of GWAS results file
@pvals = <PVALS>;

foreach ( @pvals ) {
    chomp($_);
    @each_line = split(/\s/,$_);
    if ( looks_like_number(@each_line[1] ) ) {
        $checked_linenumber = 0;
        if ( @each_line[1] < 0 || @each_line[1] > 1 ) {
            print WRITELOG "Error: Invalid p-values in $user_filename detected.\n P-values should be in between 0 to 1";
			print WRITELOG "Please check the column 2 (p-value) for @each_line\n";
			print "\nPlease check column 2 (p-value) for @each_line\n";
            die "Error: Invalid p-values in $user_filename detected.\n";
        }
    }
    else {
        print WRITELOG "Error: Invalid p-values in $user_filename detected.\n";
		print WRITELOG "Please check the column 2 (p-value) for @each_line\n";
		print "\nPlease check column 2 (p-value) for @each_line\n";
        die "Error: Invalid p-values in $user_filename detected.\n";
    }
}
close PVALS;
print "done.\n";

}

sub testingglist
{

# Check the glist file
print "\nValidating glist file...";
open(HGLIST, "tempglist"); # Location of glist file file
@hglist = <HGLIST>;

foreach ( @hglist ) {
    chomp($_);
    @each_gene = split(/\s/,$_);
	$no_of_columns = scalar @each_gene;
	if ($no_of_columns ne 4 ){
		die "Error: glist file must contain 4 columns\n";
	}

	$given_chromosome = @each_gene[0];
	$given_start = @each_gene[1];
	$given_stop = @each_gene[2];
	$given_symbol = @each_gene[3];
	$given_genesize = $given_stop - $given_start;

# Make sure column 1 is a valid chromosome column

    if ( looks_like_number( $given_chromosome ) ) {
        $checked_genenumber = 0;
		if ( $given_chromosome =~ /^\d+$/ ) { #chomosome is a whole number
			if ( $given_chromosome < 1 || $given_chromosome > 24 ) {
            	print WRITELOG "Error: Invalid chromosomes in $glist detected.\n Column 1 Chromosomes should be in between 1 to 24.\n";
				print WRITELOG "Please check the column 1 (chromosome number) for @each_gene\n";
				print "\nPlease check column 1 (chromosome number) for @each_gene\n";
            	die "Error: Invalid chromosome in $glist detected. Column 1 should be a chromsome column with number between 1 to 24.\n";
        	}
		}
		else {
        print WRITELOG "Error: Invalid chromsome in $glist detected.\n Chromosome should be a whole number.\n";
		print WRITELOG "Please check the column 1 (chromosome number) for @each_gene\n";
		print "\nPlease check column 1 (chromosome number) for @each_gene\n";
        die "Error: Invalid chromosomes in $glist detected.\n Chromosome should be a whole number.\n";
		}
	}
    else {
        print WRITELOG "Error: Invalid chromsome in $glist detected.\n Please convert X and Y to 23 and 24 respectively if any.\n";
		print WRITELOG "Please check the column 1 (chromosome number) for @each_gene\n";
		print "\nPlease check column 1 (chromosome number) for @each_gene\n";
        die "Error: Invalid chromosomes in $glist detected.\n Please convert X and Y to 23 and 24 respectively if any.\n";
    }

# Make sure column 2 and 3 are valid location columns
    if ( looks_like_number( $given_genesize ) ) {
        $checked_genenumber = 0;
		if ( $given_genesize =~ /^\d+$/ ) { #chomosome is a whole number
			if ( $given_genesize < 1 || $given_genesize > 30000000 ) {
            	print WRITELOG "Error: Invalid location in $glist detected.\n Difference between column 3 (stop) and column 2 (start) should be in between 1 base to 3 Mb.\n";
				print WRITELOG "Please check the column 3 (stop) and column 3 (start) for @each_gene\n";
				print "\nPlease check column 3 (stop) and column 2 (start) for @each_gene\n";
            	die "Error: Invalid location in $glist detected. Difference between column 3 (stop) and column 2 (start) should be in between 1 base to 3 Mb.\n";
        	}
		}
		else {
        print WRITELOG "Error: Invalid location in $glist detected.\n Start and stop values should be a whole number.\n";
		print WRITELOG "Please check the column 2 (start) and column 3 (stop) for @each_gene\n";
		print "\nPlease check column 2 (start) and column 3 (stop) for @each_gene\n";
        die "Error: Invalid location in $glist detected.\n Start and stop values should be a whole number.\n";
		}
	}
    else {
        print WRITELOG "Error: Invalid location in $glist detected.\n Make sure locations (columns 2 and 3) are numeric.\n";
		print WRITELOG "Please check column 2 (start) and column 3 (stop) for @each_gene\n";
		print "\nPlease check column 2 (start) and column 3 (stop) for @each_gene\n";
        die "Error: Invalid location in $glist detected.\n Make sure locations (columns 2 and 3) are numeric.\n";
    }

	#Now make a new genelist with upadated upper and lower limit
	$new_start = $given_start - $lower;
	if($new_start < 0){$new_start = 1;}
	$new_stop = $given_stop + $upper;

	open(PRINTRESULTS, ">>NEWglist");
	print PRINTRESULTS "$given_chromosome $new_start $new_stop $given_symbol\n";

}
close PRINTRESULTS;
close HGLIST;

print "done.\n";


}


###--- End subroutines ---###
