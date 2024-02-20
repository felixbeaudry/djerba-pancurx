#!/usr/bin/perl

#use strict;
use warnings;
use Tabix;
use Data::Dumper;
use JSON;
use Bio::DB::Fasta;
use List::Util qw/all any/;

my $rootPath = shift;
my $donor = shift;
my $sample = shift;
my $normal = shift;
my $seqType = shift;
my $aligner = shift;
my (%data, %vars);

my $outDir = shift;
if (!defined $outDir)
{
	$outDir = "./";
}

my $summaryFile = "$outDir/$sample.summary.csv";
my $variantsFile = "$outDir/$sample.variants.csv";

unless (-e "$rootPath/$donor/$sample/$seqType/$aligner")
{
	die "$rootPath/$donor/$sample/$seqType/$aligner not a valid path!\n";
}

my $hugoSynPath = "/.mounts/labs/PCSI/raw_data/hugo/HUGO_synonyms_171213.txt";

#PCSI_1014_Pa_P_merged_final_snvs.vcf.gz
my $snvPath = "$rootPath/$donor/$sample/$seqType/$aligner/final_strelka2_mutect2/$sample\_merged_final.vcf.gz";
my $indelPath = "$rootPath/$donor/$sample/$seqType/$aligner/final_indel/$sample\_merged_final_indels.vcf.gz";
my $sgvPath = "$rootPath/$donor/$sample/$seqType/$aligner/final_germline_variants/$sample\_germline_final.vcf.gz"; #PCSI_1018_Lv_M_finalannotatedSV.txt
my $svPath = "$rootPath/$donor/$sample/$seqType/$aligner/mergeSV/$sample\_finalannotatedSV.txt";
my $paramPath = "$rootPath/$donor/$sample/$seqType/$aligner/celluloidXY/v0.11.7/solution/parameters_$sample.txt";
my $segPath = "$rootPath/$donor/$sample/$seqType/$aligner/celluloidXY/v0.11.7/solution/segments_$sample.txt.sorted.bed.gz";
my $neoPaths = "$rootPath/$donor/$sample/$seqType/$aligner/netMHC/pan-2.8/polysolver/1.0/";
my $xenomePaths = "$rootPath/$donor/$sample/$seqType/xenome/1.0.1-r/*.log";
my $cosmicSignnlsPath = "$rootPath/$donor/$sample/$seqType/$aligner/cosmicSigNNLS/${sample}_signatures.txt";

#/.mounts/labs/PCSI/analysis/oicr/PCSI0218/PCSI_0218_Pa_P_526/wgs/bwa/0.7.12/tdp_score
my $tdpPath = "$rootPath/$donor/$sample/$seqType/$aligner/tdp_score/$sample.tdp_score.csv";
my $hrdetectPath = "$rootPath/$donor/$sample/$seqType/$aligner/hrdetect/${sample}.hrdetect_score.csv";

my $alexandrovFile = "/.mounts/labs/PCSI/users/rdenroche/phoenixParser/alexandrov_signature.tsv";
#my $rnaSigFile = "/.mounts/labs/PCSI/users/rdenroche/phoenixParser/RNA-subtype.txt";
#my $rnaSigFile = "/.mounts/labs/PCSI/users/rdenroche/phoenixParser/rna_signature.tsv";
my $rnaSigFile = "/.mounts/labs/PCSI/users/azhang/RNA/RNA-subtype.txt";
my $cosmicCensusFile = "/.mounts/labs/PCSI/users/rdenroche/phoenix/phoenixParser/cosmic_census_20161115.tsv";

my $dsbrList = "/.mounts/labs/PCSI/users/rdenroche/lists/dsbr_list.tsv";
my $mmrList = "/.mounts/labs/PCSI/users/rdenroche/lists/mmr_list.tsv";

##AMY externalID
my %external_ID = (
	"PCSI_0657_Pc_X_526" => "SPC123AR",
	"PCSI_0787_Pc_X_526" => "SPC285AR",
);

my %synonym;
parseHugo($hugoSynPath, \%synonym, \%vars);

#my $refSeqFile = "/oicr/data/genomes/homo_sapiens_mc/refSeq/refSeq_genes.bed";
my $refSeqFile = "/.mounts/labs/PCSI/references/GRCh38/rna/refSeq/hg38_refSeq_genes.txt";
my %genes = lookUpGenePositions($refSeqFile, \%vars, \%synonym);

my %tabix;
my $file;

$tabix{snv_annovar} = "null"; #PCSI_1018_Lv_M_merged_final_snvs.hg38_multianno.txt.gz
$file = "$rootPath/$donor/$sample/$seqType/$aligner/final_strelka2_mutect2/annovar/20170716/$sample\_merged_final.hg38_multianno.txt.gz";
if (-e $file)
{
	$tabix{snv_annovar} = Tabix->new(-data => $file, -index => "$file.tbi");
} else {
	die "can't open $file\n";
}

$tabix{indel_annovar} = "null"; #PCSI_1018_Lv_M_merged_final_indels.hg38_multianno.txt.gz
$file = "$rootPath/$donor/$sample/$seqType/$aligner/final_indel/annovar/20170716/$sample\_merged_final_indels.hg38_multianno.txt.gz";
if (-e $file)
{
	$tabix{indel_annovar} = Tabix->new(-data => $file, -index => "$file.tbi");
} else {
	die "can't open $file\n";
}

$tabix{sgv_annovar} = "null"; #PCSI_1018_Lv_M_snvindel_final
$file = "$rootPath/$donor/$sample/$seqType/$aligner/final_germline_variants/annovar/20170716/$sample\_germline_final.hg38_multianno.txt.gz";
if (-e $file)
{
	$tabix{sgv_annovar} = Tabix->new(-data => $file, -index => "$file.tbi");
} else {
	die "can't open $file\n";
}

#$tabix{cadd_snv} = "null";
#$file = "/.mounts/labs/PCSI/raw_data/cadd/whole_genome_SNVs_inclAnno.tsv.gz";
#if (-e $file)
#{
#	$tabix{cadd_snv} = Tabix->new(-data => $file, -index => "$file.tbi");
#}

#$tabix{cadd_indel} = "null";
#$file = "/.mounts/labs/PCSI/raw_data/cadd/InDels_inclAnno.tsv.gz";
#if (-e $file)
#{
#	$tabix{cadd_indel} = Tabix->new(-data => $file, -index => "$file.tbi");
#}


my $outLine;

$data{donor} = $donor;
$data{tumour} = $sample;
$data{normal} = $normal;
$data{seq_type} = $seqType;

my %coveragePaths = ("tumour", "$rootPath/$donor/$data{tumour}/$seqType/$aligner/coverage/$data{tumour}\_coverage_collapsed.metrics", "normal", "$rootPath/$donor/$data{normal}/$seqType/$aligner/coverage/$data{normal}\_coverage_collapsed.metrics");
#my %coveragePaths = ("tumour", "$rootPath/$donor/$data{tumour}/$seqType/$aligner/coverage/$data{tumour}\_coverage_collapsed.metrics", "normal", "$rootPath/PANX1219/$data{normal}/$seqType/$aligner/coverage/$data{normal}\_coverage_collapsed.metrics");
my $sexPath = "$rootPath/$donor/$data{tumour}/$seqType/$aligner/gendertype/final/$donor.final_gender.txt";

parseCosmic($cosmicCensusFile, \%vars, \%synonym);
##parseNeo($neoPaths,\%data,\%vars);
parseSex($sexPath, \%data);
parseCoverage(\%coveragePaths, \%data);
parseLIMS("http://pinery-prod.hpc.oicr.on.ca:8080/pinery-ws-miso", \%data);
print "beginning germline\n";
parseSGV($sgvPath,\%data,\%vars,\%genes, \%tabix, \%synonym); ##TODO annovar seems working, double check, then common rare?
print "beginning SNV\n";
parseSSM($snvPath,\%data,\%vars,\%genes, \%tabix, \%synonym); ##TODO just check is ok
parseSSM($indelPath,\%data,\%vars,\%genes, \%tabix, \%synonym); ##TODO, check that indels are ok for ins and del, no errors due to bp +/-1
parseSV($svPath,\%data,\%vars,\%genes, \%synonym); ##TODO is this really working?!
print "parseSV done\n";
parseParams($paramPath,\%data);
parseSegs($segPath,\%data,\%vars,\%genes, \%synonym);
parseAlexandrov($alexandrovFile, \%data);
parseRNAsigs($rnaSigFile, \%data);
print "RNA done\n";
parseCosmicSigNNLS($cosmicSignnlsPath,\%data);
print "Cosmic done\n";
parseTDP($tdpPath, \%data);
print "TDP done\n";
if (($sample =~ /_X/) or ($sample =~ /Pa_O/))
{
#	parseXenome($xenomePaths,\%data);
}


#my $genotypePath = "$rootPath/$donor/$data{normal}/$seqType/$aligner/sampleIdentity/genotype/*.genotype $rootPath/$donor/$data{tumour}/$seqType/$aligner/sampleIdentity/genotype/*.genotype";
#parseGenotype($genotypePath, \%data);

# hallmark cutoffs
$data{dsbr_snv_load_cut} = 12000;
$data{dsbr_ct_ratio_cut} = 0.3;
$data{dsbr_del4_load_cut} = 200;
$data{dsbr_del4_ratio_cut} = 0.4;
$data{dsbr_delsv_load_cut} = 50;
$data{dsbr_delsv_ratio_cut} = 0.45;
$data{dsbr_dup_load_cut} = 75;
$data{dsbr_sv_load_cut} = 200;

$data{mmr_snv_load_cut} = 35000;
$data{mmr_indel_load_cut} = 3000;

my @dsbrGenes = qw/BRCA1 BRCA2 PALB2 RAD51 RAD51B RAD51C RAD51D XRCC2 XRCC3 DMC1 BRIP1/;
my @mmrGenes = qw/MLH1 MSH2 MSH6 PMS2 EPCAM/;

doHallmarkScoring(\%data, \%vars, \@dsbrGenes, \@mmrGenes); ##depends on indel final check
print "Hallmark done\n";
#parseHRDetect($hrdetectPath, \%data);
print "parser done\n";
parseManualDeficiency(\%data,$dsbrList,"dsbr");
parseManualDeficiency(\%data,$mmrList,"mmr");

# print one line values

open (FILE, ">$summaryFile") or die "Couldn't open $summaryFile\n";

my @headers = qw/donor tumour normal external_id seq_type tumour_coverage normal_coverage mouse_content 
	tumour_reads_per_sp normal_reads_per_sp tumour_error_rate normal_error_rate tumour_soft_clip normal_soft_clip inferred_sex genotype_concordance
	snv_count del_count ins_count indel_count deleterious_count nonsyn_count stopgain_count stoploss_count splice_count noncoding_count
	frameshift_count nonframeshift_count del_frameshift_count del_nonframeshift_count ins_frameshift_count ins_nonframeshift_count
	sv_count sv_del_count sv_dup_count sv_inv_count sv_tra_count sv_del_bp_gene_count sv_dup_bp_gene_count sv_inv_bp_gene_count sv_tra_bp_gene_count
	waddell neo_antigens neo_antigens_weak neo_antigens_strong hla_types cellularity ploidy alexandrov_class alexandrov_etiology moffitt collisson bailey
	dsbr_deficient dsbr_first_hit dsbr_second_hit mmr_deficient mmr_first_hit mmr_second_hit
	csnnls_sig1 csnnls_sig2 csnnls_sig3 csnnls_sig5 csnnls_sig6 csnnls_sig8 csnnls_sig13 csnnls_sig17 csnnls_sig18 csnnls_sig20 csnnls_sig26 csnnls_residuals csnnls_n.mutations
	germline_snv_count germline_indel_count germline_titv_ratio germline_missense_count germline_nonsense_count
	snv_ca snv_cg snv_ct snv_ta snv_tc snv_tg del_1_count del_4_count ins_1_count ins_4_count
	DEL:100 DEL:1k DEL:10k DEL:100k DEL:1m DEL:10m DEL:>10m DUP:100 DUP:1k DUP:10k DUP:100k DUP:1m DUP:10m DUP:>10m INV:100 INV:1k INV:10k INV:100k INV:1m INV:10m INV:>10m
	snv_ca_aa snv_ca_ac snv_ca_ag snv_ca_at snv_ca_ca snv_ca_cc snv_ca_cg snv_ca_ct snv_ca_ga snv_ca_gc snv_ca_gg snv_ca_gt snv_ca_ta snv_ca_tc snv_ca_tg snv_ca_tt
	snv_cg_aa snv_cg_ac snv_cg_ag snv_cg_at snv_cg_ca snv_cg_cc snv_cg_cg snv_cg_ct snv_cg_ga snv_cg_gc snv_cg_gg snv_cg_gt snv_cg_ta snv_cg_tc snv_cg_tg snv_cg_tt
	snv_ct_aa snv_ct_ac snv_ct_ag snv_ct_at snv_ct_ca snv_ct_cc snv_ct_cg snv_ct_ct snv_ct_ga snv_ct_gc snv_ct_gg snv_ct_gt snv_ct_ta snv_ct_tc snv_ct_tg snv_ct_tt
	snv_ta_aa snv_ta_ac snv_ta_ag snv_ta_at snv_ta_ca snv_ta_cc snv_ta_cg snv_ta_ct snv_ta_ga snv_ta_gc snv_ta_gg snv_ta_gt snv_ta_ta snv_ta_tc snv_ta_tg snv_ta_tt
	snv_tc_aa snv_tc_ac snv_tc_ag snv_tc_at snv_tc_ca snv_tc_cc snv_tc_cg snv_tc_ct snv_tc_ga snv_tc_gc snv_tc_gg snv_tc_gt snv_tc_ta snv_tc_tc snv_tc_tg snv_tc_tt
	snv_tg_aa snv_tg_ac snv_tg_ag snv_tg_at snv_tg_ca snv_tg_cc snv_tg_cg snv_tg_ct snv_tg_ga snv_tg_gc snv_tg_gg snv_tg_gt snv_tg_ta snv_tg_tc snv_tg_tg snv_tg_tt
	dsbr_snv_load dsbr_ct_ratio dsbr_del4_load dsbr_del4_ratio dsbr_delsv_load dsbr_delsv_ratio dsbr_dup_load dsbr_sv_load dsbr_first_genes dsbr_second_genes dsbr_score
	dsbr_snv_load_cut dsbr_ct_ratio_cut dsbr_del4_load_cut dsbr_del4_ratio_cut dsbr_delsv_load_cut dsbr_delsv_ratio_cut dsbr_dup_load_cut dsbr_sv_load_cut
	mmr_snv_load mmr_indel_load mmr_first_genes mmr_second_genes mmr_score mmr_snv_load_cut mmr_indel_load_cut
	tdp_score hrdetect_score hrdetect_myriad_score hrdetect_deletion_microhomology_proportion hrdetect_snv_signature_3 hrdetect_snv_signature_8 hrdetect_rearrangement_signature_3 hrdetect_rearrangement_signature_5 hrdetect_loh hrdetect_tai hrdetect_lst
	snv_file indel_file sgv_file sv_file param_file seg_file neo_files xenome_logs rna_sig_file alexandrov_file json_files genotype_files/;

$outLine = "";
for my $type (@headers)
{
	$outLine .= ",$type";
}
$outLine =~ s/^,//;
print FILE "$outLine\n";


$outLine = "";
for my $type (@headers)
{
	if (exists $data{$type})
	{
		$outLine .= ",$data{$type}";
	}
	else
	{
		$outLine .= ",NA";
	}
}
$outLine =~ s/^,//;
print FILE "$outLine\n";

close FILE;


# print genes and variants

open (FILE, ">$variantsFile") or die "Couldn't open $variantsFile\n";

my @everyLine = qw/donor tumour normal external_id seq_type tumour_coverage normal_coverage cellularity ploidy/;
@headers = qw/gene full_name gene_chr copy_number ab_counts mutation_class mutation_type position fusion_genes base_change tumour_freq tumour_depth normal_freq normal_depth nuc_context aa_context dbsnp cosmic neoantigen cadd_phred rarity 1000G_all ExAC_all ESP6500siv2_all gnomad_all clinvar gene_position maf_mean maf_p cosmic_census_flag cosmic_census_data synonyms entrez_id ensembl_id refseq_id pubmed_ids gene_fam_id gene_fam_name/;

$outLine = "";
for my $type (@everyLine)
{
	$outLine .= ",$type";
}
for my $type (@headers)
{
	$outLine .= ",$type";
}
$outLine =~ s/^,//;
print FILE "$outLine\n";


for my $gene (sort comparePos keys %vars)
{
	if ((exists $vars{$gene}{variants}) and (scalar (keys %{ $vars{$gene}{variants} } > 0)))
	{
		if ($gene eq "PALB2")
		{
			warn "$gene\t$vars{$gene}{copy_number}\n";
		}
		for my $var (sort keys %{ $vars{$gene}{variants} })
		{
			if ($gene eq "PALB2")
			{
				warn " printing $var\n";
			}
			$outLine = "";
			for my $type (@everyLine)
			{
				if (exists $data{$type})
				{
					$outLine .= ",$data{$type}";
				}
				else
				{
					$outLine .= ",NA";
				}
			}
			for my $type (@headers)
			{
				if ((exists $vars{$gene}{variants}{$var}{$type}) and (defined $vars{$gene}{variants}{$var}{$type}))
				{
					$outLine .= ",$vars{$gene}{variants}{$var}{$type}";
				}
				elsif (exists $vars{$gene}{$type})
				{
					unless (defined $vars{$gene}{$type}) {
						if ($type =~ /(ensembl|entrez)_id/) {
							$vars{$gene}{$type} = 'NA';
						}
					}
					$outLine .= ",$vars{$gene}{$type}";
				}
				else
				{
					$outLine .= ",NA";
				}
			}
			$outLine =~ s/^,//;
			print FILE "$outLine\n";
		}
	}
	else
	{
		$outLine = "";
		for my $type (@everyLine)
		{
			if (exists $data{$type})
			{
				$outLine .= ",$data{$type}";
			}
			else
			{
				$outLine .= ",NA";
			}
		}
		for my $type (@headers)
		{
			if (exists $vars{$gene}{$type})
			{
				unless (defined $vars{$gene}{$type}) {
					if ($type =~ /(entrez|ensembl)_id/) {
					   $vars{$gene}{$type} = 'NA';
					}
				}
				$outLine .= ",$vars{$gene}{$type}";
			}
			else
			{
				$outLine .= ",NA";
			}
		}
		$outLine =~ s/^,//;
		print FILE "$outLine\n";
	}
}

close FILE;



sub parseHugo
{
	my $file = shift;
	my $synonym = shift;
	my $vars = shift;

	my $l;
	my @header;
	my @f;
	my %row;

	my $chr;

	my $g;
	my $tempg;

	open (FILE, $file) or die "Couldn't open $file\n";
	while ($l = <FILE>)
	{
		chomp $l;
		if ($l =~ /^HGNC ID/)
		{
			@header = split(/\t/, $l);
		}
		else
		{
			%row = ();
			@f = split(/\t/, $l);
			for (my $i = 0; $i < scalar(@f); $i++)
			{
				$row{$header[$i]} = $f[$i];
			}

			unless ($row{Status} =~ /Withdrawn/)
			{
				$g = $row{"Approved Symbol"};

				$vars->{$g}{gene} = $g;
				
				$vars->{$g}{full_name} = $row{"Approved Name"};
				$vars->{$g}{full_name} =~ s/,/;/g;

				$vars->{$g}{gene_chr} = $row{"Chromosome"};
				$vars->{$g}{gene_chr} =~ s/,/;/g;


				$vars->{$g}{gene_chr_for_syn} = $chr;


				$vars->{$g}{synonyms} = "$row{'Previous Symbols'}, $row{'Synonyms'}";
				$vars->{$g}{synonyms} =~ s/,/;/g;
				$vars->{$g}{synonyms} =~ s/^;//;
				$vars->{$g}{synonyms} =~ s/;$//;

				$vars->{$g}{entrez_id} = $row{"Entrez Gene ID"};
				$vars->{$g}{ensembl_id} = $row{"Ensembl Gene ID"};
				if (exists $row{"RefSeq IDs"})
				{
					$vars->{$g}{refseq_id} = $row{"RefSeq IDs"};
					$vars->{$g}{refseq_id} =~ s/,/;/g;
				}

				if (exists $row{"Pubmed IDs"})
				{
					$vars->{$g}{pubmed_ids} = $row{"Pubmed IDs"};
					$vars->{$g}{pubmed_ids} =~ s/,/;/g;
				}

				if (exists $row{"Gene Family Name"})
				{
					$vars->{$g}{gene_fam_id} = $row{"Gene Family ID"};
					$vars->{$g}{gene_fam_name} = $row{"Gene Family Name"};
					$vars->{$g}{gene_fam_name} =~ s/,/;/g;
				}

				for my $c (split(/ and /, $row{"Chromosome"}))		# to handle "Xp22.33 and Yp11.32"
				{
					if ($c =~ /^(.*?)[pq]/)
					{
						$chr = "chr$1";
					}
					else
					{
						$chr = "chr" . $c;
					}

					for my $sym (split(/;/, $vars->{$g}{synonyms}))
					{
						if (exists $synonym->{$sym}{$chr})
						{
							$synonym->{$sym}{$chr} .= ";$g";
						}
						else
						{
							$synonym->{$sym}{$chr} = $g;
						}
					}
				}

			}



		}
	}
	close FILE;

	for $g (keys %{ $vars })
	{
		$synonym->{$g}{$vars->{$g}{gene_chr_for_syn}} = $g;
	}
}








sub parseSSM
{
	my $file = shift;
	my $data = shift;
	my $vars = shift;
	my $genes = shift;
	my $tabix = shift;
	my $synonym = shift;

	my $l;
	my ($donor,$tumour,$normal,$tCov,$nCov,$gene,$consequence,$change);
	my ($chr,$pos,$id,$ref,$alt,$alts,$qual,$filter,$info,$format,$nGT,$tGT);

	my $snvCount = 0;
	my $insCount = 0;
	my $delCount = 0;

	my $deleteriousCount = 0;
	my %anno;

#	my %fastaHandles = ("path" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/fasta/");
	my %fastaHandles = ("path" => "/.mounts/labs/PCSI/references/GRCh38/dna/fasta/");
	my %referenceHash = ("chunkSize" => 10000);

	my %comp = (
	    "ac" => "tg",
	    "ag" => "tc",
	    "at" => "ta",
	    "ga" => "ct",
	    "gc" => "cg",
	    "gt" => "ca",
	);

	my %changeCount = (
		"tg" => 0,
		"tc" => 0,
		"ta" => 0,
		"ct" => 0,
		"cg" => 0,
		"ca" => 0,
	);
	
	my %contextComp = (
	    "aa" => "tt",
	    "ac" => "gt",
	    "ag" => "ct",
	    "at" => "at",
	    "ca" => "tg",
	    "cc" => "gg",
	    "cg" => "cg",
	    "ct" => "ag",
	    "ga" => "tc",
	    "gc" => "gc",
	    "gg" => "cc",
	    "gt" => "ac",
	    "ta" => "ta",
	    "tc" => "ga",
	    "tg" => "ca",
	    "tt" => "aa",
	);

	my $context;
	my %contextCount;
	for $change (qw/ca cg ct ta tc tg/)
	{
		for $context (qw/aa ac ag at ca cc cg ct ga gc gg gt ta tc tg tt/)
		{
			$contextCount{"${change}_$context"} = 0;
		}
	}


	my $del1Count = 0;
	my $del4Count = 0;
	my $ins1Count = 0;
	my $ins4Count = 0;

	my $nonsynCount = 0;
	my $stopgainCount = 0;
	my $stoplossCount = 0;
	my $spliceCount = 0;
	my $noncodingCount = 0;
	my $frameshiftCount = 0;
	my $nonframeshiftCount = 0;

	my $delFrameshiftCount = 0;
	my $delNonframeshiftCount = 0;
	my $insFrameshiftCount = 0;
	my $insNonframeshiftCount = 0;

	my $cosmicFlag;

	my %baseCounts;
	my ($varName,$chrPos,$mutClass,$mutType,$nFreq,$tFreq,$nDepth,$tDepth,$baseChange,$nucContext,$aaContext,$dbSNP,$COSMIC,$rarity,$thouG,$ExAC,$ESP6500,$caddPhred,$clinvar,$valuesExist,$popFreq,$annovar);

	my $isDeletion;
	my $isInsertion;

	my %infoHash;

	if (-e $file)
	{
		open (FILE, "zcat $file | ") or die "Couldn't open $file\n";
	
		while ($l = <FILE>)
		{
			$isDeletion = 0;
			$isInsertion = 0;

			chomp $l;
			if ($l =~ /^##normal_sample=(([A-Z0-9]+)_([0-9]+)_.*)$/) 
			{
				$normal = $1;
				$donor = $2.$3;
			}
			elsif ($l =~ /^##tumor_sample=(.*)$/) 
			{
				$tumour = $1;
			}
			elsif ($l =~ /##DCC=<analyzed_seq_coverage=(.*?)>/)
			{
				$tCov = sprintf("%0.2f", $1);
			}
			elsif ($l =~ /##DCC=<matched_seq_coverage=(.*?)>/)
			{
				$nCov = sprintf("%0.2f", $1);
			}
			elsif ($l =~ /^#/)
			{
				next;
			}
			else
			{
				($chr,$pos,$id,$ref,$alts,$qual,$filter,$info,$format,$nGT,$tGT) = split(/\t/, $l);
	
				for $alt (split(/,/, $alts))
				{
					$gene = "";
					$consequence = "";

					$varName = "";
					$chrPos = "$chr:$pos";
					$mutClass = "";
					$mutType = "";
					$nFreq = "NA";
					$nDepth = "NA";
					$tFreq = "NA";
					$tDepth = "NA";
					$baseChange = ".";
					$nucContext = ".";
					$aaContext = ".";
					$dbSNP = ".";
					$COSMIC = ".";
					$rarity = "NA";
					$thouG = "NA";
					$ExAC = "NA";
					$ESP6500 = "NA";
					$caddPhred = "NA";
					$clinvar = "NA";


					if (length($ref) == length($alt))
					{

						$snvCount++;
						$mutClass = "somatic snv";
						$annovar = "snv_annovar";

						$baseChange = "$ref>$alt";
	
						$change = lc("$ref$alt");
						unless (length($ref) > 1) { #skip context for multint substitutions
						   $context = lc(getBase($chr, $pos - 1, \%referenceHash, \%fastaHandles) . getBase($chr, $pos + 1, \%referenceHash, \%fastaHandles));
						   if (exists $comp{$change}) {
							$change = $comp{$change};
							unless ($context =~ /n/)
							{
								$context = $contextComp{$context};
							}
						   }
						   $changeCount{$change}++;
						   unless ($context =~ /n/) {
							$contextCount{"${change}_$context"}++;
						   }

						   freeMemChunks($chr, $pos, \%referenceHash);

						} #end of context for single nt substitutions
						if ($nGT =~ /^.*:.*,.*:.*:(.*):([0-9]+,[0-9]+):([0-9]+,[0-9]+):.*/) { #GT:AD:AF:DP:F1R2:F2R1:optional(PGT:PID:PS:)SB
							$nDepth = $1;
							my ($ref1, $alt1) = split /,/, $2;
							my ($ref2, $alt2) = split /,/, $3;
							$nFreq = '.';
							if ($nDepth > 0) {
								if (($alt1 + $alt2) == 0) {
									$nFreq = 0
								} else {
									$nFreq = ($alt1 + $alt2) / ($ref1 + $ref2 + $alt1 + $alt2);
								}
							}
						}
						if ($tGT =~ /^.*:.*,.*:.*:(.*):([0-9]+,[0-9]+):([0-9]+,[0-9]+):.*/) {
							$tDepth = $1;
							my ($ref1, $alt1) = split /,/, $2;
							my ($ref2, $alt2) = split /,/, $3;
							$tFreq = '.';
							if ($tDepth > 0) {
#								0/1:7,14:0.653:21:7,11:0,3:3,4,6,8
								$tFreq = ($alt1 + $alt2) / $tDepth;
							}
						}

						#DP:FDP:SDP:SUBDP:AU:CU:GU:TU
						if ($nGT =~ /^(.*?):.*:.*:.*:(.*?),.*?:(.*?),.*?:(.*?),.*?:(.*?),.*?$/)
						{
							$nDepth = $1;
							$baseCounts{A} = $2;
							$baseCounts{C} = $3;
							$baseCounts{G} = $4;
							$baseCounts{T} = $5;
							$nFreq = ".";
							if ($nDepth > 0)
							{
								$nFreq = $baseCounts{$alt} / $nDepth;
							}
						}
						if ($tGT =~ /^(.*?):.*:.*:.*:(.*?),.*?:(.*?),.*?:(.*?),.*?:(.*?),.*?$/)
						{
							$tDepth = $1;
							$baseCounts{A} = $2;
							$baseCounts{C} = $3;
							$baseCounts{G} = $4;
							$baseCounts{T} = $5;
							$tFreq = ".";
							if ($tDepth > 0)
							{
								$tFreq = $baseCounts{$alt} / $tDepth;
							}
						}

					}
					elsif (length($ref) > length($alt))
					{
						$delCount++;
						$isDeletion = 1;
						$annovar = "indel_annovar";

						$mutClass = "somatic indel";
						$baseChange = "$ref>$alt";

						if (length($ref) > 4)
						{
							$del4Count++;
						}
						else
						{
							$del1Count++;
						}

						#GT:AD:AF:DP:F1R2:F2R1:optional
						if ($nGT =~ /^.*:.*,.*:.*:(.*):([0-9]+,[0-9]+):([0-9]+,[0-9]+):.*/) { 
                                                        $nDepth = $1;
                                                        my ($ref1, $alt1) = split /,/, $2;
                                                        my ($ref2, $alt2) = split /,/, $3;
                                                        $nFreq = '.';
                                                        if ($nDepth > 0) {
                                                                $nFreq = ($alt1 + $alt2) / ($ref1 + $ref2 + $alt1 + $alt2);
                                                        }
                                                }
                                                if ($tGT =~ /^.*:.*,.*:.*:(.*):([0-9]+,[0-9]+):([0-9]+,[0-9]+):.*/) {
                                                        $tDepth = $1;
                                                        my ($ref1, $alt1) = split /,/, $2;
                                                        my ($ref2, $alt2) = split /,/, $3;
                                                        $tFreq = '.';
                                                        if ($tDepth > 0) {
                                                                $tFreq = ($alt1 + $alt2) / ($ref1 + $ref2 + $alt1 + $alt2);
                                                        }
						}
						#DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50
						if ($format =~ /DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50/ and $nGT =~ /^(.*?):.*?:.*?,.*?:(.*?),.*/)
						{
							$nDepth = $1;
							$nFreq = ".";
							if ($nDepth > 0)
							{
								$nFreq = $2 / $nDepth;
							}
						}
						if ($format =~ /DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50/ and $tGT =~ /^(.*?):.*?:.*?,.*?:(.*?),.*/)
						{
							$tDepth = $1;
							$tFreq = ".";
							if ($tDepth > 0)
							{
								$tFreq = $2 / $tDepth;
							}
						}

					}
					else
					{
						$insCount++;
						$isInsertion = 1;
						$annovar = "indel_annovar";

						$mutClass = "somatic indel";
						$baseChange = "$ref>$alt";

						if (length($alt) > 4)
						{
							$ins4Count++;
						}
						else
						{
							$ins1Count++;
						}
						#DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50
						if ($nGT =~ /^.*:.*,.*:.*:(.*):([0-9]+,[0-9]+):([0-9]+,[0-9]+):.*/) { 
                                                        $nDepth = $1;
                                                        my ($ref1, $alt1) = split /,/, $2;
                                                        my ($ref2, $alt2) = split /,/, $3;
                                                        $nFreq = '.';
                                                        if ($nDepth > 0) {
                                                                $nFreq = ($alt1 + $alt2) / ($ref1 + $ref2 + $alt1 + $alt2);
                                                        }
                                                }
                                                if ($tGT =~ /^.*:.*,.*:.*:(.*):([0-9]+,[0-9]+):([0-9]+,[0-9]+):.*/) {
                                                        $tDepth = $1;
                                                        my ($ref1, $alt1) = split /,/, $2;
                                                        my ($ref2, $alt2) = split /,/, $3;
                                                        $tFreq = '.';
                                                        if ($tDepth > 0) {
                                                                $tFreq = ($alt1 + $alt2) / ($ref1 + $ref2 + $alt1 + $alt2);
                                                        }
						}
						if ($format =~ /DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50/ and $nGT =~ /^(.*?):.*?:.*?,.*?:(.*?),.*/)
						{
							$nDepth = $1;
							$nFreq = ".";
							if ($nDepth > 0)
							{
								$nFreq = $2 / $nDepth;
							}
						}
						if ($format =~ /DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50/ and $tGT =~ /^(.*?):.*?:.*?,.*?:(.*?),.*/)
						{
							$tDepth = $1;
							$tFreq = ".";
							if ($nDepth > 0)
							{
								$tFreq = $2 / $tDepth;
							}
						}
					}
	
					unless ($info =~ /;$/) {
						$info = "$info;";
					}
					%infoHash = ();
					%infoHash = parseANNOVAR($info, $baseChange);

					if ($info =~ /COSMIC=(.*?);/)
					{
						$COSMIC = $1;
					}



					for my $t (sort keys %infoHash)
					{
						for my $g (sort keys %{ $infoHash{$t} })
						{
							$gene = $g;
							$consequence = $infoHash{$t}{$g}{consequence};
							$nucContext = $infoHash{$t}{$g}{nuc};
							$aaContext = $infoHash{$t}{$g}{aa};

							if ($consequence eq "splicing")
							{
								$mutType = "splicing";
								$spliceCount++;
							}
							elsif (($consequence =~ /^nonframeshift/) or ($consequence =~ /^nonsynonymous/) or ($consequence =~ /^stoploss/))
							{
								$deleteriousCount++;
								if ($consequence =~ /^nonframeshift/)
								{
									$nonframeshiftCount++;
									$mutType = "nonframeshift";
									if ($isDeletion == 1)
									{
										$delNonframeshiftCount++;
									}
									else
									{
										$insNonframeshiftCount++;
									}
								}
								elsif ($consequence =~ /^nonsynonymous/)
								{
									$nonsynCount++;
									$mutType = "nonsynonymous";
								}
								elsif ($consequence =~ /^stoploss/)
								{
									$stoplossCount++;
									$mutType = "stoploss";
								}
							}
							elsif (($consequence =~ /^stopgain/) or ($consequence =~ /^frameshift/))
							{
								$deleteriousCount++;
								if ($consequence =~ /^stopgain/)
								{
									$stopgainCount++;
									$mutType = "stopgain";
								}
								elsif ($consequence =~ /^frameshift/)
								{
									$frameshiftCount++;
									$mutType = "frameshift";
									if ($isDeletion == 1)
									{
										$delFrameshiftCount++;
									}
									else
									{
										$insFrameshiftCount++;
									}
								}
							}
							elsif ($consequence =~ /^unknown/)
							{
								$mutType = "unknown";		# shrug...
							}
		
		
		
							unless ($mutType eq "")
							{
								$gene = findBestGeneSymbol($gene, $vars, $synonym, $chr, "ANNOVAR (SSM)", "warn");
								%anno = ();
								doANNOVARlookup($chr,$pos,$ref,$alt,$tabix->{$annovar},\%anno,$info);
								$popFreq = -1;
								$valuesExist = 0;
								if (exists $anno{"avsnp150"}) {
									$id = $anno{"avsnp150"};
								}
								for my $value (qw/1000g2015aug_all ExAC_ALL AF/)
								{
									if (exists $anno{$value})
									{
										$valuesExist = 1;
										unless ($anno{$value} eq "NA")
										{
											if ($anno{$value} > $popFreq)
											{
												$popFreq = $anno{$value};
											}
										}
								 	} else {
                                                                                warn "did not find $value for $l\n";
                                                                                warn Dumper %anno;
                                                                                die;
                                                                        }
								}
								$cosmicFlag = "NA";
								if (exists $vars->{$gene}{cosmic_census_types}{$mutType})
								{
									$cosmicFlag = "cosmic_mutation";
								}
								
								if ($popFreq > 0.01)
								{
									$rarity = "common";
								}
								elsif ($popFreq > 0)
								{
									$rarity = "rare"
								}
								elsif ($valuesExist == 1)
								{
									$rarity = "novel";
								}
		
								$varName = "SSM $mutType $chrPos $baseChange";
								$vars->{$gene}{gene} = $gene;
								$vars->{$gene}{variants}{$varName}{mutation_class} = $mutClass;
								$vars->{$gene}{variants}{$varName}{mutation_type} = $mutType;
								$vars->{$gene}{variants}{$varName}{position} = $chrPos;
								$vars->{$gene}{variants}{$varName}{base_change} = $baseChange;
		
								$vars->{$gene}{variants}{$varName}{tumour_freq} = sprintf("%.3f",$tFreq);
								$vars->{$gene}{variants}{$varName}{tumour_depth} = $tDepth;
								$vars->{$gene}{variants}{$varName}{normal_freq} = sprintf("%.3f",$nFreq);
								$vars->{$gene}{variants}{$varName}{normal_depth} = $nDepth;
		
								$vars->{$gene}{variants}{$varName}{nuc_context} = $nucContext;
								$vars->{$gene}{variants}{$varName}{aa_context} = $aaContext;
								$vars->{$gene}{variants}{$varName}{dbsnp} = $id;
								$vars->{$gene}{variants}{$varName}{cosmic} = $COSMIC;
								$vars->{$gene}{variants}{$varName}{rarity} = $rarity;
								$vars->{$gene}{variants}{$varName}{"1000G_all"} = $anno{'1000G_ALL'};
								$vars->{$gene}{variants}{$varName}{ExAC_all} = $anno{ExAC_ALL};
								$vars->{$gene}{variants}{$varName}{gnomad_all} = $anno{AF};
#								$vars->{$gene}{variants}{$varName}{cadd_phred} = $anno{CADD_phred};
								$vars->{$gene}{variants}{$varName}{clinvar} = $anno{CLNSIG};
								$vars->{$gene}{variants}{$varName}{cosmic_census_flag} = $cosmicFlag;
							}
						}
					}
	
					$mutType = "";
					#funseq
					if ($info =~ /GENE=(.*?)\(Promoter\)/)
					{
						warn "found a promoter, $info\n";
						die;
						$gene = $1;
						$noncodingCount++;
						$mutType = "altered promoter";
					}
					if ($info =~ /NCENC=(.*)/)
					{
						for my $g (split(/,/, $1))
						{
							if ($g =~ /(.*?)\((.*)\|.*\)/)
							{
								if ($1 eq "Enhancer")
								{
									$gene = $2;
									$noncodingCount++;
									$mutType = "altered enhancer";
								}
							}
						}
					}

					if (($mutType ne "") and ($gene ne "chmm/segway") and ($gene ne "drm") and (exists $vars->{$gene}))
					{
						%anno = ();
						doANNOVARlookup($chr,$pos,$ref,$alt,$tabix->{$annovar},\%anno, $info);
						$varName = "SSM $mutType $chrPos $baseChange";

						$popFreq = -1;
						$valuesExist = 0;
						if (exists $anno{avsnp150}) {
							$id = $anno{avsnp150};
						}
						for my $value (qw/1000g2015aug_all ExAC_ALL AF/)
						{
							if (exists $anno{$value})
							{
								$valuesExist = 1;
								unless ($anno{$value} eq "NA")
								{
									if ($anno{$value} > $popFreq)
									{
										$popFreq = $anno{$value};
									}
								}
							}
						}
						
						if ($popFreq > 0.01)
						{
							$rarity = "common";
						}
						elsif ($popFreq > 0)
						{
							$rarity = "rare"
						}
						elsif ($valuesExist == 1)
						{
							$rarity = "novel";
						}

						$varName = "SSM $mutType $chrPos $baseChange";
						$vars->{$gene}{variants}{$varName}{mutation_class} = $mutClass;
						$vars->{$gene}{variants}{$varName}{mutation_type} = $mutType;
						$vars->{$gene}{variants}{$varName}{position} = $chrPos;
						$vars->{$gene}{variants}{$varName}{base_change} = $baseChange;
						$vars->{$gene}{variants}{$varName}{tumour_freq} = sprintf("%.3f",$tFreq);
						$vars->{$gene}{variants}{$varName}{tumour_depth} = $tDepth;
						$vars->{$gene}{variants}{$varName}{normal_freq} = sprintf("%.3f",$nFreq);
						$vars->{$gene}{variants}{$varName}{normal_depth} = $nDepth;
						$vars->{$gene}{variants}{$varName}{nuc_context} = $nucContext;
						$vars->{$gene}{variants}{$varName}{aa_context} = $aaContext;
						$vars->{$gene}{variants}{$varName}{dbsnp} = $id;
						$vars->{$gene}{variants}{$varName}{cosmic} = $COSMIC;
						$vars->{$gene}{variants}{$varName}{rarity} = $rarity;
						$vars->{$gene}{variants}{$varName}{"1000G_all"} = $anno{'1000g2015aug_all'};
						$vars->{$gene}{variants}{$varName}{ExAC_all} = $anno{ExAC_ALL};
						$vars->{$gene}{variants}{$varName}{"gnomad_AF"} = $anno{"AF"};
#						$vars->{$gene}{variants}{$varName}{cadd_phred} = $anno{CADD_phred};
						$vars->{$gene}{variants}{$varName}{clinvar} = $anno{CLNSIG};
					}

				}									

			}
		}
		close FILE;
		
		$data->{normal} = $normal;

		if ($file eq $snvPath) {	
		   $data->{snv_file} = $file;
		   $data->{snv_count} = $snvCount;
		   $data->{deleterious_count} = $deleteriousCount;
		   $data->{nonsyn_count} = $nonsynCount;
		   $data->{stopgain_count} = $stopgainCount;
		   $data->{stoploss_count} = $stoplossCount;
		   $data->{splice_count} = $spliceCount;
		   $data->{noncoding_count} = $noncodingCount;
		   for $change (keys %changeCount) {
			$data->{"snv_$change"} = $changeCount{$change};
		   }
		   for $context (keys %contextCount) {
			$data->{"snv_$context"} = $contextCount{$context};
		   }
#		$data->{tumour_coverage} = $tCov; 
#		$data->{normal_coverage} = $nCov;
		} elsif (exists $data->{snv_count}) {	
		   $data->{indel_file} = $file;
		   $data->{deleterious_count} += $deleteriousCount;
		   $data->{ins_count} = $insCount;
		   $data->{del_count} = $delCount;
		   $data->{indel_count} = $insCount + $delCount;
	
		   $data->{del_1_count} = $del1Count;
		   $data->{del_4_count} = $del4Count;
		   $data->{ins_1_count} = $ins1Count;
		   $data->{ins_4_count} = $ins4Count;
	
		   #nonsyn_count,stopgain_count,stoploss_count,frameshift_count,nonframeshift_count
		   $data->{frameshift_count} = $frameshiftCount;
		   $data->{nonframeshift_count} = $nonframeshiftCount;

		   $data->{del_frameshift_count} = $delFrameshiftCount;
		   $data->{del_nonframeshift_count} = $delNonframeshiftCount;
		   $data->{ins_frameshift_count} = $insFrameshiftCount;
		   $data->{ins_nonframeshift_count} = $insNonframeshiftCount;
		} else {
		   $data->{indel_file} = $file;
		   $data->{deleterious_count} = $deleteriousCount;
		}
		
	}
	else
	{
		warn "$file doesn't exist\n";
	}
	#print Dumper $data;
}


sub parseSGV
{
	my $file = shift;
	my $data = shift;
	my $vars = shift;
	my $genes = shift;
	my $tabix = shift;
	my $synonym = shift;

	$data->{sgv_file} = $file;

	# germline_snv_count germline_indel_count germline_common_count germline_rare_count germline_novel_count germline_titv_ratio germline_missense_count germline_nonsense_count
	my $snvCount = 0;
	my $insCount = 0;
	my $delCount = 0;
	my $commonCount = 0;
	my $rareCount = 0;
	my $novelCount = 0;
	my $titvRatio = "NA";
	my $tiCount = 0;
	my $missenseCount = 0;
	my $nonsenseCount = 0;
	my $noncodingCount = 0;
	my $spliceCount = 0;

	my %infoHash;

	my $cosmicFlag;

	my %transition = (
		"AG" => 1,
		"GA" => 1,
		"CT" => 1,
		"TC" => 1,
	);


	my $l;
	my @f;
	my ($chr,$pos,$id,$ref,$alt,$alts,$qual,$filter,$info,$format,$firstGT,$secondGT,$tumIsFirst,$nGT,$tGT);

	my ($gene, $consequence, $isIndel, $popFreq, $valuesExist);

	my %anno;
	my %cadd;

	# qw/gene copy_number ab_counts mutation_class mutation_type position fusion_genes base_change nuc_context aa_context dbsnp cosmic cadd_phred clinvar rarity 1000G_all ExAC_all ESP6500siv2_all gene_position maf_mean maf_p/
	my ($varName, $chrPos, $nFreq, $nDepth, $tFreq, $tDepth, $mutClass, $mutType, $baseChange, $nucContext, $aaContext, $dbSNP, $COSMIC, $rarity);

	if (-e $file)
	{
		open (FILE, "zcat $file | ") or die "Couldn't open $file\n";

		while ($l = <FILE>)
		{
			chomp $l;

			if ($l =~ /^#CHROM/) {
				if ($l =~ /$normal\t$sample/) {
					$firstGT = "normal";
				} else {
					$firstGT = "tumour";
				}

			}
			unless ($l =~ /^#/)
			{
				if ($firstGT eq "normal") {
					($chr,$pos,$id,$ref,$alts,$qual,$filter,$info,$format,$nGT,$tGT) = split(/\t/, $l);
				} else {
					($chr,$pos,$id,$ref,$alts,$qual,$filter,$info,$format,$tGT,$nGT) = split(/\t/, $l);
				}

				for $alt (split(/,/, $alts))
				{
					$gene = "";
					$consequence = "";
					$isIndel = "true";

					$varName = "";
					$chrPos = "$chr:$pos";
					$mutClass = "";
					$mutType = "";
					$nFreq = "NA";
					$nDepth = "NA";
					$tFreq = "NA";
					$tDepth = "NA";
					$baseChange = ".";
					$nucContext = ".";
					$aaContext = ".";
					$dbSNP = ".";
					$COSMIC = ".";
					$rarity = "NA";

					if (length($ref) == length($alt))
					{
						$snvCount++;
						$isIndel = "false";
						$mutClass = "germline snp";

						if (exists $transition{"$ref$alt"})
						{
							$tiCount++;
						}
					}
					elsif (length($ref) > length($alt))
					{
						$delCount++;
						$mutClass = "germline indel";
					}
					else
					{
						$insCount++;
						$mutClass = "germline indel";
					}
					
					#GT:AD:DP:GQ:PL  1/1:4,19:57:33.40:3682,33,0    
					if ($nGT =~ /^.*:.*,(.*):(.*):.*:.*/)
					{
						$nDepth = $2;
						$nFreq = ".";
						if ($nDepth > 0)
						{
							$nFreq = $1 / $nDepth;
						}
					}
					if ($tGT =~ /^.*:.*,(.*):(.*):.*:.*/)
					{
						$tDepth = $2;
						$tFreq = ".";
						if ($tDepth > 0)
						{
							$tFreq = $1 / $tDepth;
						}
					}

					unless ($info =~ /;$/) {
						$info = "$info;";
					}
					%infoHash = parseANNOVAR($info);
					if ($info =~ /COSMIC=(.*?);/)
					{
						$COSMIC = $1;
					}

					for my $t (sort keys %infoHash)
					{
						for my $g (sort keys %{ $infoHash{$t} })
						{
							$gene = $g;
							$consequence = $infoHash{$t}{$g}{consequence};
							$nucContext = $infoHash{$t}{$g}{nuc};
							$aaContext = $infoHash{$t}{$g}{aa};

							if ($consequence eq "splicing")
							{
								$mutType = "splicing";
								$spliceCount++;
							}
							elsif (($consequence =~ /^nonframeshift/) or ($consequence =~ /^nonsynonymous/) or ($consequence =~ /^stoploss/))
							{
								$missenseCount++;
		
								if ($consequence =~ /^nonframeshift/)
								{
									$mutType = "nonframeshift";
								}
								elsif ($consequence =~ /^nonsynonymous/)
								{
									$mutType = "nonsynonymous";
								}
								elsif ($consequence =~ /^stoploss/)
								{
									$mutType = "stoploss";
								}
							}
							elsif (($consequence =~ /^stopgain/) or ($consequence =~ /^frameshift/))
							{
								$nonsenseCount++;
								if ($consequence =~ /^stopgain/)
								{
									$mutType = "stopgain";
								}
								elsif ($consequence =~ /^frameshift/)
								{
									$mutType = "frameshift";
								}
							}
							elsif ($consequence =~ /^splicing/)
							{
								$mutType = "splicing";
							}
		
							unless ($mutType eq "")
							{
								$gene = findBestGeneSymbol($gene, $vars, $synonym, $chr, "ANNOVAR (SGV)", "warn");
								# extended annovar
								%anno = ();
								if ($alt eq '*') {
									next;
								}
								doANNOVARlookup($chr,$pos,$ref,$alt,$tabix->{sgv_annovar},\%anno,$l);
			
								$popFreq = -1;
								$valuesExist = 0;
								for my $value (qw/1000g2015aug_all ExAC_ALL AF/)
								{
									if (exists $anno{$value})
									{
										$valuesExist = 1;
										unless ($anno{$value} eq "NA")
										{
											if ($anno{$value} > $popFreq)
											{
												$popFreq = $anno{$value};
											}
										}
									} else {
										warn "did not find $value for $chr,$pos,$ref,$alt,$l\n";
										if ($l =~ /32530184/) {
										#	warn "$chr\t$pos\t$ref\t$alt\n";
										}
#										warn Dumper %anno;
#										die;
									}
								}
								
								if ($popFreq > 0.01)
								{
									$commonCount++;
									$rarity = "common";
								}
								elsif ($popFreq > 0)
								{
									$rareCount++;
									$rarity = "rare";
								}
								elsif ($valuesExist == 1)
								{
									$novelCount++;
									$rarity = "novel";
								}
		
								$cosmicFlag = "NA";
								if (exists $vars->{$gene}{cosmic_census_types}{$mutType})
								{
									$cosmicFlag = "cosmic_mutation";
								}
		
								$baseChange = "$ref>$alt";

								my @clinvarheaders = grep /^CLNSIG/, sort keys %anno;
								my $clinvar = "";
								for (my $i = 0; $i < @clinvarheaders; $i++) {
									my $header = $clinvarheaders[$i];
									my $clininfo = $anno{$header};
									$clinvar .= "$header=$clininfo;";
								}
		
								$varName = "SGV $mutType $chrPos $baseChange";
								$vars->{$gene}{gene} = $gene;
								$vars->{$gene}{variants}{$varName}{mutation_class} = $mutClass;
								$vars->{$gene}{variants}{$varName}{mutation_type} = $mutType;
								$vars->{$gene}{variants}{$varName}{position} = $chrPos;
								$vars->{$gene}{variants}{$varName}{base_change} = $baseChange;
		
								$vars->{$gene}{variants}{$varName}{tumour_freq} = sprintf("%.3f",$tFreq);
								$vars->{$gene}{variants}{$varName}{tumour_depth} = $tDepth;
								$vars->{$gene}{variants}{$varName}{normal_freq} = sprintf("%.3f",$nFreq);
								$vars->{$gene}{variants}{$varName}{normal_depth} = $nDepth;
		
								$vars->{$gene}{variants}{$varName}{nuc_context} = $nucContext;
								$vars->{$gene}{variants}{$varName}{aa_context} = $aaContext;
								$vars->{$gene}{variants}{$varName}{dbsnp} = $id;
								if ($id eq "rs80357389") {
									$anno{CLINSIG} = "CLINSIG=pathogenic";
								}
								$vars->{$gene}{variants}{$varName}{cosmic} = $COSMIC;
								$vars->{$gene}{variants}{$varName}{rarity} = $rarity;
								$vars->{$gene}{variants}{$varName}{"1000G_all"} = $anno{"1000g2015aug_all"};
								$vars->{$gene}{variants}{$varName}{ExAC_all} = $anno{ExAC_ALL};
								$vars->{$gene}{variants}{$varName}{"gnomad_all"} = $anno{AF};
#								$vars->{$gene}{variants}{$varName}{ESP6500siv2_all} = $anno{ESP6500siv2_ALL};
								$vars->{$gene}{variants}{$varName}{cadd_phred} = $anno{CADD_phred};
#								$vars->{$gene}{variants}{$varName}{clinvar} = $anno{clinvar_20150330};
								$vars->{$gene}{variants}{$varName}{clinvar} = $anno{CLINSIG};
								$vars->{$gene}{variants}{$varName}{cosmic_census_flag} = $cosmicFlag;
							
							}
						}
					}
	
					$mutType = "";
					#funseq
					if ($info =~ /GENE=(.*?)\(Promoter\)/)
					{
						$gene = $1;
						$noncodingCount++;
						$mutType = "altered promoter";
					}
					if ($info =~ /NCENC=(.*)/)
					{
						for my $g (split(/,/, $1))
						{
							if ($g =~ /(.*?)\((.*)\|.*\)/)
							{
								if ($1 eq "Enhancer")
								{
									$gene = $2;
									$noncodingCount++;
									$mutType = "altered enhancer";
								}
							}
						}
					}

					if (($mutType ne "") and ($gene ne "chmm/segway") and ($gene ne "drm") and (exists $vars->{$gene}))
					{
						%anno = ();
						doANNOVARlookup($chr,$pos,$ref,$alt,$tabix->{sgv_annovar},\%anno,$l);
						$varName = "SGV $mutType $chrPos $baseChange";

						$popFreq = -1;
						$valuesExist = 0;
						for my $value (qw/1000g2015aug_all ExAC_ALL AF/)
						{
							if (exists $anno{$value})
							{
								$valuesExist = 1;
								unless ($anno{$value} eq "NA")
								{
									if ($anno{$value} > $popFreq)
									{
										$popFreq = $anno{$value};
									}
								}
							}
						}
						
						if ($popFreq > 0.01)
						{
							$rarity = "common";
						}
						elsif ($popFreq > 0)
						{
							$rarity = "rare"
						}
						elsif ($valuesExist == 1)
						{
							$rarity = "novel";
						}

						$baseChange = "$ref>$alt";

						my @clinvarheaders = grep /^CLNSIG/, sort keys %anno;
						my $clinvar = "";
						for (my $i = 0; $i < @clinvarheaders; $i++) {
							my $header = $clinvarheaders[$i];
							my $clininfo = $anno{$header};   
							$clinvar .= "$header=$clininfo;";
						}
						$varName = "SGV $mutType $chrPos $baseChange";
						$vars->{$gene}{variants}{$varName}{mutation_class} = $mutClass;
						$vars->{$gene}{variants}{$varName}{mutation_type} = $mutType;
						$vars->{$gene}{variants}{$varName}{position} = $chrPos;
						$vars->{$gene}{variants}{$varName}{base_change} = $baseChange;

						$vars->{$gene}{variants}{$varName}{tumour_freq} = $tFreq;
						unless ($tFreq eq "NA")
						{
							$vars->{$gene}{variants}{$varName}{tumour_freq} = sprintf("%.3f",$tFreq);
						}

						$vars->{$gene}{variants}{$varName}{tumour_depth} = $tDepth;
						$vars->{$gene}{variants}{$varName}{normal_freq} = sprintf("%.3f",$nFreq);
						$vars->{$gene}{variants}{$varName}{normal_depth} = $nDepth;
						$vars->{$gene}{variants}{$varName}{nuc_context} = $nucContext;
						$vars->{$gene}{variants}{$varName}{aa_context} = $aaContext;
						$vars->{$gene}{variants}{$varName}{dbsnp} = $id;
						$vars->{$gene}{variants}{$varName}{cosmic} = $COSMIC;
						$vars->{$gene}{variants}{$varName}{rarity} = $rarity;
						$vars->{$gene}{variants}{$varName}{"1000G_all"} = $anno{'1000g2015aug_all'};
						$vars->{$gene}{variants}{$varName}{ExAC_all} = $anno{ExAC_ALL};
#						$vars->{$gene}{variants}{$varName}{ESP6500siv2_all} = $anno{ESP6500siv2_ALL};
						$vars->{$gene}{variants}{$varName}{"gnomad_all"} = $anno{AF};
#						$vars->{$gene}{variants}{$varName}{cadd_phred} = $anno{CADD_phred};
						$vars->{$gene}{variants}{$varName}{clinvar} = $anno{CLINSIG};
						
						
					}


				}
			}
		}

		close FILE;
		#germline_snv_count germline_indel_count germline_common_count germline_rare_count germline_novel_count germline_titv_ratio germline_missense_count germline_nonsense_count
		$data->{germline_snv_count} = $snvCount;
		$data->{germline_indel_count} = $delCount + $insCount;
		$data->{germline_common_count} = $commonCount;
		$data->{germline_rare_count} = $rareCount;
		$data->{germline_novel_count} = $novelCount;
		if (($snvCount - $tiCount) > 0)
		{
			$data->{germline_titv_ratio} = $tiCount / ($snvCount - $tiCount);
		}
		$data->{germline_missense_count} = $missenseCount;
		$data->{germline_nonsense_count} = $nonsenseCount;
		$data->{germline_noncoding_count} = $noncodingCount;

	}
	else
	{
		warn "$file doesn't exist\n";
	}
	#print Dumper $data;
}



sub parseSV
{
	my $file = shift;
	my $data = shift;
	my $vars = shift;
	my $genes = shift;
	my $synonym = shift;

	my ($l,@f,$gene,$CT,$ALT,$chr1,$pos1,$chr2,$pos2,$type,$split1,$near1,$split2,$near2,$varName,$uniq1,$uniq2);
	my (%reported,%hash);

	my $svCount = 0;

	my $delCount = 0;
	my $dupCount = 0;
	my $invCount = 0;
	my $traCount = 0;

	my $delGeneBPCount = 0;
	my $dupGeneBPCount = 0;
	my $invGeneBPCount = 0;
	my $traGeneBPCount = 0;

	my @chroms = qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY/;		# no M

	my %prettyType = (
		"DEL" => "deletion",
		"DUP" => "duplication",
		"INV" => "inversion",
		"TRA" => "translocation"
	);

	my %svPerChr = (
		"chr1" => 0,
		"chr2" => 0,
		"chr3" => 0,
		"chr4" => 0,
		"chr5" => 0,
		"chr6" => 0,
		"chr7" => 0,
		"chr8" => 0,
		"chr9" => 0,
		"chr10" => 0,
		"chr11" => 0,
		"chr12" => 0,
		"chr13" => 0,
		"chr14" => 0,
		"chr15" => 0,
		"chr16" => 0,
		"chr17" => 0,
		"chr18" => 0,
		"chr19" => 0,
		"chr20" => 0,
		"chr21" => 0,
		"chr22" => 0,
		"chrX" => 0,
		"chrY" => 0,
		"chrM" => 0
	);

	my %chrLength = ( ##AMY updated to GRCh38 https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38, also saved in /.mounts/labs/PCSI/references/GRCh38/dna/size/GRCh38_size.bed
		"chr1" => 248956422, #249250621,
		"chr2" => 242193529, #243199373,
		"chr3" => 198295559, #198022430,
		"chr4" => 190214555, #191154276,
		"chr5" => 181538259, #180915260,
		"chr6" => 170805979, #171115067,
		"chr7" => 159345973, #159138663,
		"chrX" => 156040895, #155270560,
		"chr8" => 145138636, #146364022,
		"chr9" => 138394717, #141213431,
		"chr10" => 133797422, #135534747,
		"chr11" => 135086622, #135006516,
		"chr12" => 133275309, #133851895,
		"chr13" => 114364328, #115169878,
		"chr14" => 107043718, #107349540,
		"chr15" => 101991189, #102531392,
		"chr16" => 90338345, #90354753,
		"chr17" => 83257441, #81195210,
		"chr18" => 80373285, #78077248,
		"chr20" => 64444167, #63025520,
		"chrY" => 57227415, #59373566,
		"chr19" => 58617616, #59128983,
		"chr22" => 50818468, #51304566,
		"chr21" => 46709983, #48129895,
		"chrM" => 16571
	);

	my $size;
	my %svSizeBins;
	for my $type (qw/DEL DUP INV/)
	{
		for my $bin (qw/100 1k 10k 100k 1m 10m >10m/)
		{
			$svSizeBins{$type}{$bin} = 0;
		}
	}

	if (-e $file)
	{
		open (FILE, $file) or die "Couldn't open $file\n";
	
		while ($l = <FILE>)
		{
			chomp $l;
			if ($l =~ /^chrom1/)
			{
			}
			else
			{
				@f = split(/\t/, $l);
				
				$svCount++;

				$chr1 = $f[0];
				$pos1 = $f[1];
				$chr2 = $f[2];
				$pos2 = $f[3];
				$type = $f[4];
				$split1 = $f[7];
				$split2 = $f[8];
				$near1 = $f[9];
				$near2 = $f[10];
	
				$svPerChr{$chr1}++;
				$svPerChr{$chr2}++;
	
				if ($type eq 'BND' and $chr1 ne $chr2) {
					$type = "TRA";
				}	
				
				if ($type eq 'INS') {
					next;
				}

				if ($type eq "DEL")
				{
					$delCount++;
					if (($split1 ne ".") or ($split2 ne "."))
					{
						$delGeneBPCount++;
					}
				}
				elsif ($type eq "DUP")
				{
					$dupCount++;
					if (($split1 ne ".") or ($split2 ne "."))
					{
						$dupGeneBPCount++;
					}
				}
				elsif ($type eq "INV")
				{
					$invCount++;
					if (($split1 ne ".") or ($split2 ne "."))
					{
						$invGeneBPCount++;
					}
				}
				elsif ($type eq "TRA")
				{
					$traCount++;
					if (($split1 ne ".") or ($split2 ne "."))
					{
						$traGeneBPCount++;
					}
				}
	
				unless ($type eq "TRA")
				{
					$size = abs($pos2 - $pos1);
					if ($size <= 100)
					{
						$svSizeBins{$type}{100}++;
					}
					elsif ($size <= 1000)
					{
						$svSizeBins{$type}{'1k'}++;
					}
					elsif ($size <= 10000)
					{
						$svSizeBins{$type}{'10k'}++;
					}
					elsif ($size <= 100000)
					{
						$svSizeBins{$type}{'100k'}++;
					}
					elsif ($size <= 1000000)
					{
						$svSizeBins{$type}{'1m'}++;
					}
					elsif ($size <= 10000000)
					{
						$svSizeBins{$type}{'10m'}++;
					}
					else
					{
						$svSizeBins{$type}{'>10m'}++;
					}
				}


				$varName = "SV $type $chr1:$pos1-$chr2:$pos2";

				$split1 =~ s/,/|/g;
				$split2 =~ s/,/|/g;
				$near1 =~ s/,/|/g;
				$near2 =~ s/,/|/g;

	
				for $gene (split(/\|/, $split1))
				{
					unless ($gene eq ".")
					{
						$gene = findBestGeneSymbol($gene, $vars, $synonym, $chr1, "RefSeq (SV)", "warn");

						$vars->{$gene}{variants}{$varName}{mutation_class} = "somatic sv";
						$vars->{$gene}{variants}{$varName}{mutation_type} = "$prettyType{$type} breakpoint";
						$vars->{$gene}{variants}{$varName}{position} = "$chr1:$pos1-$chr2:$pos2";
						
						$uniq2 = "";
						%hash = ();
						for my $fusion (split(/\|/, "$split2|$near2"))
						{
							unless (($gene eq $fusion) or ($fusion eq "."))
							{
								$hash{$fusion}++;
							}
						}
						for my $fusion (sort keys %hash)
						{
							$uniq2 .= "|$fusion";
							if (exists $vars->{$gene}{cosmic_census_types}{translocation})
							{
								if (exists $vars->{$gene}{cosmic_census_fusions}{fusions}{$fusion})
								{
									$vars->{$gene}{variants}{$varName}{cosmic_census_flag} = "cosmic_mutation";
								}
							}
						}
						$uniq2 =~ s/^\|//;
						if ($uniq2 eq "")
						{
							$uniq2 = ".";
						}
						$vars->{$gene}{variants}{$varName}{fusion_genes} = $uniq2;

						$reported{$varName}++;
					}
				}
	
				for $gene (split(/\|/, $split2))
				{
					unless ($gene eq ".")
					{
						$gene = findBestGeneSymbol($gene, $vars, $synonym, $chr2, "RefSeq (SV)", "warn");

						$vars->{$gene}{variants}{$varName}{mutation_class} = "somatic sv";
						$vars->{$gene}{variants}{$varName}{mutation_type} = "$prettyType{$type} breakpoint";
						$vars->{$gene}{variants}{$varName}{position} = "$chr1:$pos1-$chr2:$pos2";

						$uniq1 = "";
						%hash = ();
						for my $fusion (split(/\|/, "$split1|$near1"))
						{
							unless (($gene eq $fusion) or ($fusion eq "."))
							{
								$hash{$fusion}++;
							}
						}
						for my $fusion (sort keys %hash)
						{
							$uniq1 .= "|$fusion";
							if (exists $vars->{$gene}{cosmic_census_types}{translocation})
							{
								if (exists $vars->{$gene}{cosmic_census_fusions}{fusions}{$fusion})
								{
									$vars->{$gene}{variants}{$varName}{cosmic_census_flag} = "cosmic_mutation";
								}
							}
						}
						$uniq1 =~ s/^\|//;
						if ($uniq1 eq "")
						{
							$uniq1 = ".";
						}
						$vars->{$gene}{variants}{$varName}{fusion_genes} = $uniq1;



						$reported{$varName}++;
					}
				}

				if (($near1 ne ".") and ($near2 ne ".") and (!exists $reported{$varName}) and ($near1 ne $near2))
				{
					for $gene (split(/\|/, $near1))
					{
						$gene = findBestGeneSymbol($gene, $vars, $synonym, $chr1, "RefSeq (SV)", "warn");

						$vars->{$gene}{variants}{$varName}{mutation_class} = "somatic sv";
						$vars->{$gene}{variants}{$varName}{mutation_type} = "$prettyType{$type} potential fusion";
						$vars->{$gene}{variants}{$varName}{position} = "$chr1:$pos1-$chr2:$pos2";
						$vars->{$gene}{variants}{$varName}{fusion_genes} = $near2;
						if (exists $vars->{$gene}{cosmic_census_types}{translocation})
						{
							for my $fusion (split(/\|/, $near2))
							{
								if (exists $vars->{$gene}{cosmic_census_fusions}{fusions}{$fusion})
								{
									$vars->{$gene}{variants}{$varName}{cosmic_census_flag} = "cosmic_mutation";
								}
							}
						}
					}
					for $gene (split(/\|/, $near2))
					{
						$gene = findBestGeneSymbol($gene, $vars, $synonym, $chr2, "RefSeq (SV)", "warn");
						
						$vars->{$gene}{variants}{$varName}{mutation_class} = "somatic sv";
						$vars->{$gene}{variants}{$varName}{mutation_type} = "$prettyType{$type} potential fusion";
						$vars->{$gene}{variants}{$varName}{position} = "$chr1:$pos1-$chr2:$pos2";
						$vars->{$gene}{variants}{$varName}{fusion_genes} = $near1;
						if (exists $vars->{$gene}{cosmic_census_types}{translocation})
						{
							for my $fusion (split(/\|/, $near1))
							{
								if (exists $vars->{$gene}{cosmic_census_fusions}{fusions}{$fusion})
								{
									$vars->{$gene}{variants}{$varName}{cosmic_census_flag} = "cosmic_mutation";
								}
							}
						}
					}
				}
	
			}
		}
	
		my $waddellClass;
		my @rankedSvPerChr;
		my $svPerMinimum;
	
		if ($svCount < 50)
		{
			$waddellClass = "stable";
		}
		elsif ($svCount > 200)
		{
			$waddellClass = "unstable";
		}
		else	# scattered or locally rearranged
		{
			for my $chr (@chroms)
			{
				$svPerChr{$chr} = $svPerChr{$chr} / $chrLength{$chr};
				push(@rankedSvPerChr, $svPerChr{$chr});
			}
	
			@rankedSvPerChr = sort { $a <=> $b } @rankedSvPerChr;
	
	
			my $svPerMinimum = 5 * (($rankedSvPerChr[17] + $rankedSvPerChr[18]) / 2);
	
			$waddellClass = "scattered";
			for my $chr (@chroms)
			{
				if ($svPerChr{$chr} > $svPerMinimum)
				{
					$waddellClass = "locally rearranged";
				}
			}
		}
	
	
		$data->{sv_file} = $file;
		$data->{sv_count} = $svCount;
		$data->{waddell} = $waddellClass;
	
		for my $type (qw/DEL DUP INV/)
		{
			for my $bin (qw/100 1k 10k 100k 1m 10m >10m/)
			{
				$data->{"$type:$bin"} = $svSizeBins{$type}{$bin};
			}
		}
	
		$data->{sv_del_count} = $delCount;
		$data->{sv_dup_count} = $dupCount;
		$data->{sv_inv_count} = $invCount;
		$data->{sv_tra_count} = $traCount;

		$data->{sv_del_bp_gene_count} = $delGeneBPCount;
		$data->{sv_dup_bp_gene_count} = $dupGeneBPCount;
		$data->{sv_inv_bp_gene_count} = $invGeneBPCount;
		$data->{sv_tra_bp_gene_count} = $traGeneBPCount;
	}
	else
	{
		warn "$file doesn't exist\n";
	}
	#print Dumper $data;
}


sub parseSegs
{
	my $file = shift;
	my $data = shift;
	my $vars = shift;
	my $genes = shift;
	my $synonym = shift;

#	my $ploidyFoldAmp = ($data->{ploidy} * 4) - 0.5;
	my $ploidyFoldAmp = ($data->{ploidy} * 3);
	my $multiSeg;
	my $homdel;

	my $position;

	if (-e $file)
	{
		my $iter;
		my $tabix = Tabix->new(-data => $file, -index => "$file.tbi");
	
		my ($l,@f);
		my ($chr,$start,$end);
	
	
		for my $gene (keys %{ $genes })
		{
			$genes->{$gene} =~ /(.*):(.*)-(.*)/;
			$chr = $1;
			$start = $2;
			$end = $3;

			$multiSeg = 0;
			$homdel = 0;
			$position = "";
	
			$gene = findBestGeneSymbol($gene, $vars, $synonym, $chr, "RefSeq (SEG)", "no warn");
	
			$iter = $tabix->query($chr,$start,$end);
			if (defined $iter->{"_"})
			{
				while ($l = $tabix->read($iter))
				{
						# chr1    23918001        33076000        9158000 p       PCSI_0652_Pa_P_526      0.875   0.981165109381573       FALSE   0.4221  0.37697585636584        1.61472711055046        1.1

					@f = split(/\t/, $l);
					unless ($f[11] eq "NA")
					{
						if ($f[11] <= 0)
						{
							$homdel = 1;
						}
					}
					if (exists $vars->{$gene}{copy_number})
					{

						$vars->{$gene}{copy_number} .= "|" . sprintf("%.3f",$f[11]);
						$vars->{$gene}{maf_mean} .= "|" . sprintf("%.3f",$f[9]);
						unless ($f[10] eq "NA")
						{
							$vars->{$gene}{maf_p} .= "|" . sprintf("%.3f",$f[10]);
						}
						else
						{
							$vars->{$gene}{maf_p} .= "|NA";
						}
						$vars->{$gene}{ab_counts} .= "|" . $f[12];

						$position .= "|$f[0]:$f[1]-$f[2]";
						$multiSeg = 1;
					}
					else
					{
						$vars->{$gene}{copy_number} = sprintf("%.3f",$f[11]);

						unless ($f[9] eq "NA")
						{
							$vars->{$gene}{maf_mean} = sprintf("%.3f",$f[9]);
						}
						else
						{
							$vars->{$gene}{maf_mean} = $f[9];
						}
						unless ($f[10] eq "NA")
						{
							$vars->{$gene}{maf_p} = sprintf("%.3f",$f[10]);
						}
						else
						{
							$vars->{$gene}{maf_p} = $f[10];
						}
						$vars->{$gene}{ab_counts} = $f[12];
						$position = "$f[0]:$f[1]-$f[2]";
					}
				}

				if ((exists $vars->{$gene}{ab_counts}) and (exists $vars->{$gene}{copy_number}))
				{
					if ($multiSeg == 0)
					{
						if ($vars -> {$gene}{copy_number} > $ploidyFoldAmp and $vars -> {$gene}{copy_number} > 5.5)
						{
							$vars->{$gene}{variants}{"CNV strong amplification"}{mutation_class} = "somatic cnv";
							$vars->{$gene}{variants}{"CNV strong amplification"}{mutation_type} = "strong amplification";
							$vars->{$gene}{variants}{"CNV strong amplification"}{position} = $position;
							if (exists $vars->{$gene}{cosmic_census_types}{"strong amplification"})
							{
								$vars->{$gene}{variants}{"CNV strong amplification"}{cosmic_census_flag} = "cosmic_mutation";
							}
						}
						if (($vars->{$gene}{ab_counts} =~ /0\.0/) or ($homdel == 1) or ($vars->{$gene}{copy_number} < 0.2)) ##added so copy numbers < 0.2 will count as homozygous deletion as well 
						{
							$vars->{$gene}{variants}{"CNV homozygous deletion"}{mutation_class} = "somatic cnv";
							$vars->{$gene}{variants}{"CNV homozygous deletion"}{mutation_type} = "homozygous deletion";
							$vars->{$gene}{variants}{"CNV homozygous deletion"}{position} = $position;
							if (exists $vars->{$gene}{cosmic_census_types}{"homozygous deletion"})
							{
								$vars->{$gene}{variants}{"CNV homozygous deletion"}{cosmic_census_flag} = "cosmic_mutation";
							}
						}
					}
					###COPY NUMBER HER2
					elsif ($multiSeg == 1) {
						my @copy_numbers = split /\|/, $vars->{$gene}{copy_number};
						my @ab_counts = split /\|/, $vars->{$gene}{ab_counts};
						if ( all { $_ > $ploidyFoldAmp } @copy_numbers and all { $_ > 5.5 } @copy_numbers ){
							$vars->{$gene}{variants}{"CNV strong amplification"}{mutation_class} = "somatic cnv";
							$vars->{$gene}{variants}{"CNV strong amplification"}{mutation_type} = "strong amplification";
							$vars->{$gene}{variants}{"CNV strong amplification"}{position} = $position;
							if (exists $vars->{$gene}{cosmic_census_types}{"strong amplification"})
							{
								$vars->{$gene}{variants}{"CNV strong amplification"}{cosmic_census_flag} = "cosmic_mutation";
							}
							
						}
						if ( all { $_ =~ /0\.0/ } @ab_counts or ($homdel == 1) or all { $_ < 0.2 } @copy_numbers) ##added so copy numbers < 0.2 will count as homozygous deletion as well 
						{
							$vars->{$gene}{variants}{"CNV homozygous deletion"}{mutation_class} = "somatic cnv";
							$vars->{$gene}{variants}{"CNV homozygous deletion"}{mutation_type} = "homozygous deletion";
							$vars->{$gene}{variants}{"CNV homozygous deletion"}{position} = $position;
							if (exists $vars->{$gene}{cosmic_census_types}{"homozygous deletion"})
							{
								$vars->{$gene}{variants}{"CNV homozygous deletion"}{cosmic_census_flag} = "cosmic_mutation";
							}
						}

					}
				}
			}
			else
			{
			}
		}
		$data->{seg_file} = $file;
	}
	else
	{
		warn "$file doesn't exist\n";
	}

}

sub parseParams
{
	my $file = shift;
	my $data = shift;

	my $l;

	$data->{param_file} = $file;

	if (-e $file)
	{
		open (FILE, $file) or die "Couldn't open $file\n";
	
		$l = <FILE>;	# header
		$l = <FILE>;
		chomp $l;
		my @f = split(/ /, $l);
		$data->{cellularity} = sprintf("%.3f",$f[3]);
		$data->{ploidy} = sprintf("%.3f",$f[4]);
		close FILE;
	}
	else
	{
		warn "$file doesn't exist\n";
	}
}

sub parseNeo
{
	my $path = shift;
	my $data = shift;
	my $vars = shift;

	$data->{neo_path} = $path;

	my $l;

	my %neoList;
	my %varList;

	my $ls = `ls $path/$sample.*.txt`;
	chomp $ls;

	for my $file (split(/\n/, $ls))
	{
		open (FILE, $file) or die "Couldn't open $file\n";
	
		while ($l = <FILE>)
		{
			#     0  HLA-A*02:01    AAERQELGG         PEPLIST         0.008     45729.24    50.00
			if ($l =~ /.*?(HLA-.*?)  *(.*?)  *?PEPLIST.* <=/)
			{
				$neoList{$2}{$1}++;
			}
			elsif ($l =~ /^Protein PEPLIST\. Allele (.*?)\. Number of high binders (.*?)\. Number of weak binders (.*?)\. Number of peptides .*?$/)
			{
				$data->{neo_antigens} += $2 + $3;
				$data->{neo_antigens_weak} += $3;
				$data->{neo_antigens_strong} += $2;
				$data->{hla_types} .= "$1|";
			}
		}

		close FILE;
	}

	# read peptide map to collect variants for neo-antigens
	open (FILE, "$path/$sample.peptideMap") or die "Couldn't open $path/$sample.peptideMap\n";

	my ($chr,$pos,$ref,$alt,$pep);
	while ($l = <FILE>)
	{
		chomp $l;
		($chr,$pos,$ref,$alt,$pep) = split(/\t/, $l);

		if (exists $neoList{$pep})
		{
			for my $hla (sort keys %{ $neoList{$pep} })
			{
				$varList{"$chr:$pos\t$ref>$alt"}{$pep}{$hla}++;
			}
		}
	}

	close FILE;
	
	
			if (exists $vars->{ARID1A}{variants})
			{
				warn "ARID1A variants exist before matching neo\n";
			}

	# go through variant hash (somatic noncoding, specifically) and add field for predicted neo-antigens to those that match %neoList and the peptide map
	
	for my $gene (keys %{ $vars })
	{
		if (exists $vars->{$gene}{variants})
		{
			for my $varName (keys %{ $vars->{$gene}{variants} })
			{
				if (($vars->{$gene}{variants}{$varName}{mutation_class} eq "somatic snv") or ($vars->{$gene}{variants}{$varName}{mutation_class} eq "somatic indel"))
				{
					$pos = $vars->{$gene}{variants}{$varName}{position};
					$alt = $vars->{$gene}{variants}{$varName}{base_change};
	
					if (exists $varList{"$pos\t$alt"})
					{
						for $pep (sort keys %{ $varList{"$pos\t$alt"} })
						{
							for my $hla (sort keys %{ $varList{"$pos\t$alt"}{$pep} })
							{
								$vars->{$gene}{variants}{$varName}{neoantigen} .= "$hla>$pep;";
							}
						}
						$vars->{$gene}{variants}{$varName}{neoantigen} =~ s/;$//;
					}
				}
			}
		}
	}




	$data->{hla_types} =~ s/\|$//;

}

sub parseXenome
{
	my $files = shift;
	my $data = shift;

	$data->{xenome_logs} = $files;

	my $l;

	my $ls = `ls $files`;
	chomp $ls;

	my ($sum,$count);

	for my $file (split(/\n/, $ls))
	{
		open (FILE, $file) or die "Couldn't open $file\n";

		while ($l = <FILE>)
		{
			if ($l =~ /^.*\t(.*?)\tmouse$/)
			{
				$sum += $1;
				$count++;
			}
		}
	}

	if ($count > 0)
	{
		$data->{mouse_content} = $sum / $count;
	}
	else
	{
		$data->{mouse_content} = "NA";
	}
	return;
}

# my %genePos = lookUpGenePositions($refSeqFile, \@genesOfInterest);
sub lookUpGenePositions
{
	my $file = shift;
	my $vars = shift;
	my $synonym = shift;
	my %genePos;

	open (FILE, $file) or die "Couldn't open $file\n";

	my ($l, $gene, $chr, $start, $end, $info);

	while ($l = <FILE>)
	{
		chomp $l;
		($chr, $start, $end, $gene) = split(/\t/, $l);

		$gene = findBestGeneSymbol($gene, $vars, $synonym, $chr, "RefSeq", "no warn");

		$genePos{$gene} = "$chr:$start-$end";
		$vars->{$gene}{gene} = $gene;
		$vars->{$gene}{gene_position} = "$chr:$start-$end";
	}

	return %genePos;
}


sub findBestGeneSymbol
{
	my $g = shift;
	my $vars = shift;
	my $synonym = shift;
	my $chr = shift;
	my $origin = shift;
	my $doWarn = shift;
	my $synChr;

	my $lastSynChr;

	if (exists $vars->{$g})
	{
		return $g;
	}
	elsif (exists $synonym->{$g})
	{
		for $synChr (keys %{ $synonym->{$g} })
		{
			$lastSynChr = $synChr;
			if ($synChr eq $chr)
			{
				unless ($synonym->{$g}{$synChr} =~ /;/)
				{
					return $synonym->{$g}{$synChr};
				}
				else
				{
					if ($doWarn eq "warn")
					{
						warn " $origin: multiple same chr synonyms for $g on $chr: $synonym->{$g}{$synChr}\n";
					}
					return $synonym->{$g}{$synChr};
				}
			}
		}
		if ($doWarn eq "warn")
		{
			warn " $origin: no synonyms for $g on $chr, using original name\n";
		}
		return $g;
	}
	else
	{
		if ($doWarn eq "warn")
		{
			warn " $origin: no match or synonym entry for $g, using original name\n";
		}
		return $g;
	}
}


sub parseCosmicSigNNLS
{
	my $file = shift;
	my $data = shift;

	$data->{cosmicSigNNLS_file} = $file;

	my $l;

	my @header;
	my %row;

	my %keep;

	if (-e $file)
	{
		open (FILE, $file) or die "Couldn't open $file\n";
	
		$l = <FILE>;
		@header = split(/ /, $l);
	
		my $i = 0;
		$l = <FILE>;
		for my $f (split(/ /, $l))
		{
			$row{$header[$i]} = $f;
			$i++;
		}
		close FILE;
	
	
		my $prefix = "csnnls";
		%keep = (
			"Signature.1" => "${prefix}_sig1",
			"Signature.2" => "${prefix}_sig2",
			"Signature.3" => "${prefix}_sig3",
			"Signature.5" => "${prefix}_sig5",
			"Signature.6" => "${prefix}_sig6",
			"Signature.8" => "${prefix}_sig8",
			"Signature.13" => "${prefix}_sig13",
			"Signature.17" => "${prefix}_sig17",
			"Signature.18" => "${prefix}_sig18",
			"Signature.20" => "${prefix}_sig20",
			"Signature.26" => "${prefix}_sig26",
			"Residuals" => "${prefix}_residuals",
			"N.Mutations" => "${prefix}_n.mutations",
		);
	}

	for my $f (@header)
	{
		if (exists $row{$f})
		{
			if (exists $keep{$f})
			{
				$data->{$keep{$f}} = $row{$f};
			}
		}
	}


}




sub doANNOVARlookup
{
	my $chr = shift;
	my $pos = shift;
	my $ref = shift;
	my $alt = shift;
	my $tabix = shift;
	my $dataRef = shift;
	my $info = shift;

	my @header = qw/Chr Start End Ref Alt Func.refGene Gene.refGene GeneDetail.refGene ExonicFunc.refGene AAChange.refGene avsnp150 1000g2015aug_all 1000g2015aug_afr 1000g2015aug_eas 1000g2015aug_eur 1000g2015aug_sas ExAC_ALL ExAC_AFR ExAC_AMR ExAC_EAS ExAC_FIN ExAC_NFE ExAC_OTH ExAC_SAS AF AF_popmax AF_male AF_female AF_raw AF_afr AF_sas AF_amr AF_eas AF_nfe AF_fin AF_asj AF_oth non_topmed_AF_popmax non_neuro_AF_popmax non_cancer_AF_popmax controls_AF_popmax CLNALLELEID CLNDN CLNDISDB CLNREVSTAT CLNSIG Otherinfo/;

	unless ($tabix eq "null")
	{
		my $iter = $tabix->query($chr,$pos-5,$pos+10);
		my $l;
		my @f;
	
		if (defined $iter->{"_"}) {
		   while ($l = $tabix->read($iter)) {
			@f = split(/\t/, $l);
			if (length($alt) > length($ref) and $ref ne '-') {
			   $ref = "-";
			   $alt = substr($alt, 1);
			} 
			if (length($ref) > length($alt) and $alt ne '-') {
			   $alt = "-";
			   $ref = substr($ref, 1);
			   $pos = $pos + 1;
			}
			if (length($ref) > length($alt) and substr($ref, 0, 1) eq substr($alt, 0, 1)) {
			   $alt = "-";
			   $ref = substr($ref, 1);
			   $pos = $pos + 1;
			}
			if (($f[0] eq $chr) and ($f[1] eq $pos) and ($f[3] eq $ref or $f[52] eq $alt) and ($f[4] eq $alt or $f[53] eq $alt)) {
			   for (my $i = 5; $i < scalar(@header); $i++){
				if ($f[$i] eq ".") {
				   $f[$i] = "NA";
				}
				$f[$i] =~ s/,/./g;
				$dataRef->{$header[$i]} = $f[$i];
			   }
			} else {
#				print "found $f[0], $f[1], $f[52], $f[53] does not match search $chr,$pos,$ref,$alt\n";
			}
		  }		
		} else {
			my $search = $pos - 5;
			my $search2 = $pos + 10;
			print "did not find variant $chr, $pos, $ref, $alt when searching from $search to $search2\traw:$info\n";
			die;
		}
	   
	}

	return;
}



sub doCADDlookup
{

	my $chr = shift;
	my $pos = shift;
	my $ref = shift;
	my $alt = shift;
	my $tabix = shift;
	my $dataRef = shift;

	unless ($tabix eq "null")
	{

		$chr =~ s/chr//;
	
		my @header = qw/Chrom  Pos     Ref     Anc     Alt     Type    Length  isTv    isDerived       AnnoType        Consequence     ConsScore       ConsDetail      GC      CpG     mapAbility20bp  mapAbility35bp  scoreSegDup     priPhCons       mamPhCons       verPhConspriPhyloP        mamPhyloP       verPhyloP       GerpN   GerpS   GerpRS  GerpRSpval      bStatistic      mutIndex        dnaHelT dnaMGW  dnaProT dnaRoll mirSVR-Score    mirSVR-E        mirSVR-Aln      targetScan      fitCons cHmmTssA        cHmmTssAFlnk      cHmmTxFlnk      cHmmTx  cHmmTxWk        cHmmEnhG        cHmmEnh cHmmZnfRpts     cHmmHet cHmmTssBiv      cHmmBivFlnk     cHmmEnhBiv      cHmmReprPC      cHmmReprPCWk    cHmmQuies       EncExp  EncH3K27Ac      EncH3K4Me1      EncH3K4Me3      EncNucleo EncOCC  EncOCCombPVal   EncOCDNasePVal  EncOCFairePVal  EncOCpolIIPVal  EncOCctcfPVal   EncOCmycPVal    EncOCDNaseSig   EncOCFaireSig   EncOCpolIISig   EncOCctcfSig    EncOCmycSig     Segway  tOverlapMotifs  motifDist       motifECount     motifEName        motifEHIPos     motifEScoreChng TFBS    TFBSPeaks       TFBSPeaksMax    isKnownVariant  ESP_AF  ESP_AFR ESP_EUR TG_AF   TG_ASN  TG_AMR  TG_AFR  TG_EUR  minDistTSS      minDistTSE      GeneID  FeatureID       CCDS    GeneName        cDNApos   relcDNApos      CDSpos  relCDSpos       protPos relProtPos      Domain  Dst2Splice      Dst2SplType     Exon    Intron  oAA     nAA     Grantham        PolyPhenCat     PolyPhenVal     SIFTcat SIFTval RawScore        PHRED/;
	
		my $iter = $tabix->query($chr,$pos,$pos+1);
	
		my $l;
		my @f;
	
		$dataRef = {};
	
		if (defined $iter->{"_"})
		{
			while ($l = $tabix->read($iter))
			{
				@f = split(/\t/, $l);
				if (($f[0] eq $chr) and ($f[1] eq $pos) and ($f[2] eq $ref) and ($f[4] eq $alt))
				{
					for (my $i = 5; $i < scalar(@header); $i++)
					{
						$dataRef->{$header[$i]} = $f[$i];
					}
				}
			}
		}
	}

	return;
}



sub comparePos
{
	my ($chr1, $chr2, $start1, $start2);

	if ((exists $vars{$a}{gene_position}) and (exists $vars{$b}{gene_position}))
	{
		if ($vars{$a}{gene_position} =~ /^chr(.*?):(.*?)-.*$/)
		{
			$chr1 = $1;
			$start1 = $2;
		}
		if ($vars{$b}{gene_position} =~ /^chr(.*?):(.*?)-.*$/)
		{
			$chr2 = $1;
			$start2 = $2;
		}
	
		if ($chr1 eq "X")
		{
			$chr1 = 23;
		}
		if ($chr2 eq "X")
		{
			$chr2 = 23;
		}
	
		if ($chr1 eq "Y")
		{
			$chr1 = 24;
		}
		if ($chr2 eq "Y")
		{
			$chr2 = 24;
		}
	

		if ($chr1 eq 'M' or $chr2 eq 'M') {
			return 0;
		}
		if ($chr1 > $chr2)
		{
			return 1;
		}
		elsif ($chr1 < $chr2)
		{
			return -1;
		}
		elsif ($start1 > $start2)
		{
			return 1;
		}
		elsif ($start1 < $start2)
		{
			return -1;
		}
		else
		{
			return 0;
		}
	}
	else
	{
		return 0;
	}
}


sub parseHRDetect
{
	my $file = shift;
	my $data = shift;

	$data->{hrdetect_file} = $file;

	my ($l, $i);
	my @header = qw/patient sample deletion_microhomology_proportion snv_signature_3 snv_signature_8 rearrangement_signature_3 rearrangement_signature_5 hrd_score loh tai lst HRDetect/;
	my %row;

	if (-e $file)
	{
		open (FILE, $file) or die "Couldn't open $file\n";

		while ($l = <FILE>)
		{
			chomp $l;
			unless ($l =~ /^sample_name/)
			{
				%row = ();
				$i = 0;
				for my $var (split(/,/, $l))
				{
					$row{$header[$i]} = $var;
					$i++;
				}

				$data->{"hrdetect_score"} = $row{HRDetect};
				$data->{"hrdetect_myriad_score"} = $row{hrd_score};

				for my $var (qw/deletion_microhomology_proportion snv_signature_3 snv_signature_8 rearrangement_signature_3 rearrangement_signature_5 loh tai lst/)
				{
					$data->{"hrdetect_$var"} = $row{$var};
				}
			}
		}
	}
	return;
}


sub parseTDP
{
	my $file = shift;
	my $data = shift;

	$data->{tdp_file} = $file;

	my ($l, $i);
	my @header = qw/sample_name chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 total score/;
	my %row;

	if (-e $file)
	{
		open (FILE, $file) or die "Couldn't open $file\n";

		while ($l = <FILE>)
		{
			chomp $l;
			unless ($l =~ /^sample_name/)
			{
				%row = ();
				$i = 0;
				for my $var (split(/,/, $l))
				{
					$row{$header[$i]} = $var;
					$i++;
				}

				$data->{"tdp_score"} = $row{score};
			}
		}
	}
	return;
}



sub parseAlexandrov
{
	my $file = shift;
	my $data = shift;

	$data->{alexandrov_file} = $file;

	my ($l, $i);
	my @header = qw/sample sig1 sig2 sig3 sig5 sig6 sig8 sig13 sig17 sig18 sig20 sig26 accuracy class etiology/;
	my %row;

	for my $var (qw/sig1 sig2 sig3 sig5 sig6 sig8 sig13 sig17 sig18 sig20 sig26 accuracy class etiology/)
	{
		$data->{"alexandrov_$var"} = "NA";
	}

	if (-e $file)
	{
		open (FILE, $file) or die "Couldn't open $file\n";

		while ($l = <FILE>)
		{
			chomp $l;
			unless ($l =~ /^Sample/)
			{
				%row = ();
				$i = 0;
				for my $var (split(/\t/, $l))
				{
					$row{$header[$i]} = $var;
					$i++;
				}

				if ($row{sample} eq $sample)
				{
					for my $var (qw/sig1 sig2 sig3 sig5 sig6 sig8 sig13 sig17 sig18 sig20 sig26 accuracy class etiology/)
					{
						$data->{"alexandrov_$var"} = $row{$var};
					}
				}
			}
		}
	}
	return;
}


sub parseRNAsigs
{
	my $file = shift;
	my $data = shift;

	$data->{rna_sig_file} = $file;

	my ($l, $i);
	my @header = qw/library donor compass moffitt collisson bailey pre-pass study.id/;
	my %row;

	my $sample = $sample;
	$sample =~ s/_526//;

	my $rnaSample;

	if (-e $file)
	{
		open (FILE, $file) or die "Couldn't open $file\n";

		while ($l = <FILE>)
		{
			chomp $l;
			unless ($l =~ /^Name/ or $l !~ /[a-zA-Z0-9]/)
			{
				%row = ();
				$i = 0;
				for my $var (split(/\t/, $l))
				{
					$row{$header[$i]} = $var;
					$i++;
				}

				$rnaSample = $row{library};
				$rnaSample =~ s/_526//;
				
				if ($rnaSample eq $sample)
				{
					for my $var (qw/moffitt collisson bailey/) {
						if ($data->{tumour} =~ /526[82]/ or $l !~ /Pre-Pass/ or $var !~ /moffitt/) {			
								warn "can get $var for $sample\n";
							if (exists $data->{$var}) {
								print "Also found $data->{$var} for $sample ($row{library})\n";
	#							$data->{$var} .= "|$row{$var}";
							}
							else
							{
								$data->{$var} = $row{$var};
								print "Found $data->{$var} for $sample ($row{library})\n";
							}
						}
					}
				}
			}
		}
	}
	return;
}


sub parseCoverage
{
	my $path = shift;
	my $data = shift;

	my %coverageHash;
	my $l;

	for my $samp (keys %{ $path })
	{
		if (-e $path->{$samp})
		{
			open (FILE, $path->{$samp}) or die "Couldn't open $path->{$samp}\n";
			my $l = `cat $path->{$samp}`;
			chomp $l;
			my @lines = split /\n/, $l;
			for (my $i = 0; $i < scalar @lines; $i++) {
			   $l = $lines[$i];
			   if ($l =~ /GENOME_TERRITORY/) {			
				my ($genome, $mean_coverage, $sd, $median_coverage, @extra) = split /\t/, $lines[$i+1];
				print "$samp\t$median_coverage\n";
				$data{"$samp\_coverage"} = $median_coverage;
#			$data{"${samp}_error_rate"} = sprintf("%0.5f", $error);
#			$data{"${samp}_soft_clip"} = sprintf("%0.5f", $soft);
#			$data{"${samp}_reads_per_sp"} = sprintf("%0.5f", $rpsp);
			   }
			}

			close FILE;
		}
	}
	return;
}


sub parseSex
{
	my $file = shift;
	my $data = shift;

	my $inferredSex = "NA";

	open (FILE, $file) or die "Couldn't open $file\n";
	my $l = <FILE>;
	chomp $l;
	if ($l eq 'M') {
		$inferredSex = "Male";
	} elsif ($l eq 'F') {
		$inferredSex = "Female";
	}
	close FILE;
	
	$data->{inferred_sex} = $inferredSex;
	return;

}

sub parseGenotype
{
	my $paths = shift;
	my $data = shift;

	my $ls = `ls $paths`;
	chomp $ls;
	my ($l, $chr, $pos, $id, $call, $file1, $file2);
	my %calls;

	my @files = split(/\n/, $ls);

	for my $file (@files)
	{
		open (FILE, $file) or die "Couldn't open $file\n";
		while ($l = <FILE>)
		{
			chomp $l;
			($chr, $pos, $id, $call) = split(/\t/, $l);
			$calls{"$chr,$pos,$id"}{$file} = $call;
		}
		close FILE;
	}

	my $assayedCount;
	my $matchedCount;
	my $concordance;
	my $minConcordance = 999;

	for my $f1 (@files)
	{
		for my $f2 (@files)
		{
			$assayedCount = 0;
			$matchedCount = 0;
			for my $pos (keys %calls)
			{
				if ((exists $calls{$pos}{$f1}) and (exists $calls{$pos}{$f2}))
				{
					if (($calls{$pos}{$f1} ne "") and ($calls{$pos}{$f2} ne ""))
					{
						$assayedCount++;
						if ($calls{$pos}{$f1} eq $calls{$pos}{$f2})
						{
							$matchedCount++;
						}
					}
				}
			}

			if ($assayedCount > 0)
			{
				$concordance = $matchedCount / $assayedCount;
				if ($concordance < $minConcordance)
				{
					$minConcordance = $concordance;
				}
			}
		}
	}

	unless ($minConcordance == 999)
	{
		$data->{genotype_concordance} = $minConcordance;
	}
	return;

}


sub parseCosmic
{
	my $file = shift;
	my $vars = shift;
	my $synonym = shift;

	my ($l,$foundGene, $chr);
	my (@f, @header);
	my (%row, %synonyms);

	open (FILE, $file) or die "Couldn't open $file\n";

	my %mutTypes = (
		"A" => "strong amplification",
		"D" => "homozygous deletion",
		"F" => "frameshift",
		"Mis" => "nonsynonymous,nonframeshift",
		"N" => "stopgain",
		"S" => "splicing",
		"T" => "translocation",
	);

	while ($l = <FILE>)
	{
		chomp $l;

		# expecting: Gene Symbol	Name	Entrez GeneId	Genome Location	Chr Band	Somatic	Germline	Tumour Types(Somatic)	Tumour Types(Germline)	Cancer Syndrome	Tissue Type	Molecular Genetics	Role in Cancer	Mutation Types	Translocation Partner	Other Germline Mut	Other Syndrome	Synonyms
		#
		if ($l =~ /^Gene Symbol/)
		{
			@header = split(/\t/, $l);
		}
		else
		{
			$foundGene = 0;
			%row = ();
			@f = split(/\t/, $l);

			for (my $i = 0; $i < scalar(@f); $i++)
			{
				$row{$header[$i]} = $f[$i];
			}

			%synonyms = ();

			$synonyms{$row{'Gene Symbol'}}++;
			if (exists $row{Synonyms})
			{
				for my $name (split(/,/, $row{Synonyms}))
				{
					$synonyms{$name}++;
				}
			}

			$chr = $row{'Genome Location'};
			$chr =~ s/:.*//;
			$chr = "chr$chr";

			for my $name (keys %synonyms)
			{
				$name = findBestGeneSymbol($name, $vars, $synonym, $chr, "COSMIC Census", "no warn");
				
				if (exists $vars->{$name})
				{
					$foundGene = 1;

					# add cosmic data to gene
					$vars->{$name}{cosmic_census_data} = "";
					for my $h (@header)
					{
						if (exists $row{$h})
						{
							$vars->{$name}{cosmic_census_data} .= "$h:$row{$h}|";
						}
					}
					$vars->{$name}{cosmic_census_data} =~ s/\|$//;
					$vars->{$name}{cosmic_census_data} =~ s/,/;/g;

					for my $t (split(/, /, $row{'Mutation Types'}))
					{
						if (exists $mutTypes{$t})
						{
							for my $type (split(/,/, $mutTypes{$t}))
							{
								$vars->{$name}{cosmic_census_types}{$type}++;
							}
						}
					}

					for my $g (split(/, /, $row{'Translocation Partner'}))
					{
						$vars->{$name}{cosmic_census_fusions}{$g}++;
					}
				}
			}

			if ($foundGene == 0)
			{
				if (exists $row{Synonyms})
				{
					warn "Couldn't find HUGO gene for Cosmic census gene with the following names: $row{'Gene Symbol'},$row{Synonyms}\n";
				}
				else
				{
					warn "Couldn't find HUGO gene for Cosmic census gene with the following names: $row{'Gene Symbol'}\n";
				}
			}

		}
	}
}


sub parseLIMS
{
	my $url = shift;
	my $data = shift;

	my $donor = $data{donor};
	my $tumour = $data{tumour};
	unless ($donor =~ /EPPIC/) {
	$donor =~ s/^([A-Z0-9]+)(....)/$1_$2/;
	} else {
	$donor =~ s/EPPIC/EPPIC_/;
	}
	print "$donor\n";

	my $externalID = "NA";
	my $idRef;
	my $l;

	my $limsDump = "/.mounts/labs/PCSI/raw_data/oicr/lims/lims_identity_dump_20170927.json";
	my $prettyJson = "";

	open (FILE, $limsDump) or die "Couldn't open $limsDump\n";
	while ($l = <FILE>)
	{
		chomp $l;
		$prettyJson .= $l;
	}
	$idRef = decode_json($prettyJson);
	print "donor is $donor\n";

	for (my $i = 0; $i < scalar(@$idRef); $i++)
	{
		if ($idRef->[$i]{name} eq $donor)
		{
			for (my  $j = 0; $j < scalar(@{ $idRef->[$i]{attributes} }); $j++)
			{
				if ($idRef->[$i]{attributes}[$j]{name} eq "External Name")
				{
					$externalID = $idRef->[$i]{attributes}[$j]{value};
					$externalID =~ s/,/ /g;
				}
			}
		}
	}

	if ($donor eq "PCSI1146") {
		$externalID = "COMP-0268-M";
	}
	if (exists $external_ID{$tumour}) {
		$externalID = $external_ID{$tumour};
	}

	if ($externalID eq "NA")		# only check lims if the ID isn't in the dump
	{
		warn "Checking LIMS\n";
#		$prettyJson = `curl -X GET http://pinery-prod.hpc.oicr.on.ca:8080/pinery-ws-miso/samples?type=Identity`;
		if ($donor =~ /BTC/) {
			$prettyJson = `curl -X GET http://pinery.gsi.oicr.on.ca/samples?type=Identity&project=BTC`;
		} else {
			$prettyJson = `curl -X GET http://pinery.gsi.oicr.on.ca/samples?type=Identity`;
		}
		$prettyJson =~ s/\n//g;
		$idRef = decode_json($prettyJson);

		for (my $i = 0; $i < scalar(@$idRef); $i++)
		{
			if ($idRef->[$i]{name} eq $donor)
			{
				for (my  $j = 0; $j < scalar(@{ $idRef->[$i]{attributes} }); $j++)
			{
					if ($idRef->[$i]{attributes}[$j]{name} eq "External Name")
					{
						print "found external name\n";
						$externalID = $idRef->[$i]{attributes}[$j]{value};
						$externalID =~ s/,/ /g;
					}
				}
			}
		}
	}

	print "external $externalID\n";
	$data{external_id} = $externalID;
	#print Dumper $data;
}





# doHallmarkScoring(\%data, \%vars, \@dsbrGenes, \@mmrGenes);
sub doHallmarkScoring
{
	my $data = shift;
	my $vars = shift;
	my $dsbrGenes = shift;
	my $mmrGenes = shift;

	#     dsbr_snv_load dsbr_ct_ratio dsbr_del4_load dsbr_del4_ratio dsbr_delsv_load dsbr_delsv_ratio dsbr_dup_load dsbr_sv_load dsbr_first_hit dsbr_second_hit
	#     dsbr_snv_load_cut dsbr_ct_ratio_cut dsbr_del4_load_cut dsbr_del4_ratio_cut dsbr_delsv_load_cut dsbr_delsv_ratio_cut dsbr_dup_load_cut dsbr_sv_load_cut dsbr_first_hit dsbr_second_hit dsbr_score
	#     mmr_snv_load mmr_indel_load mmr_first_hit mmr_second_hit mmr_snv_load_cut mmr_indel_load_cut mmr_score
	
	#print Dumper $data;
	$data->{dsbr_score} = 0;
	$data->{dsbr_snv_load} = $data->{snv_count};
	$data->{dsbr_ct_ratio} = ($data->{snv_count} > 0 ? ($data->{snv_ct} / $data->{snv_count}) : 0);
	$data->{dsbr_del4_load} = $data->{del_4_count};
	if ($data{del_count} == 0) {
		$data->{dsbr_del4_ratio} = "NA";
	} else {
		$data->{dsbr_del4_ratio} = $data->{del_4_count} / $data->{del_count};
	}
	$data->{dsbr_delsv_load} = $data->{'DEL:100'} + $data->{'DEL:1k'} + $data->{'DEL:10k'};
	if ((exists $data{sv_count}) and ($data{sv_count} > 0))
	{
		$data->{dsbr_delsv_ratio} = ($data->{'DEL:100'} + $data->{'DEL:1k'} + $data->{'DEL:10k'}) / $data->{sv_count};
	}
	else
	{
		$data->{dsbr_delsv_ratio} = "NA";
	}
	$data->{dsbr_dup_load} = $data->{'DUP:10k'} + $data->{'DUP:100k'};
	$data->{dsbr_sv_load} = $data->{sv_count};
	$data->{dsbr_first_genes} = "";
	$data->{dsbr_second_genes} = "";

	$data->{mmr_score} = 0;
	$data->{mmr_snv_load} = $data->{snv_count};
	$data->{mmr_indel_load} = $data->{indel_count};
	$data->{mmr_first_genes} = "";
	$data->{mmr_second_genes} = "";

	for my $type (qw/dsbr_snv_load dsbr_del4_load dsbr_del4_ratio dsbr_delsv_load dsbr_delsv_ratio dsbr_dup_load dsbr_sv_load/)
	{
		if ($data->{$type} > $data->{"${type}_cut"})
		{
			$data->{dsbr_score}++;
		}
	}
	for my $type (qw/dsbr_ct_ratio/)
	{
		if ($data->{$type} < $data->{"${type}_cut"})
		{
			$data->{dsbr_score}++;
		}
	}

	for my $type (qw/mmr_snv_load mmr_indel_load/)
	{
		if ($data->{$type} > $data->{"${type}_cut"})
		{
			$data->{mmr_score}++;
		}
	}

	my %geneHits;
	my %geneHitTypes;

	my %somaticHits = (
		"nonsynonymous" => "sNS",
		"stopgain" => "sSG",
		"stoploss" => "sSL",
		"splicing" => "sSP",
		"nonframeshift" => "sNF",
		"frameshift" => "sFS",
		"deletion breakpoint" => "sSV",
		"duplication breakpoint" => "sSV",
		"inversion breakpoint" => "sSV",
		"translocation breakpoint" => "sSV",
		"homozygous deletion" => "sHD",
	);

	my %germlineHits = (
		"stopgain" => "gSG",
		"frameshift" => "gFS",
	);

	my @allGenes;
	push(@allGenes, @dsbrGenes);
	push(@allGenes, @mmrGenes);

	for my $g (@allGenes)
	{
		$geneHits{$g} = 0;
		$geneHitTypes{$g} = "";
		
		for my $v (keys %{ $vars->{$g}{variants} })
		{
			if ($vars->{$g}{variants}{$v}{mutation_class} =~ /germline/)
			{
				unless ($vars->{$g}{variants}{$v}{rarity} eq "common")
				{
					if (exists $germlineHits{ $vars->{$g}{variants}{$v}{mutation_type} })
					{
						$geneHits{$g}++;
						$geneHitTypes{$g} .= $germlineHits{ $vars->{$g}{variants}{$v}{mutation_type} } . "|";
					}
					elsif ($vars->{$g}{variants}{$v}{clinvar} =~ /^CLINSIG=pathogenic/)
					{
						$geneHits{$g}++;
						print "pathogenic $g,$v\n";
						$geneHitTypes{$g} .= $germlineHits{ $vars->{$g}{variants}{$v}{mutation_type} } . "|";
					}
				}
			}
			elsif ($vars->{$g}{variants}{$v}{mutation_class} =~ /somatic/)
			{
				if (exists $somaticHits{ $vars->{$g}{variants}{$v}{mutation_type} })
				{
					$geneHits{$g}++;
					$geneHitTypes{$g} .= $somaticHits{ $vars->{$g}{variants}{$v}{mutation_type} } . "|";
				}
			}
		}

		# loh
		if ($vars->{$g}{ab_counts} =~ /^0\./)
		{
			$geneHits{$g}++;
			$geneHitTypes{$g} .= "sLH|";
		}

		$geneHitTypes{$g} =~ s/\|$//;
	}

	my @hits;
	my $hitString;

	for my $g (@{ $dsbrGenes })
	{
		@hits = split(/\|/, $geneHitTypes{$g});
		if ($geneHits{$g} > 0)
		{
			$data->{dsbr_first_genes} .= "$g($hits[0])|";
		}
		if ($geneHits{$g} > 1)
		{
			$hitString = "";
			for (my $i = 1; $i < scalar(@hits); $i++)
			{
				$hitString .= "$hits[$i];";
			}
			$hitString =~ s/;$//;
			$data->{dsbr_second_genes} .= "$g($hitString)|";
		}
	}

	$data->{dsbr_first_genes} =~ s/\|$//;
	$data->{dsbr_second_genes} =~ s/\|$//;

	if ($data->{dsbr_first_genes} ne "")
	{
		$data->{dsbr_score}++;
	}
	if ($data->{dsbr_second_genes} ne "")
	{
		$data->{dsbr_score}++;
	}

	for my $g (@{ $mmrGenes })
	{
		@hits = split(/\|/, $geneHitTypes{$g});
		if ($geneHits{$g} > 0)
		{
			$data->{mmr_first_genes} .= "$g($hits[0])|";
		}
		if ($geneHits{$g} > 1)
		{
			$hitString = "";
			for (my $i = 1; $i < scalar(@hits); $i++)
			{
				$hitString .= "$hits[$i];";
			}
			$hitString =~ s/;$//;
			$data->{mmr_second_genes} .= "$g($hitString)|";
		}
	}

	$data->{mmr_first_genes} =~ s/\|$//;
	$data->{mmr_second_genes} =~ s/\|$//;

	if ($data->{mmr_first_genes} ne "")
	{
		$data->{mmr_score}++;
	}
	if ($data->{mmr_second_genes} ne "")
	{
		$data->{mmr_score}++;
	}


	print "finished hallmark score\n";

}



sub parseManualDeficiency
{
	my $data = shift;
	my $file = shift;
	my $type = shift;

	my $l;
	my ($donor, $tumour, $normal, $external, $first, $second);

	open (FILE, $file) or die "Couldn't open $file\n";

	$l = <FILE>;		# header
	while ($l = <FILE>)
	{
		chomp $l;
		($donor, $tumour, $normal, $external, $first, $second) = split(/\t/, $l);

		if ($tumour eq $sample)
		{
			$data->{"${type}_deficient"} = "${type}_deficient";
			$data->{"${type}_first_hit"} = $first;
			$data->{"${type}_second_hit"} = $second;
		}
	}

}




# getBase returns the base at a specific position in the reference.  It will initialize fasta handles and pull new chunks of reference if necessary
# input is chromosome, position, reference hash and fasta handles
# output is a single base (and the reference hash and fasta handles may be modified)
sub getBase
{
    my $chr = $_[0];
    my $pos = $_[1];
    my $reference = $_[2];
    my $fastaHandles = $_[3];

    my $chunkStart = int(($pos - 1) / $reference->{"chunkSize"}) * $reference->{"chunkSize"} + 1;       # +1 because the first base in the reference is 1, $pos - 1 so that multiples of chunk size resolve to the correct chunk
    my $chunkEnd = $chunkStart + $reference->{"chunkSize"} - 1;

    unless (exists $reference->{$chr}{$chunkStart}{$pos})       # if the position isn't in our hash, we need to get a new chunk from the reference
    {
        unless (exists ($fastaHandles->{$chr}))     # create a handle for the chromosome fasta, if it doesn't exist
        {
            $fastaHandles->{$chr} = Bio::DB::Fasta->new("$fastaHandles->{path}/$chr.fa");
        }

	my $newChunk;
	$newChunk = uc($fastaHandles->{$chr}->seq($chr, $chunkStart, $chunkEnd));
        
	my $i = $chunkStart;
        for my $base (split("", $newChunk))
        {
            $reference->{$chr}{$chunkStart}{$i} = $base;
            $i++;
        }
    }
	if (exists $reference->{$chr}{$chunkStart}{$pos})
	{
		return $reference->{$chr}{$chunkStart}{$pos};
	}
	else
	{
		return "N";
	}
}

# getRange returns a string of bases from the reference in the specified range by calling getBase
# input is chromosome, start pos, end pos, reference hash and fasta handles
# output is a string of bases
sub getRange
{
    my $chr = $_[0];
    my $start = $_[1];
    my $end = $_[2];
    my $reference = $_[3];
    my $fastaHandles = $_[4];

    my $seq = "";

    for (my $p = $start; $p <= $end; $p++)
    {
        $seq .= getBase($chr, $p, $reference, $fastaHandles);
    }

    return $seq;
}

# freeMemChunks tests if the next indel to be processed is on a different chromosome or more than a chunk away from the reference sequences currently in memory
#   if there exist chunks that we won't need again (assuming the input is sorted) the chunks will be removed from the reference hash
# input is the chromosome and position of the current indel, and the reference hash
# there is no output
sub freeMemChunks
{
    my $chr = $_[0];
    my $pos = $_[1];
    my $reference = $_[2];

    # delete chunks from non-current chromosomes
    for my $refChr (keys %$reference)
    {
        if (($refChr ne $chr) and ($refChr ne "chunkSize"))
        {
            delete $reference->{$refChr};
        }
    }

    # delete chunks if they are more than 1.5 chunks away from the current indel
    # 1.5 so that we are at least in the middle of the current chunk before dropping the previous one
    for my $chunkPos (keys %{ $reference->{$chr} })
    {
        if ($chunkPos < ($pos - (1.5 * $reference->{"chunkSize"})))
        {
            delete $reference->{$chr}{$chunkPos};
        }
    }

    return;
}




sub parseANNOVAR 
{
	my $info = shift;
	my $baseChange = shift;
	my %infoHash;
	my $gene = ".";
	my $consequence = ".";
	my $mutType;
	my ($t, $g, $nucContext, $aaContext);
	my (@types, @genes);
	my %annoHash;

	if ($info =~ /intergenic/i or $info =~ /ANNOVAR=(upstream|downstream)/) {
		#next; #ignore
	} elsif ($info =~ /ANNOVAR=([a-zA-Z0-9]+\|?[a-z]*\|?)\|([a-zA-Z0-9-|()\+\->_:\.]*)/) {
#splicing|NOC2L(NM_015658:exon8:c.888+3T>G)
##TODO: exonic single gene is working!!
##next: add splicing should allow for brackets is working!!
##then check for multi ANNOVAR, add types and genes and for these, should somehow add exonic annotation from ANNOVAR_EXONIC in as well...
#ANNOVAR_EXONIC=synonymous-SNV|OR2T27:NM_001001824:exon1:c.C60T:p.N20N|;ANNOVAR=exonic|OR2T27;superdup
#ANNOVAR_EXONIC=stoploss|TGOLN2:NM_001368095:exon4:c.T1363C:p.X455Q|;
#ANNOVAR=exonic|splicing|TGOLN2|TGOLN2(NM_006464:exon4:c.1309-4T>C|NM_001206844:exon5:c.1135-4T>C)
#rs4240199
		@types = split(/\|/, $1);
		my $part2 = $2;


		if ( all {$_ =~ /UTR|intergenic|intronic/ } @types ) { ##AMY
			return(%infoHash);
		}
		elsif (scalar @types == 1) {
			my $gene = $part2;
			my @variants = split /\|/, $part2;
			if (all {$_ =~ /UTR|intergenic|intronic/ } @variants ) {
				return(%infoHash);
			} elsif ($types[0] eq 'exonic') {
			   @genes = split /\|/, $gene;
			} elsif ($types[0] eq 'splicing') {
			   @genes = split /\)/, $gene;
			   for (my $i = 0; $i < scalar @genes; $i++) {
				$genes[$i] =~ s/^\|//;
				if ($genes[$i] =~ /\(/) {
				   $genes[$i] .= ")";
				}
			   }
			} else {
				@genes = ($gene);
			}
		} elsif (scalar @types >= 2) {
		     if ($part2 =~ /(([a-zA-Z0-9\-]+\|)+)([a-zA-Z0-9>\+\-_\.:|()]+)$/) {
			@genes = split /\|/, $1;
			my @rest = split /\)/, $3;
			for (my $i = 0; $i < scalar @rest; $i++) {
                                $rest[$i] =~ s/^\|//;
                                if ($rest[$i] =~ /\(/) {
                                   $rest[$i] .= ")";
                                }
                        }
			@genes = (@genes, @rest);
		     } else {
			warn "parse error: @types, $part2, $info\n";
			die;
		     }

		}
		for (my $i = 0; $i < scalar @types; $i++)
		{
			$t = $types[$i];
			if ($t eq "exonic")
			{
				for $g (split(/\|/, $genes[$i]))
				{
					unless ($g eq 'splicing' or $g =~ /\(/) {
					   $annoHash{$t}{$g}++;
					}
				}
			}
			elsif ($t eq "splicing")
			{
				unless (defined $genes[$i] ) {
					warn "@types, @genes, $info\n";
					die;
				}
				if ($genes[$i] =~ /^(.*?)\((.*?)\)/) ##may not capture the whole () expression if multiple contexts split by |
				{
					my $gene = $1;
					my $context = $2;
					$annoHash{$t}{$gene} = $context;
					$genes[$i] =~ s/^.*?\(.*?\)//;
					$genes[$i] =~ s/^,//;
				} elsif (defined $baseChange) {
					$annoHash{$t}{$genes[$i]} = "$baseChange";
				}
			}
		}


		for $t (sort keys %annoHash)
		{
			if ($t eq "splicing")
			{
				for $g (sort keys %{ $annoHash{$t} })
				{
					$nucContext = "";
					for my $isoform (split (/\|/, $annoHash{$t}{$g}))
					{
						if ($isoform =~ /^(.*?):.*?:(.*?)$/)
						{
							$nucContext .= "|$1:$2";
						}
					}
					$nucContext =~ s/^\|//;
					$aaContext = "NA";

					$infoHash{$t}{$g}{consequence} = "splicing";
					$infoHash{$t}{$g}{nuc} = $nucContext;
					$infoHash{$t}{$g}{aa} = $aaContext;

									
				
				}
			}
			elsif ($t eq "exonic")
			{
				for $g (sort keys %{ $annoHash{$t} })
				{
					$nucContext = "";
					$aaContext = "";
					$consequence = ".";
					#nonframeshift-deletion_LINGO4:NM_001004432:exon2:c.404_424del:p.135_142del_;
					if ($info =~ /ANNOVAR_EXONIC=([a-z]+-?([a-zA-Z]+)?)\|(.*:NM_[0-9]+:.*);ANNOVAR/)
					{
						$consequence = "$1,$3";		# seems like there is never more than one kind of exonic consequence for a variant, even when multiple genes are affected?
						$consequence =~ s/\|$//;
						for my $isoform (split (/\|/, $consequence))
						{
							if ($isoform =~ /^([a-z]+-?[a-zA-Z]+,)?$g:(NM_[0-9]+):.*:(.*?):(.*?)\|?$/)
							{
								$nucContext .= "|$2:$3";
								$aaContext .= "|$2:$4";
							} elsif ($isoform =~ /^([a-z]+-?[a-zA-Z]+,)?$g:(NM_[0-9]+):wholegene/) {
								$nucContext .= "|$2:wholegene";
								$aaContext .= "|$2:wholegene";
							} else {
								my @othergenes = grep !/^$g$/, @genes;
								if (scalar @othergenes > 0) {
#									warn "$g, @othergenes\n";
								   for (my $i = 0; $i < scalar @othergenes; $i++) {
									my $gene = $othergenes[$i];
									if ($isoform =~ /^([a-z]+-?[a-zA-Z]+,)?$gene:(NM_[0-9]+):.*:(.*?):(.*?)(\|)?$/) {
										$nucContext .= "|$2:$3";
										$aaContext .= "|$2:$4";
									} elsif ($isoform =~ /^([a-z]+-?[a-zA-Z]+,)?$gene:(NM_[0-9]+):wholegene/) {
										$nucContext .= "|$2:wholegene";
										$aaContext .= "|$2:wholegene";
									}
								   }
								   if ($nucContext !~ /[a-zA-Z0-9]/) {
									warn "isoform mismatch1: $isoform\t$consequence\t$info\n";
									warn "$g: @othergenes, types: @types\n";
									for (my $i = 0; $i < scalar @othergenes; $i++) {
                                      	my $gene = $othergenes[$i];
										if ($isoform =~ /^([a-z]+-?[a-zA-Z]+,)?$gene:(NM_[0-9]+):.*:(.*?):(.*?)(\|)?$/) {
                                        	$nucContext .= "|$2:$3";
                                        	$aaContext .= "|$2:$4";
                                        } else {
											warn "$gene not the right one\n";
										}
									}
									die;
								   }
								} else {
								print "isoform mismatch2: $isoform, $consequence, $info\n";
								print "$g: @genes, types: @types\n";
								}
					
								
							}
						}
						$nucContext =~ s/^\|//;
						$aaContext =~ s/^\|//;
						
						$infoHash{$t}{$g}{consequence} = $consequence;
						$infoHash{$t}{$g}{nuc} = $nucContext;
						$infoHash{$t}{$g}{aa} = $aaContext;
					}
					else
					{
						%infoHash = ();
					}
				}
			}
		}


	}
	else
	{
		if ($info =~ /ANNOVAR=/ and $info !~ /ncRNA/ and $info !~/intronic/) {
			print "gene mismatch: $info\n";
		}
		%infoHash = ();
	}
		
	return %infoHash;

}
