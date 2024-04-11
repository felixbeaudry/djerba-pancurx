#!/usr/bin/perl

use strict;
use warnings;

use JSON;

use lib '/.mounts/labs/PCSI/users/rdenroche/minerva/';
use minerva;
use Data::Dumper;


my $variantFile = $ARGV[0];
my $expressionFile = $ARGV[1];

my $reportPath = set_out_path();

#my $geneKB = "1yNkkcMyaSz5pabtiJOMbz2w7QmCDa_JVyrGOUJZquyo";
#my $varKB = "1htEqHqvrS19v7S1pjOxam80z2uH44HCZgsXvtL5yFfc";
#my $patientKB = "12Iq61c_dw4083hvx3OtNxN1O2SvPf9TbcSMwH81UMs0";
my $geneKB = "11Zf6gRbVaNyiD88iMeUC4RxVIsNTcIM8_ZFBeW2L64g";
my $varKB = "1J9ap1R4Mtxdof9DIHs2pNlkNcSFIMl_K6jAQQ5DhEuI";
my $patientKB = "1lX0J9afUQbeZoBhoSlHbDOzQSAEbbvZJPyQTKDeQRm0";

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon += 1;
$mon = sprintf("%02d", $mon);
$mday = sprintf("%02d", $mday);
$hour = sprintf("%02d", $hour);
$min = sprintf("%02d", $min);
$sec = sprintf("%02d", $sec);

my $date = "$year${mon}${mday}_${hour}h${min}m${sec}s";

my $version = "v0.0";

my %includeGenes;
my %excludeGenes;


my %data;

my ($donor, $sample, $comp);

my ($l, $uniq);
my (@f, @lines, @header);
my (%row, %geneKB, %varKB);

# read variant.csv into gene structure


open (FILE, $variantFile) or die "Couldn't open $variantFile\n";

while ($l = <FILE>)
{
	chomp $l;
	if ($l =~ /^donor/)
	{
		@header = split (/,/, $l);
	}
	else
	{
		%row = ();
		@f = split(/,/, $l);
		for (my $i = 0; $i < scalar(@f); $i++)
		{
			$row{$header[$i]} = $f[$i];
		}

		for my $h (qw/copy_number ab_counts gene_position maf_mean maf_p full_name gene_chr cosmic_census_flag cosmic_census_data/)
		{
			$data{$row{gene}}{$h} = $row{$h};
		}
		unless (($row{mutation_type} eq "altered promoter") or ($row{mutation_type} eq "NA"))
		{
			for my $h (qw/mutation_class mutation_type position fusion_genes base_change tumour_freq tumour_depth normal_freq normal_depth nuc_context aa_context dbsnp cosmic cadd_phred rarity 1000G_all ExAC_all ESP6500siv2_all clinvar/)
			{
				$data{$row{gene}}{variants}{"$row{mutation_type},$row{position},$row{base_change}"}{$h} = $row{$h};
			}
		}
	}
}
close FILE;

# hope the last row isn't messed up
$donor = $row{donor};
$sample = $row{tumour};





# read RNA expressions & ranks?

if ((defined $expressionFile) and (-e $expressionFile))
{
	# Gene.ID Gene.Name       Reference       Strand  Start   End     Coverage        FPKM    TPM     foldchange_cohort
	
	open (FILE, $expressionFile) or die "Couldn't open $expressionFile\n";

	while ($l = <FILE>)
	{
		chomp $l;
		%row = ();

		if ($l =~ /^Gene\.ID\t/)
		{
			@header = split(/\t/, $l);
		}
		else
		{
			@f = split(/\t/, $l);
			for (my $i = 0; $i < scalar(@f); $i++)
			{
				$row{$header[$i]} = $f[$i];
			}

			if (exists $data{ $row{"Gene.Name"} })
			{
				$data{$row{"Gene.Name"}}{rna}{fpkm} = $row{FPKM};
				$data{$row{"Gene.Name"}}{rna}{tpm} = $row{TPM};
				if (exists $row{foldchange_cohort})
				{
					$data{$row{"Gene.Name"}}{rna}{foldchange} = $row{foldchange_cohort};
					$data{$row{"Gene.Name"}}{rna}{cohort} = "cohort";
				}
				elsif (exists $row{foldchange_compass})
				{
					$data{$row{"Gene.Name"}}{rna}{foldchange} = $row{foldchange_compass};
					$data{$row{"Gene.Name"}}{rna}{cohort} = "compass";
				}
			}
			else
			{
				warn "no dna match for rna gene " . $row{"Gene.Name"} . "\n";
			}
		}
	}

	close FILE;



#	# old fpkm code below
#	# tracking_id     class_code      nearest_ref_id  gene_id gene_short_name tss_id  locus   length  coverage        FPKM    FPKM_conf_lo    FPKM_conf_hi    FPKM_status
#	
#	open (FILE, $expressionFile) or die "Couldn't open $expressionFile\n";
#
#	while ($l = <FILE>)
#	{
#		chomp $l;
#		%row = ();
#
#		if ($l =~ /^tracking_id\t/)
#		{
#			@header = split(/\t/, $l);
#		}
#		else
#		{
#			@f = split(/\t/, $l);
#			for (my $i = 0; $i < scalar(@f); $i++)
#			{
#				$row{$header[$i]} = $f[$i];
#			}
#
#			if (exists $data{ $row{gene_short_name} })
#			{
#				$data{$row{gene_short_name}}{rna}{fpkm} = $row{FPKM};
#				$data{$row{gene_short_name}}{rna}{status} = $row{FPKM_status};
#				$data{$row{gene_short_name}}{rna}{fpkm_conf_lo} = $row{FPKM_conf_lo};
#				$data{$row{gene_short_name}}{rna}{fpkm_conf_hi} = $row{FPKM_conf_hi};
#			}
#			else
#			{
#				warn "no dna match for rna gene $row{gene_short_name}\n";
#			}
#
#
#		}
#	}
#	close FILE;
}

# read RNA fusions?





# read gene KB
warn "\nReading COMPASS Gene Knowledgebase\n";
#my $doc = `curl "https://docs.google.com/spreadsheets/d/$geneKB/export?format=tsv" | tr '\r' '\n' | grep "	"`;
my $doc = `curl https://docs.google.com/spreadsheets/d/$geneKB/gviz/tq?tqx=out:csv | sed 's/^"//' | sed 's/","/	/g' | sed 's/"\$//' | tr '\r' '\n' | grep "	"`;

chomp $doc;
@lines = split(/\n/, $doc);

for $l (@lines)
{
	%row = ();
	if ($l =~ /^gene/)
	{
		@header = split(/\t/, $l);
	}
	else
	{
		@f = split(/\t/, $l);
		$row{gene} = "";
		$row{type} = "";
		$row{source} = "";
		for (my $i = 0; $i < scalar(@f); $i++)
		{
			$row{$header[$i]} = $f[$i];
		}

		$uniq = "$row{gene},$row{type},$row{source}";

		for my $h (qw/type preferred_isoform oncogene supressor report_text source/)
		{
			$geneKB{$row{gene}}{$uniq}{$h} = $row{$h};
		}
	}
	if (defined $row{gene} and $row{gene} eq "KRAS")
	{
#		print "$l\n";
	}
}

# read variant KB
warn "\nReading COMPASS Variant Knowledgebase\n";

#$doc = `curl "https://docs.google.com/spreadsheets/d/$varKB/export?format=tsv" | tr '\r' '\n' | grep "	"`;
$doc = `curl https://docs.google.com/spreadsheets/d/$varKB/gviz/tq?tqx=out:csv | sed 's/^"//' | sed 's/","/	/g' | sed 's/"\$//' | tr '\r' '\n' | grep "	"`;
chomp $doc;
@lines = split(/\n/, $doc);

for $l (@lines)
{
	%row = ();
	if ($l =~ /^gene/)
	{
		@header = split(/\t/, $l);
	}
	else
	{
		@f = split(/\t/, $l);
		for (my $i = 0; $i < scalar(@f); $i++)
		{
			$row{$header[$i]} = $f[$i];
		}

		$uniq = "$row{gene},$row{alteration},$row{source}";
		for my $h (qw/refseq alteration protein_change oncogenicity effect pmids abstracts source report_text/)
		{
			$varKB{$row{gene}}{$uniq}{$h} = $row{$h};
		}
	}
}


# read patient KB (for gene lists)
warn "\nReading COMPASS Patient Knowledgebase\n";
#$doc = `curl "https://docs.google.com/spreadsheets/d/$patientKB/export?format=tsv" | tr '\r' '\n' | grep "	"`;
$doc = `curl https://docs.google.com/spreadsheets/d/$patientKB/gviz/tq?tqx=out:csv | sed 's/^"//' | sed 's/","/	/g' | sed 's/"\$//' | tr '\r' '\n' | grep "	"`;

chomp $doc;
@lines = split(/\n/, $doc);

for $l (@lines)
{
	%row = ();
	if ($l =~ /^pcsi_id/)
	{
		@header = split(/\t/, $l);
	}
	else
	{
		@f = split(/\t/, $l);
		for (my $i = 0; $i < scalar(@f); $i++)
		{
			$row{$header[$i]} = $f[$i];
		}

		if ((defined $row{sample_id}) and ($row{sample_id} eq $sample))
		{
			# grab include and exclude genes
			$row{exclude_genes} =~ s/ //g;
			for my $g (split(/,/, $row{exclude_genes}))
			{
				$excludeGenes{$g}++;
			}

			$row{include_genes} =~ s/ //g;
			for my $g (split(/,/, $row{include_genes}))
			{
				$includeGenes{$g}++;
			}
		}
	}
}

# add geneKB and varKB info to genes in %data

for my $g (keys %data)
{
	if (exists $geneKB{$g})
	{
		# key is gene,type,source
		for my $u (keys %{ $geneKB{$g} })
		{
			if (exists $geneKB{$g}{$u}{preferred_isoform}) {
				$data{$g}{genekb}{ $geneKB{$g}{$u}{type} }{ $geneKB{$g}{$u}{source} }{preferred_isoform} = $geneKB{$g}{$u}{preferred_isoform};
			}
			$data{$g}{genekb}{ $geneKB{$g}{$u}{type} }{ $geneKB{$g}{$u}{source} }{oncogene} = $geneKB{$g}{$u}{oncogene};
			$data{$g}{genekb}{ $geneKB{$g}{$u}{type} }{ $geneKB{$g}{$u}{source} }{supressor} = $geneKB{$g}{$u}{supressor};
			$data{$g}{genekb}{ $geneKB{$g}{$u}{type} }{ $geneKB{$g}{$u}{source} }{report_text} = $geneKB{$g}{$u}{report_text};

			if ($g eq "KRAS")
			{
#				print "$geneKB{$g}{$u}{preferred_isoform}, $geneKB{$g}{$u}{oncogene}, $geneKB{$g}{$u}{supressor}, $geneKB{$g}{$u}{report_text}\n";
			}
		}
	}
}

my $aa;
my %aaHash;
my $matched;

my %truncating = (
	"stopgain" => 1,
	"splicing" => 1,
	"frameshift" => 1,
	"deletion breakpoint" => 1,
	"duplication breakpoint" => 1,
	"inversion breakpoint" => 1,
	"translocation breakpoint" => 1,
);

for my $g (keys %data)
{
	if (exists $varKB{$g})
	{
		for my $v (keys %{ $data{$g}{variants} })
		{
			if (($data{$g}{variants}{$v}{mutation_class} =~ /somatic/) or (($data{$g}{variants}{$v}{mutation_class} =~ /germline/) and (exists $data{$g}{genekb}{germline})))
			{
				%aaHash = ();
				for $aa (split(/\|/, $data{$g}{variants}{$v}{aa_context}))
				{
					$aa =~ s/.*\.//;
					$aaHash{$aa}++;
				}

				for my $u (keys %{ $varKB{$g} })
				{
					$matched = 0;

					# match splice - not yet implemented!!!
					

					# match aa change
					unless ($data{$g}{variants}{$v}{aa_context} eq "NA")
					{
						if (exists $aaHash{ $varKB{$g}{$u}{alteration} })
						{
							$matched = 1;
						}
					}
					
					# match truncating
					if ($varKB{$g}{$u}{alteration} eq "Truncating Mutations")
					{
						if (exists $truncating{ $data{$g}{variants}{$v}{mutation_type} })
						{
							$matched = 1;
						}
					}

					# match amplification
					if ($varKB{$g}{$u}{alteration} eq "Amplification")
					{
						if ($data{$g}{variants}{$v}{mutation_type} eq "strong amplification")
						{
							$matched = 1;
						}
					}

					# match deletion
					if ($varKB{$g}{$u}{alteration} eq "Deletion")
					{
						if ($data{$g}{variants}{$v}{mutation_type} eq "homozygous deletion")
						{
							$matched = 1;
						}
					}



					if ($matched > 0)
					{
						for my $f (qw/alteration oncogenicity effect pmids abstracts report_text/)
						{
							$data{$g}{variants}{$v}{varkb}{ $varKB{$g}{$u}{source} }{$f} = $varKB{$g}{$u}{$f};
						}
					}
				}
			}
		}
	}
}


# determine which genes are activated or inactivated
my ($lofCount, $gofCount, $loh, $amped, $ab);
my $germLof;
my @abs;
my ($kbLoss,$kbGain);

for my $g (keys %data)
{
	$lofCount = 0;
	$gofCount = 0;
	$germLof = 0;
	$amped = 0;
	$loh = 0;

	if (defined $data{$g}{ab_counts})
	{
		@abs = split(/\|/, $data{$g}{ab_counts});
		for $ab (@abs)
		{
			if ($ab =~ /^0/)
			{
				$loh = 1;
			}
		}
	}

	$lofCount += $loh;

	for my $v (keys %{ $data{$g}{variants} })
	{
		if ($data{$g}{variants}{$v}{mutation_class} =~ /somatic/)
		{
			if ($data{$g}{variants}{$v}{mutation_type} eq "homozygous deletion")
			{
				$lofCount += 2;
			}
			# check for a varkb entry
			elsif (exists $truncating{ $data{$g}{variants}{$v}{mutation_type} })
			{
				$lofCount += 1;
			}
			elsif ($data{$g}{variants}{$v}{mutation_type} eq "strong amplification")
			{
				$amped = 1;
			}
			elsif (exists $truncating{ $data{$g}{variants}{$v}{mutation_type} })
			{
				if (exists $data{$g}{variants}{$v}{varkb})
				{
					$kbLoss = 0;
					$kbGain = 0;
					for my $kb (keys %{ $data{$g}{variants}{$v}{varkb} })
					{
						if (($data{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /Loss-of-function/) or ($data{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /Pathogenic/))
						{
							$kbLoss = 1;
						}
					}
					$lofCount += $kbLoss;
				}
				else	# assume LOF for truncating if not in KB
				{
					$lofCount += 1;
				}

				# handle fusion partners here??
			}
			else
			{
				if (exists $data{$g}{variants}{$v}{varkb})
				{
					$kbLoss = 0;
					$kbGain = 0;
					
					for my $kb (keys %{ $data{$g}{variants}{$v}{varkb} })
					{
						if (($data{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /Loss-of-function/) or ($data{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /Pathogenic/))
						{
							$kbLoss = 1;
						}
						elsif ($data{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /Gain-of-function/)
						{
							$kbGain = 1;
						}
					}
					$lofCount += $kbLoss;
					$gofCount += $kbGain;
				}
				else
				{
					# do what with missense not in kb??
					# also, skip non-coding...
				}
			}

		}
		elsif ($data{$g}{variants}{$v}{mutation_class} =~ /germline/)
		{
			if (exists  $data{$g}{genekb}{germline})
			{
				# check for a varkb entry
				if (exists $data{$g}{variants}{$v}{varkb})
				{
					for my $kb (keys %{ $data{$g}{variants}{$v}{varkb} })
					{
						if (($data{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /Loss-of-function/) or ($data{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /Pathogenic/))
						{
							$germLof = 1;
						}
					}
				}

				# check if germline variant has been lost due to somatic LOH
				if (($loh == 1) and ($data{$g}{variants}{$v}{tumour_freq} ne "NA"))
				{
					if ($data{$g}{variants}{$v}{tumour_freq} < 0.4)
					{
						$germLof = 0;
						$data{$g}{variants}{$v}{lost_in_tumour} = 1;
					}
				}
			}
		}
	}
	$lofCount += $germLof;


	if ($lofCount >= 2)
	{
		$data{$g}{biallelic_lof} = 1;
	}
	if ($gofCount > 0)
	{
		$data{$g}{gof} = 1;
	}




}


# print dated & versioned gene structured json, including coding and protien change, sequencing metrics on variants, symbols, isoforms, copy number, expression & series rank (resected & advanced), etc...

my $outPath = "$reportPath/$donor/$sample/";

if (!-e $outPath)
{
	`mkdir $outPath/archive -p`;
}

my $outFile = "$outPath/$sample-genes.json";
my $outFileArchive = "$outPath/archive/$sample-genes_${version}_${date}.json";

open (my $fh, ">$outFile") or die "Couldn't open >$outFile\n";
open (my $fha, ">$outFileArchive") or die "Couldn't open >$outFileArchive\n";

my $json = JSON->new->allow_nonref;

my $string = 'translocation breakpoint,chr5:81108716-chr4:71561531,NA';
print "$data{RASGRF2}{variants}{$string}{mutation_type}\n";
#die;
my $pretty = $json->pretty->encode(\%data);

print $fh $pretty;
print $fha $pretty;

close $fh;
close $fha;

