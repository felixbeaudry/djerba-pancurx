#!/usr/bin/perl

package minerva;

# use lib '/.mounts/labs/PCSI/users/rdenroche/minerva/';
# use minerva;



require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/set_out_path set_common_path get_sample_kb get_clinical_kb read_meta read_genes get_gene_list short_gene short_legend full_gene classify_metric get_best_isoform
 convert_lims_tissue print_table pass_icon fail_icon pass_icon_size fail_icon_size pass_or_fail_icon pass_or_fail_icon_size print_text_box commify get_all_genes_with_variants gene_has_variant
 ssm_variant_row get_rna_mark_class get_cn_mark_class slide_gene variant_to_report print_mof_gene_box full_gene_slide get_preferred_isoform short_legend_page full_gene_page variant_row_tables/;

use strict;
use warnings;
use JSON;




sub set_out_path
{
	return "/.mounts/labs/PCSI/reports/minerva/";
}

sub set_common_path
{
	return "../../common/";
}

sub get_sample_kb
{
	#return "12Iq61c_dw4083hvx3OtNxN1O2SvPf9TbcSMwH81UMs0";
	return "1lX0J9afUQbeZoBhoSlHbDOzQSAEbbvZJPyQTKDeQRm0";
}

sub get_clinical_kb
{
#	return "1Uv4NppBhLG7O8llNg6uiy-ebhcE1AzzUtA-yO6121yA";
	return "1oSbZcgLz-Z-b4f3QdROnNXJFEAQ48tzy5PVJXgCIPfc";
}

sub read_meta
{
	my $file = shift @_;

	my $l;
	my (@header, @f);
	my %data;

	open (FILE, $file) or die "Couldn't open $file\n";

	while ($l = <FILE>)
	{
		chomp $l;

		if ($l =~ /^donor/)
		{
			@header = split(/,/, $l);
		}
		else
		{
			@f = split(/,/, $l);

			for (my $i = 0; $i < scalar(@f); $i++)
			{
				$data{$header[$i]} = $f[$i];
			}
		}
	}
	close FILE;

	return %data;

}

sub read_genes
{
	my $file = shift @_;

	my $l;
	my $json = "";

	open (FILE, $file) or die "Couldn't open $file\n";

	while ($l = <FILE>)
	{
		chomp $l;
		$json .= $l;
	}
	close FILE;

	my $ref = decode_json($json);
	return %$ref;

}

sub convert_lims_tissue
{
	my $tissue = shift @_;

	my %limsTissue = (
	    "Ad" => "Adipose",
	    "Ap" => "Appendix",
	    "Ag" => "Adrenal Gland",
	    "As" => "Ascites",
	    "Bm" => "Bone Marrow",
	    "Bn" => "Brain",
	    "Br" => "Breast",
	    "Bu" => "Buccal cells",
	    "Cb" => "Cored Blood",
	    "Cn" => "Central Nervous System",
	    "Du" => "Duodenum",
	    "Ep" => "Esophagus",
	    "Es" => "Esophagus",
	    "Fs" => "Foreskin",
	    "Gb" => "Gallbladder",
	    "Hr" => "Heart",
	    "Ki" => "Kidney",
	    "Le" => "Leukocyte",
	    "Li" => "Large Intestine",
	    "Ln" => "Lymph Node",
	    "Lu" => "Lung",
	    "Lv" => "Liver",
	    "Lx" => "Larynx",
	    "Ly" => "Lymphocyte",
	    "Md" => "Mediastinum",
	    "Me" => "Mesenchyme",
	    "Nk" => "Neck",
	    "Oc" => "Oral Cavity",
	    "Om" => "Omentum",
	    "Ov" => "Ovary",
	    "Pa" => "Pancreas",
	    "Pb" => "peripheral blood",
	    "Pr" => "Prostate",
	    "Sa" => "Saliva",
	    "Sg" => "Salivary Gland",
	    "Si" => "Small Intestine",
	    "Sk" => "Skin",
	    "Sm" => "Skeletal Muscle",
	    "Sp" => "Spleen",
	    "St" => "Stomach",
	    "Ta" => "Tail (Referring to Models)",
	    "Tr" => "Trachea",
	    "Mu" => "Muscle",
	    "Wm" => "Worm (Referring to Models)",
	    "nn" => "Unknown"
	);

	if (exists $limsTissue{$tissue})
	{
		return $limsTissue{$tissue};
	}
	else
	{
		return $tissue;
	}
}

sub is_truncating
{
}

sub classify_metric
{
	my $metric = shift;
	my $val = shift;

	if ($metric eq "tmb")
	{
		if ($val > 10)
		{
			return "high";
		}
		elsif ($val > 5)
		{
			return "elevated";
		}
		elsif ($val > 1)
		{
			return "typical";
		}
		else
		{
			return "low";
		}
	}

	if ($metric eq "snv")
	{
		if ($val > 20000)
		{
			return "high";
		}
		elsif ($val > 10000)
		{
			return "elevated";
		}
		elsif ($val > 2500)
		{
			return "typical";
		}
		else
		{
			return "low";
		}
	}

	if ($metric eq "indel")
	{
		if ($val > 2000)
		{
			return "high";
		}
		elsif ($val > 1000)
		{
			return "elevated";
		}
		elsif ($val > 250)
		{
			return "typical";
		}
		else
		{
			return "low";
		}
	}

	if ($metric eq "sv")
	{
		if ($val > 280)
		{
			return "high";
		}
		elsif ($val > 140)
		{
			return "elevated";
		}
		elsif ($val > 35)
		{
			return "typical";
		}
		else
		{
			return "low";
		}
	}

	if ($metric eq "ploidy")
	{
		if ($val > 2.3)
		{
			return "polyploid";
		}
		else
		{
			return "diploid";
		}
	}

	return "metric not implemented";
}


sub filename_date
{
}

sub printing_date
{
}

sub get_mut_letter
{
	my $var = shift;

	my %mutLetter = (
		"nonsynonymous" => "n",
		"stopgain" => "x",
		"stoploss" => "l",
		"frameshift" => "f",
		"nonframeshift" => "i",
		"splicing" => "s",
		"strong amplification" => "a",
		"homozygous deletion" => "d",
		"deletion breakpoint" => "r",
		"duplication breakpoint" => "r",
		"inversion breakpoint" => "r",
		"translocation breakpoint" => "r",
		"deletion potential fusion" => "r",
		"duplication potential fusion" => "r",
		"inversion potential fusion" => "r",
		"translocation potential fusion" => "r",
		"unknown" => "u",	# shrug
    );

	if (exists $mutLetter{$var})
	{
		return $mutLetter{$var};
	}
	else
	{
		warn "No letter match for $var\n";
		return $var;
	}

}

sub get_rna_mark_class
{
	my $fold = shift;

	my $rnaClass = "rna0";

	if (defined $fold)
	{
		if ($fold eq "-Inf")
		{
		$rnaClass = "rna-6";
		}
		elsif ($fold < -6)
		{
			$rnaClass = "rna-6";
		}
		elsif ($fold < -4)
		{
			$rnaClass = "rna-4";
		}
		elsif ($fold < -3)
		{
			$rnaClass = "rna-3";
		}
		elsif ($fold < -2)
		{
			$rnaClass = "rna-2";
		}
		elsif ($fold > 6)
		{
			$rnaClass = "rna6";
		}
		elsif ($fold > 4)
		{
			$rnaClass = "rna4";
		}
		elsif ($fold > 3)
		{
			$rnaClass = "rna3";
		}
		elsif ($fold > 2)
		{
			$rnaClass = "rna2";
		}
	}

	return $rnaClass;
}

sub get_cn_mark_class
{
	my $cn = shift;
	my $ploidy = shift;

	my $cnClass = "cn_typical";

	if ($cn < 0.5)
	{
		$cnClass = "cn_loss";
	}
	elsif ($cn < 1.5)
	{
		$cnClass = "cn_loh";
	}
	elsif ($cn > ($ploidy * 3))
	{
		$cnClass = "cn_very_gain";
	}
	elsif ($cn > ($ploidy * 1.5))
	{
		$cnClass = "cn_gain";
	}

	return $cnClass;
}




sub get_all_genes_with_variants
{
	my $genes = shift;

	my $include;
	my @list;

	for my $g (keys %{ $genes })
	{
		$include = 0;
		if (exists $genes->{$g}{variants})
		{
			for my $v (keys %{ $genes->{$g}{variants} })
			{
				if ($genes->{$g}{variants}{$v}{mutation_class} =~ /somatic/)
				{
					$include = 1;
				}
			}
		}

		if ($include == 1)
		{
			push(@list, $g);
		}
	}

	return sort @list;

}




sub short_legend
{
	my $html = "<h4 class=\"legend_title\">Legend:</h4>";
    $html .= "<table style=\"width:100%\">\n";
    $html .= "<tr><td>Nonsynonymous SNV: n</td><td>Splice altering SNV: s</td><td>Rearrangement Breakpoint: r</td><td>Gain-of-function: <mark class=\"gof\"><i>GENE</i></mark></td></tr>\n";
    $html .= "<tr><td>Stopgain SNV: x</td><td>Frameshift Indel: f</td><td>Strong Amplification: a</td><td>Biallelic Inactivation: <mark class=\"biallelic\"><i>GENE</i></mark></td></tr>\n";
    $html .= "<tr><td>Stoploss SNV: l</td><td>In-frame Indel: i</td><td>Homozygous Deletion: d</td><td>Germline Alteration: g_</td></tr>\n";
    $html .= "</table>\n";

	return $html;

}

sub short_legend_page
{
	my $html = "";
    $html .= "<table style=\"width:100%;font-size:0.75em\">\n";
    $html .= "<tr><td>Nonsynonymous SNV: n</td><td>Frameshift Indel: f</td><td>Rearrangement Breakpoint: r</td></tr>\n";
    $html .= "<tr><td>Stopgain SNV: x</td><td>In-frame Indel: i</td><td>Gain-of-function: <mark class=\"gof\"><i>GENE</i></mark></td></tr>\n";
    $html .= "<tr><td>Stoploss SNV: l</td><td>Strong Amplification: a</td><td>Biallelic Inactivation: <mark class=\"biallelic\"><i>GENE</i></mark></td></tr>\n";
    $html .= "<tr><td>Splice altering SNV: s</td><td>Homozygous Deletion: d</td><td>Germline Variant: g_</td></tr>\n";
    $html .= "</table>\n";

	return $html;

}



sub gene_has_variant
{
	my $g = shift;
	my $genes = shift;
	my $spec = shift;	# "no_germ" or "no_som" to exclude germline or somatic variants

	unless (defined $spec)
	{
		$spec = "all";
	}

	my $col = "typical";
	my $report;
	my $germ;
	my $vars = "";

	if (exists $genes->{$g}{biallelic_lof})
	{
		$col = "biallelic";
	}
	elsif (exists $genes->{$g}{gof})
	{
		$col = "gof";
	}

	if (exists $genes->{$g}{variants})
	{
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			$report = 0;
			$germ = 0;
			if ($genes->{$g}{variants}{$v}{mutation_class} =~ /^germline/)
			{
				unless ($spec eq "no_germ")
				{
					# only report germline if they're confirmed in the varkb
					if (exists $genes->{$g}{variants}{$v}{varkb})
					{
						for my $kb (sort keys %{ $genes->{$g}{variants}{$v}{varkb} })
						{
							if ($genes->{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /function/)
							{
								$germ = 1;
								$report = 1;
							}
						}
					}
				}
			}
			elsif ($genes->{$g}{variants}{$v}{mutation_class} =~ /^somatic/)
			{
				unless ($spec eq "no_som")
				{
					$report = 1;
				}
			}

			if ($report == 1)
			{
				return 1;
			}
		}

	}
	return 0;
}

sub ssm_variant_row
{
	my $g = shift;
	my $v = shift;
	my $genes = shift;

	# gene, variant type, position, ref>alt, nuc change, aa change, IDs
	
	my $html = "";

	my ($iso, $nu, $aa);

	($iso,$nu) = get_best_isoform($g,$v,$genes->{$g}{variants}{$v}{nuc_context},$genes);
	($iso,$aa) = get_best_isoform($g,$v,$genes->{$g}{variants}{$v}{aa_context},$genes);

	$html .= "<tr><td><a href=\"#$g\"><i>$g</i></a></td>";
	$html .= "<td>$genes->{$g}{variants}{$v}{mutation_type}</td>";
	$html .= "<td>$genes->{$g}{variants}{$v}{position}</td>";
	$html .= "<td>$genes->{$g}{variants}{$v}{base_change}</td>";
	
	$html .= "<td>$nu</td>";
	$html .= "<td>$aa</td>";

	my $ids = "$iso, ";
	for my $i (qw/dbsnp cosmic/)
	{
		if ((defined $genes->{$g}{variants}{$v}{$i}) and ($genes->{$g}{variants}{$v}{$i} ne "") and ($genes->{$g}{variants}{$v}{$i} ne "."))
		{
			$ids .= "$genes->{$g}{variants}{$v}{$i}, "
		}
	}
	$ids =~ s/, $//;
	
	$html .= "<td>$ids</td></tr>\n";


	return $html;
}


sub sv_variant_row
{
	my $g = shift;
	my $v = shift;
	my $genes = shift;

	# gene, variant type, breakpoints, potential fusion partner
	my $html = "";

	$html .= "<tr><td><a href=\"#$g\"><i>$g</i></a></td>";
	$html .= "<td>$genes->{$g}{variants}{$v}{mutation_type}</td>";
	$html .= "<td>$genes->{$g}{variants}{$v}{position}</td>";
	$html .= "<td>$genes->{$g}{variants}{$v}{fusion_genes}</td>";
	$html .= "</tr>\n";

	return $html;

}

sub cnv_variant_row
{
	my $g = shift;
	my $v = shift;
	my $genes = shift;

	# gene, variant type, copy number, AB count, expression fold change
	
	my $html = "";
	$html .= "<tr><td><a href=\"#$g\"><i>$g</i></a></td>";
	$html .= "<td>$genes->{$g}{variants}{$v}{mutation_type}</td>";
	$html .= "<td>$genes->{$g}{variants}{$v}{position}</td>";
	$html .= "<td>$genes->{$g}{copy_number}</td>";
	$html .= "<td>$genes->{$g}{ab_counts}</td>";
	if ((defined $genes->{$g}{rna}{foldchange}) and ($genes->{$g}{rna}{foldchange} ne "NA"))
	{
		$html .= "<td>$genes->{$g}{rna}{foldchange}</td>";
	}
	else
	{
		$html .= "<td>NA</td>";
	}
	$html .= "</tr>\n";

	return $html;
}


sub variant_row_tables
{
	my $geneList = shift;
	my $toPrint = shift;	# comma separated string e.g. "ssm,sv,cnv"
	my $genes = shift;
	my $chr = shift;
	my $germOrSom = shift;

	my $firstVar = 1;
    my %ssmTypes = (
        "nonsynonymous" => 1,
        "splicing" => 1,
        "stopgain" => 1,
        "unknown" => 1,
        "splicing" => 1,
        "nonframeshift" => 1,
        "frameshift" => 1,
        "stoploss" => 1,
    );

    my %svTypes = (
        "deletion breakpoint" => 1,
        "duplication breakpoint" => 1,
        "inversion breakpoint" => 1,
        "translocation breakpoint" => 1,
	);

    my %cnvTypes = (
        "homozygous deletion" => 1,
        "strong amplification" => 1,
	);

	my $posChr;

	my $classTest = "";  # should match both somatic and germline

	if ($germOrSom eq "som")
	{
		$classTest = "somatic";
	}
	elsif ($germOrSom eq "germ")
	{
		$classTest = "germline";
	}

	my $html = "";

	if ($toPrint =~ /ssm/)
	{
		for my $g (@{$geneList})
		{
			if (exists $genes->{$g}{variants})
			{
				for my $v (sort keys %{ $genes->{$g}{variants} })
				{
					if ($genes->{$g}{variants}{$v}{mutation_class} =~ /$classTest/)
					{
						if (exists  $ssmTypes{ $genes->{$g}{variants}{$v}{mutation_type} })
						{
							$posChr = $genes->{$g}{gene_position};
							$posChr =~ s/\:.*$//;
							if (($chr eq "allchr") or ($chr eq $posChr))
							{
								if ($firstVar == 1)
								{
									$html .= "<hr>\n";
									$html .= "<table class=\"regular\"><th>Gene</th><th>Variant Type</th><th>Position</th><th>Ref>Alt</th><th>Nuc Change</th><th>AA Change</th><th>IDs</th></tr>\n";
									$firstVar = 0;
								}
	
								$html .= ssm_variant_row($g,$v,$genes);
							}
						}
					}
				}
			}
	
		}
		if ($firstVar == 0)
		{
			$html .= "</table>\n";
		}
	}

	$firstVar = 1;

	if ($toPrint =~ /sv/)
	{
		for my $g (@{$geneList})
		{
			if (exists $genes->{$g}{variants})
			{
				for my $v (sort keys %{ $genes->{$g}{variants} })
				{
					if ($genes->{$g}{variants}{$v}{mutation_class} =~ /$classTest/)
					{
						if (exists  $svTypes{ $genes->{$g}{variants}{$v}{mutation_type} })
						{
							$posChr = $genes->{$g}{gene_position};
							$posChr =~ s/\:.*$//;
							if (($chr eq "allchr") or ($chr eq $posChr))
							{
								if ($firstVar == 1)
								{
									$html .= "<hr>\n";
									$html .= "<table class=\"regular\"><th>Gene</th><th>Variant Type</th><th>Position</th><th>Potential Fusion Partner</th></tr>\n";

									$firstVar = 0;
								}

								$html .= sv_variant_row($g,$v,$genes);

							}
						}
					}
				}
			}
		}
		if ($firstVar == 0)
		{
			$html .= "</table>\n";
		}
	}

	$firstVar = 1;

	if ($toPrint =~ /cnv/)
	{
		for my $g (@{$geneList})
		{
			if (exists $genes->{$g}{variants})
			{
				for my $v (sort keys %{ $genes->{$g}{variants} })
				{
					if ($genes->{$g}{variants}{$v}{mutation_class} =~ /$classTest/)
					{
						if (exists  $cnvTypes{ $genes->{$g}{variants}{$v}{mutation_type} })
						{
							$posChr = $genes->{$g}{gene_position};
							$posChr =~ s/\:.*$//;
							if (($chr eq "allchr") or ($chr eq $posChr))
							{
								if ($firstVar == 1)
								{
									$html .= "<hr>\n";
									$html .= "<table class=\"regular\"><th>Gene</th><th>Variant Type</th><th>Position</th><th>Copy Number</th><th>AB Counts</th><th>Expression Fold Change</th></tr>\n";

									$firstVar = 0;
								}

								$html .= cnv_variant_row($g,$v,$genes);
							}

						}
					}
				}
			}
		}
		if ($firstVar == 0)
		{
			$html .= "</table>\n";
		}
	}

	return $html;

}


sub short_gene
{
	my $g = shift;
	my $genes = shift;
	my $spec = shift;	# "no_germ" or "no_som" to exclude germline or somatic variants

	unless (defined $spec)
	{
		$spec = "all";
	}

	my $col = "typical";
	my $report;
	my $germ;
	my $vars = "";

	if (exists $genes->{$g}{biallelic_lof})
	{
		$col = "biallelic";
	}
	elsif (exists $genes->{$g}{gof})
	{
		$col = "gof";
	}

	if (exists $genes->{$g}{variants})
	{
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			$report = 0;
			$germ = 0;
			if ($genes->{$g}{variants}{$v}{mutation_class} =~ /^germline/)
			{
				unless ($spec eq "no_germ")
				{
					# only report germline if they're confirmed in the varkb
					if (exists $genes->{$g}{variants}{$v}{varkb})
					{
						for my $kb (sort keys %{ $genes->{$g}{variants}{$v}{varkb} })
						{
							if ($genes->{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /function/)
							{
								$germ = 1;
								$report = 1;
							}
						}
					}
				}
			}
			elsif ($genes->{$g}{variants}{$v}{mutation_class} =~ /^somatic/)
			{
				unless ($spec eq "no_som")
				{
					$report = 1;
				}
			}

			if ($report == 1)
			{
				if ($germ == 1)
				{
					$vars .= ",g" . get_mut_letter($genes->{$g}{variants}{$v}{mutation_type});
				}
				else
				{
					$vars .= "," . get_mut_letter($genes->{$g}{variants}{$v}{mutation_type});
				}
			}
		}

	}

	$vars =~ s/^,//;

	if ($vars eq "")
	{
		return "<a href=\"#$g\" style=\"text-decoration: none\"><mark class=\"$col wt\"><i>$g</i><sub>$vars</sub></mark></a>";
	}
	else
	{
		return "<a href=\"#$g\" style=\"text-decoration: none\"><mark class=\"$col\"><i>$g</i><sub>$vars</sub></mark></a>";
	}
}





# need to rewite this so works if spec or varType aren't passed
sub slide_gene
{
	my $g = shift;
	my $genes = shift;
	my $spec = shift;	# "no_germ" or "no_som" to exclude germline or somatic variants, or "all"
	my $ploidy = shift;
	my $varType = shift; # slide for slide format, page for page formatted mutation names

	unless (defined $spec)
	{
		$spec = "all";
	}

	unless (defined $varType)
	{
		$varType = "slide";	# default to slide
	}

	my $col = "typical";
	my $report;
	my $germ;
	my %vars;
	my $varText;

	if (exists $genes->{$g}{biallelic_lof})
	{
		$col = "biallelic";
	}
	elsif (exists $genes->{$g}{gof})
	{
		$col = "gof";
	}

	if (exists $genes->{$g}{variants})
	{
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			$report = 0;
			$germ = 0;
			if ($genes->{$g}{variants}{$v}{mutation_class} =~ /^germline/)
			{
				unless ($spec eq "no_germ")
				{
					# only report germline if they're confirmed in the varkb
					if (exists $genes->{$g}{variants}{$v}{varkb})
					{
						for my $kb (sort keys %{ $genes->{$g}{variants}{$v}{varkb} })
						{
							if ($genes->{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /function/)
							{
								$germ = 1;
								$report = 1;
							}
						}
					}
				}
			}
			elsif ($genes->{$g}{variants}{$v}{mutation_class} =~ /^somatic/)
			{
				unless ($spec eq "no_som")
				{
					$report = 1;
				}
			}

			if ($report == 1)
			{
				if ($varType eq "page")
				{
					$varText = get_mut_page($g, $v, $genes);
					if (($germ == 1) and ($spec ne "no_som"))
					{
						$varText = "germline $varText";
					}
				}
				else
				{
					$varText = get_mut_slide($g, $v, $genes);
					if (($germ == 1) and ($spec ne "no_som"))
					{
						$varText = "$varText(g)";
					}
				}
				$vars{$varText}++;
			}
		}

	}

	my $cn = "NA";
	my $cn_col = "cn_typical";
	my $low_cn = 99999;

	if ((defined $genes->{$g}{copy_number}) and ($genes->{$g}{copy_number} ne "NA"))
	{
		for my $cns (split(/\|/, $genes->{$g}{copy_number}))
		{
			if (($cns ne "NA") and ($cns < $low_cn))
			{
				$low_cn = $cns;
			}
		}

		$cn = sprintf("%.0f", $low_cn);
		$cn =~ s/^-0$/0/;
		$cn_col = get_cn_mark_class($low_cn, $ploidy);
	}


	my $rna = "NA";
	my $rna_col = "rna0";

	if ((defined $genes->{$g}{rna}{foldchange}) and ($genes->{$g}{rna}{foldchange} ne "NA"))
	{
		$rna = sprintf("%.1f", $genes->{$g}{rna}{foldchange});
		$rna_col = get_rna_mark_class($genes->{$g}{rna}{foldchange});
	}


	my $loh = 0;
	for my $ab (split(/\|/, $genes->{$g}{ab_counts}))
	{
		if ($ab =~ /^0/)
		{
			$loh = 1;
		}
	}

	my $variants = "";
	for my $v (sort keys %vars)
	{
		if (($v =~ /^c\./) or ($v =~ /^p\./))
		{
			$variants .= $v . ", ";
			$vars{$v} = 0;
		}
	}
	for my $v (sort keys %vars)
	{
		if ($v =~ /SV/)
		{
			$variants .= $v . ", ";
			$vars{$v} = 0;
		}
	}
	for my $v (sort keys %vars)
	{
		if ($vars{$v} > 0)
		{
			$variants .= $v . ", ";
			$vars{$v} = 0;
		}
	}

	$variants =~ s/, $//;

	if (($loh > 0) and !($variants =~ /homDel/) and !($variants =~ /homozygous deletion/))
	{
		$variants .= ", LOH";
	}

	$variants =~ s/^, //;


	if ($variants eq "")
	{
		$variants = "WT";
	}


	return "<tr><td class=\"$col\"><a href=\"#$g\" style=\"text-decoration: none\"><mark class=\"$col\"><b><i>$g</i></b></mark></a></td><td>$variants</td><td class=\"$cn_col\" align=\"center\">$cn</td><td class=\"$rna_col\" align=\"center\">$rna</td></tr>\n";
}


sub variant_to_report
{
	my $g = shift;
	my $genes = shift;
	my $meta = shift;
	my $class = shift;	# germ or som

	my $report = 0;

	# check if forced include or forced exclude
	
	for my $include_gene (split(/\|/, $meta->{include_genes}))
	{
		if ($include_gene eq $g)
		{
			$report = 1;
			return $report;
		}
	}

	for my $exclude_gene (split(/\|/, $meta->{exclude_genes}))
	{
		if ($exclude_gene eq $g)
		{
			$report = 0;
			return $report;
		}
	}

	# check if variant exists
	if (exists $genes->{$g}{variants})
	{
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			if ($genes->{$g}{variants}{$v}{mutation_class} =~ /^germline/)
			{
				unless ($class eq "som")
				{
					# only report germline if they're confirmed in the varkb
					if (exists $genes->{$g}{variants}{$v}{varkb})
					{
						for my $kb (sort keys %{ $genes->{$g}{variants}{$v}{varkb} })
						{
							if ($genes->{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /function/)
							{
								$report = 1;
							}
						}

					}
				}
			}
			elsif ($genes->{$g}{variants}{$v}{mutation_class} =~ /^somatic/)
			{
				unless ($class eq "germ")
				{
					$report = 1;
				}
			}

		}
	}

	return $report;

}



sub get_mut_page
{
	my $g = shift;
	my $v = shift;
	my $genes = shift;

	my $mut = "";

	my $type;

	my ($isoform,$change);
	

	if (($genes->{$g}{variants}{$v}{mutation_class} =~ /snp/) or ($genes->{$g}{variants}{$v}{mutation_class} =~ /indel/) or ($genes->{$g}{variants}{$v}{mutation_class} =~ /snv/))
	{
		$type = "NA";
		if ($genes->{$g}{variants}{$v}{mutation_class} =~ /snp/)
		{
			$type = "SNP";
		}
		elsif ($genes->{$g}{variants}{$v}{mutation_class} =~ /indel/)
		{
			$type = "indel";
		}
		elsif ($genes->{$g}{variants}{$v}{mutation_class} =~ /snv/)
		{
			$type = "SNV";
		}

		if ($genes->{$g}{variants}{$v}{mutation_type} eq "splicing")
		{
			($isoform,$change) = get_best_isoform($g,$v,$genes->{$g}{variants}{$v}{nuc_context},$genes);
			$change = "splice altering $type ($change)";
		}
		else
		{
			($isoform,$change) = get_best_isoform($g,$v,$genes->{$g}{variants}{$v}{aa_context},$genes);

			if ($genes->{$g}{variants}{$v}{mutation_type} eq "stopgain")
			{
				$change = "stopgain $type ($change)";
			}
			elsif ($genes->{$g}{variants}{$v}{mutation_type} eq "stoploss")
			{
				$change = "stoploss $type ($change)";
			}
			elsif ($genes->{$g}{variants}{$v}{mutation_type} eq "frameshift")
			{
				$change = "frameshift $type ($change)";
			}
			elsif ($genes->{$g}{variants}{$v}{mutation_type} eq "nonframeshift")
			{
				$change = "non-frameshift $type ($change)";
			}
			elsif ($genes->{$g}{variants}{$v}{mutation_type} eq "nonsynonymous")
			{
				$change = "nonsynonymous $type ($change)";
			}

		}

		$mut = $change;
	}
	elsif ($genes->{$g}{variants}{$v}{mutation_class} eq "somatic sv")
	{
		if ($genes->{$g}{variants}{$v}{mutation_type} eq "translocation breakpoint")
		{
			$mut = "translocation SV";
		}
		elsif ($genes->{$g}{variants}{$v}{mutation_type} eq "deletion breakpoint")
		{
			$mut = "deletion SV";
		}
		elsif ($genes->{$g}{variants}{$v}{mutation_type} eq "duplication breakpoint")
		{
			$mut = "duplication SV";
		}
		elsif ($genes->{$g}{variants}{$v}{mutation_type} eq "inversion breakpoint")
		{
			$mut = "inversion SV";
		}
	}
	elsif ($genes->{$g}{variants}{$v}{mutation_class} eq "somatic cnv")
	{
		if ($genes->{$g}{variants}{$v}{mutation_type} eq "homozygous deletion")
		{
			$mut = "homozygous deletion";
		}
		elsif ($genes->{$g}{variants}{$v}{mutation_type} eq "strong amplification")
		{
			$mut = "strong amplification"
		}
	}

	return $mut;
}




sub get_mut_slide
{
	my $g = shift;
	my $v = shift;
	my $genes = shift;

	my $mut = "";

	my ($isoform,$change);
	

	if (($genes->{$g}{variants}{$v}{mutation_class} =~ /snp/) or ($genes->{$g}{variants}{$v}{mutation_class} =~ /indel/) or ($genes->{$g}{variants}{$v}{mutation_class} =~ /snv/))
	{
		if ($genes->{$g}{variants}{$v}{mutation_type} eq "splicing")
		{
			($isoform,$change) = get_best_isoform($g,$v,$genes->{$g}{variants}{$v}{nuc_context},$genes);
		}
		else
		{
			($isoform,$change) = get_best_isoform($g,$v,$genes->{$g}{variants}{$v}{aa_context},$genes);
		}

		$mut = $change;
	}
	elsif ($genes->{$g}{variants}{$v}{mutation_class} eq "somatic sv")
	{
		if ($genes->{$g}{variants}{$v}{mutation_type} eq "translocation breakpoint")
		{
			$mut = "traSV";
		}
		elsif ($genes->{$g}{variants}{$v}{mutation_type} eq "deletion breakpoint")
		{
			$mut = "delSV";
		}
		elsif ($genes->{$g}{variants}{$v}{mutation_type} eq "duplication breakpoint")
		{
			$mut = "dupSV";
		}
		elsif ($genes->{$g}{variants}{$v}{mutation_type} eq "inversion breakpoint")
		{
			$mut = "invSV";
		}
	}
	elsif ($genes->{$g}{variants}{$v}{mutation_class} eq "somatic cnv")
	{
		if ($genes->{$g}{variants}{$v}{mutation_type} eq "homozygous deletion")
		{
			$mut = "homDel";
		}
		elsif ($genes->{$g}{variants}{$v}{mutation_type} eq "strong amplification")
		{
			$mut = "cnGain"
		}
	}

	return $mut;
}





sub full_gene
{
	my $g = shift;
	my $genes = shift;

	my $pathToGeneReport = "../../../genes/g/";

	my $html = "";

	# gene info (name, full, synonyms
	
	$html .= "<table id=\"$g\" class=\"full_gene_header\"><tr><td align=\"left\"><i>$g</i></td><td align=\"center\">$genes->{$g}{full_name}</td><td align=\"right\">$genes->{$g}{gene_chr}</td></tr></table>\n";

	my @table;
	if (exists $genes->{$g}{biallelic_lof})
	{
		push(@table, "<b>Status:</b> <mark class=\"biallelic\">Inactivated</mark>");
	}
	elsif (exists $genes->{$g}{gof})
	{
		push(@table, "<b>Status:</b> <mark class=\"gof\">Activated</mark>");
	}
	else
	{
		push(@table, "<b>Status:</b> Intact");
	}
	push(@table, "<b>Position:</b> " . $genes->{$g}{gene_position});

	my $kbTags = "";
	if (exists $genes->{$g}{genekb})
	{
		for my $kb (sort keys %{ $genes->{$g}{genekb} })
		{
			$kbTags .= "$kb, ";
		}
	}
	$kbTags =~ s/, $//;

	push(@table, "<b>KB Tags:</b> $kbTags");

	$html .= print_table(3,0,"regular",\@table);
	@table = ();

	push(@table, "<b>Copy Number:</b> " . $genes->{$g}{copy_number});

	my $ab = $genes->{$g}{ab_counts};
	my $loh = "typical";

	my $rnaClass;
	my $cnClass;

	$ab =~ s/\./-/g;

	if ($ab =~ /^0/)
	{
		$loh = "loh";
	}
	push(@table, "<b>A-B Count:</b> <mark class=\"$loh\">$ab</mark>");

	if (exists $genes->{$g}{exp_rank})
	{
		push(@table, "<b>Expression Rank:</b> " . $genes->{$g}{exp_rank});
	}
	else
	{
		push(@table, "<b>Expression Rank:</b> NA");
	}

	if ((exists $genes->{$g}{rna}) and (defined $genes->{$g}{rna}{foldchange}))
	{
		$rnaClass = get_rna_mark_class($genes->{$g}{rna}{foldchange});
		push(@table, "<b>Expression Fold Change:</b> <mark class=\"$rnaClass\">" . sprintf("%.3f", $genes->{$g}{rna}{foldchange}) . "</mark>");
	}
	else
	{
		push(@table, "<b>Expression Fold Change:</b> NA");
	}

	$html .= print_table(4,0,"regular",\@table);
	@table = ();

	$html .= "<p><a href=\"$pathToGeneReport/$g.html\">Cohort report for $g</a>\n";

	my $report;
	my $first;
	my $ids;
	my ($isoform, $change);

	my %varText = ();
	my $varTextCounter = 1;
	my $varTextInsert;

	my $noVariants = 1;
	
	$html .= "<hr>\n";
	# variant list (if exists)
	if (exists $genes->{$g}{variants})
	{
		$first = 1;
		@table = ();
	
		# do all snps and indels first
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			if (($genes->{$g}{variants}{$v}{mutation_class} =~ /snp/) or ($genes->{$g}{variants}{$v}{mutation_class} =~ /indel/) or ($genes->{$g}{variants}{$v}{mutation_class} =~ /snv/))
			{
				$report = 0;
				if ($genes->{$g}{variants}{$v}{mutation_class} =~ /^germline/)
				{
					# only report germline if they're confirmed in the varkb
					if (exists $genes->{$g}{variants}{$v}{varkb})
					{
						for my $kb (sort keys %{ $genes->{$g}{variants}{$v}{varkb} })
						{
							if ($genes->{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /function/)
							{
								$report = 1;

								if ((defined $genes->{$g}{variants}{$v}{varkb}{$kb}{report_text}) and ($genes->{$g}{variants}{$v}{varkb}{$kb}{report_text} ne ""))
								{
									$varText{$varTextCounter} .= $genes->{$g}{variants}{$v}{varkb}{$kb}{report_text} . "  ";
								}
							}
						}
					}
				}
				else
				{
					$report = 1;

					if (exists $genes->{$g}{variants}{$v}{varkb})
					{
						for my $kb (sort keys %{ $genes->{$g}{variants}{$v}{varkb} })
						{
							if ((defined $genes->{$g}{variants}{$v}{varkb}{$kb}{report_text}) and ($genes->{$g}{variants}{$v}{varkb}{$kb}{report_text} ne ""))
							{
								$varText{$varTextCounter} .= $genes->{$g}{variants}{$v}{varkb}{$kb}{report_text} . "  ";
							}
						}

					}
				}

				if ($report == 1)
				{
					$noVariants = 0;
					if ($first == 1)
					{
						$first = 0;
						# print header
						push(@table, "Variant Type", "Position", "Ref>Alt", "Nuc Change", "AA Change", "T Freq (Depth)", "N Freq (Depth)", "IDs");
					}

					$varTextInsert = "";
					if (exists $varText{$varTextCounter})
					{
						$varTextInsert = "($varTextCounter)";
						$varTextCounter++;
					}


					push (@table, "$varTextInsert $genes->{$g}{variants}{$v}{mutation_class} $genes->{$g}{variants}{$v}{mutation_type}");
					push (@table, $genes->{$g}{variants}{$v}{position});
					push (@table, $genes->{$g}{variants}{$v}{base_change});
					($isoform,$change) = get_best_isoform($g,$v,$genes->{$g}{variants}{$v}{nuc_context},$genes);
					push (@table, $change);
					($isoform,$change) = get_best_isoform($g,$v,$genes->{$g}{variants}{$v}{aa_context},$genes);
					push (@table, $change);

					push (@table, "$genes->{$g}{variants}{$v}{tumour_freq}% ($genes->{$g}{variants}{$v}{tumour_depth})");
					push (@table, "$genes->{$g}{variants}{$v}{normal_freq}% ($genes->{$g}{variants}{$v}{normal_depth})");
					
					$ids = "$isoform, ";
					for my $i (qw/dbsnp cosmic/)
					{
						if ((defined $genes->{$g}{variants}{$v}{$i}) and ($genes->{$g}{variants}{$v}{$i} ne "") and ($genes->{$g}{variants}{$v}{$i} ne "."))
						{
							$ids .= "$genes->{$g}{variants}{$v}{$i}, "
						}
					}
					$ids =~ s/, $//;
					push (@table, $ids);



				}
			}
		}
		unless ($noVariants == 1)
		{
			$html .= print_table(8,1,"variant",\@table);
		}

		# sv breakpoints next
		$first = 1;
		@table = ();
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			if ($genes->{$g}{variants}{$v}{mutation_class} eq "somatic sv")
			{
				if ($first == 1)
				{
					$first = 0;
					unless ($noVariants == 1)
					{
						$html .= "<hr>\n";
					}
					push(@table, "Rearrangement Type", "Positions", "Gene(s) at Other Breakpoint");
				}

				$noVariants = 0;

				push(@table, "$genes->{$g}{variants}{$v}{mutation_class} $genes->{$g}{variants}{$v}{mutation_type}");
				push(@table, $genes->{$g}{variants}{$v}{position});
				push(@table, $genes->{$g}{variants}{$v}{fusion_genes});

			}
		}
		unless ($noVariants == 1)
		{
			$html .= print_table(3,1,"variant",\@table);
		}


		# cnv last
		$first = 1;
		@table = ();
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			if ($genes->{$g}{variants}{$v}{mutation_class} eq "somatic cnv")
			{
				if ($first == 1)
				{
					$first = 0;
					unless ($noVariants == 1)
					{
						$html .= "<hr>\n";
					}
					push(@table, "Copy Number Variation Type", "Position");
				}

				$noVariants = 0;

				push(@table, "$genes->{$g}{variants}{$v}{mutation_class} $genes->{$g}{variants}{$v}{mutation_type}");
				push(@table, $genes->{$g}{variants}{$v}{position});

			}
		}
		unless ($noVariants == 1)
		{
			$html .= print_table(2,1,"variant",\@table);
		}
	}


	if ($noVariants == 1)
	{
		$html .= "<p>No variants of interest detected.</p>\n";
	}


	my $text = "";

	# kb info
	if (exists $genes->{$g}{genekb})
	{
		for my $tag (sort keys %{ $genes->{$g}{genekb} })
		{
			for my $kb (sort keys %{ $genes->{$g}{genekb}{$tag} })
			{
				if (exists $genes->{$g}{genekb}{$tag}{$kb}{report_text})
				{
					$text .= $genes->{$g}{genekb}{$tag}{$kb}{report_text} . " ";
				}
			}
		}
	}
	unless ($text eq "")
	{
		$text .= "<br>\n";
	}

	if (exists $varText{1})
	{
		for my $i (sort keys %varText)
		{
			$text .= "<br>";
			$text .= "($i) $varText{$i}";
		}
	}

	unless ($text eq "")
	{
		$html .= print_text_box($text);
	}

	return $html;

}





sub full_gene_slide
{
	my $g = shift;
	my $genes = shift;

	my $pathToGeneReport = "../../../genes/g/";

	my $html = "";

	# gene info (name, full, synonyms
	
	$html .= "<table id=\"$g\" class=\"full_gene_header\"><tr><td align=\"left\"><i>$g</i></td><td align=\"center\">$genes->{$g}{full_name}</td><td align=\"right\">$genes->{$g}{gene_chr}</td></tr></table>\n";

	my @table;
	if (exists $genes->{$g}{biallelic_lof})
	{
		push(@table, "<b>Status:</b> <mark class=\"biallelic\">Inactivated</mark>");
	}
	elsif (exists $genes->{$g}{gof})
	{
		push(@table, "<b>Status:</b> <mark class=\"gof\">Activated</mark>");
	}
	else
	{
		push(@table, "<b>Status:</b> Intact");
	}
	push(@table, "<b>Position:</b> " . $genes->{$g}{gene_position});

#	my $kbTags = "";
#	if (exists $genes->{$g}{genekb})
#	{
#		for my $kb (sort keys %{ $genes->{$g}{genekb} })
#		{
#			$kbTags .= "$kb, ";
#		}
#	}
#	$kbTags =~ s/, $//;

#	push(@table, "<b>KB Tags:</b> $kbTags");

	$html .= print_table(2,0,"regular",\@table);
	@table = ();

	push(@table, "<b>Copy Number:</b> " . $genes->{$g}{copy_number});

	my $ab = $genes->{$g}{ab_counts};
	my $loh = "typical";

	my $rnaClass;
	my $cnClass;

	$ab =~ s/\./-/g;

	if ($ab =~ /^0/)
	{
		$loh = "loh";
	}
	push(@table, "<b>A-B Count:</b> <mark class=\"$loh\">$ab</mark>");

#	if (exists $genes->{$g}{exp_rank})
#	{
#		push(@table, "<b>Expression Rank:</b> " . $genes->{$g}{exp_rank});
#	}
#	else
#	{
#		push(@table, "<b>Expression Rank:</b> NA");
#	}

	if ((exists $genes->{$g}{rna}) and ($genes->{$g}{rna} ne ""))
	{
		$rnaClass = get_rna_mark_class($genes->{$g}{rna}{foldchange});
		push(@table, "<b>Expression Fold Change:</b> <mark class=\"$rnaClass\">" . sprintf("%.3f", $genes->{$g}{rna}{foldchange}) . "</mark>");
	}
	else
	{
		push(@table, "<b>Expression Fold Change:</b> NA");
	}

	$html .= print_table(3,0,"regular",\@table);
	@table = ();

	my $report;
	my $first;
	my $ids;
	my ($isoform, $change);

	my %varText = ();
	my $varTextCounter = 1;
	my $varTextInsert;

	my $noVariants = 1;
	
	$html .= "<hr>\n";
	# variant list (if exists)
	if (exists $genes->{$g}{variants})
	{
		$first = 1;
		@table = ();
	
		# do all snps and indels first
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			if (($genes->{$g}{variants}{$v}{mutation_class} =~ /snp/) or ($genes->{$g}{variants}{$v}{mutation_class} =~ /indel/) or ($genes->{$g}{variants}{$v}{mutation_class} =~ /snv/))
			{
				$report = 0;
				if ($genes->{$g}{variants}{$v}{mutation_class} =~ /^germline/)
				{
					# only report germline if they're confirmed in the varkb
					if (exists $genes->{$g}{variants}{$v}{varkb})
					{
						for my $kb (sort keys %{ $genes->{$g}{variants}{$v}{varkb} })
						{
							if ($genes->{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /function/)
							{
								$report = 1;

								if ((defined $genes->{$g}{variants}{$v}{varkb}{$kb}{report_text}) and ($genes->{$g}{variants}{$v}{varkb}{$kb}{report_text} ne ""))
								{
									$varText{$varTextCounter} .= $genes->{$g}{variants}{$v}{varkb}{$kb}{report_text} . "  ";
								}
							}
						}
					}
				}
				else
				{
					$report = 1;

					if (exists $genes->{$g}{variants}{$v}{varkb})
					{
						for my $kb (sort keys %{ $genes->{$g}{variants}{$v}{varkb} })
						{
							if ((defined $genes->{$g}{variants}{$v}{varkb}{$kb}{report_text}) and ($genes->{$g}{variants}{$v}{varkb}{$kb}{report_text} ne ""))
							{
								$varText{$varTextCounter} .= $genes->{$g}{variants}{$v}{varkb}{$kb}{report_text} . "  ";
							}
						}

					}
				}

				if ($report == 1)
				{
					$noVariants = 0;
					if ($first == 1)
					{
						$first = 0;
						# print header
						push(@table, "Variant Type", "Position", "Nuc Change", "AA Change", "T Freq (Depth)", "N Freq (Depth)");
					}

					$varTextInsert = "";
					if (exists $varText{$varTextCounter})
					{
						$varTextInsert = "($varTextCounter)";
						$varTextCounter++;
					}


					push (@table, "$varTextInsert $genes->{$g}{variants}{$v}{mutation_class} $genes->{$g}{variants}{$v}{mutation_type}");
					push (@table, $genes->{$g}{variants}{$v}{position});
#					push (@table, $genes->{$g}{variants}{$v}{base_change});
					($isoform,$change) = get_best_isoform($g,$v,$genes->{$g}{variants}{$v}{nuc_context},$genes);
					push (@table, $change);
					($isoform,$change) = get_best_isoform($g,$v,$genes->{$g}{variants}{$v}{aa_context},$genes);
					push (@table, $change);

					push (@table, "$genes->{$g}{variants}{$v}{tumour_freq}% ($genes->{$g}{variants}{$v}{tumour_depth})");
					push (@table, "$genes->{$g}{variants}{$v}{normal_freq}% ($genes->{$g}{variants}{$v}{normal_depth})");
					
					$ids = "$isoform, ";
					for my $i (qw/dbsnp cosmic/)
					{
						if ((defined $genes->{$g}{variants}{$v}{$i}) and ($genes->{$g}{variants}{$v}{$i} ne "") and ($genes->{$g}{variants}{$v}{$i} ne "."))
						{
							$ids .= "$genes->{$g}{variants}{$v}{$i}, "
						}
					}
					$ids =~ s/, $//;
#					push (@table, $ids);



				}
			}
		}
		unless ($noVariants == 1)
		{
			$html .= print_table(6,1,"variant",\@table);
		}

		# sv breakpoints next
		$first = 1;
		@table = ();
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			if ($genes->{$g}{variants}{$v}{mutation_class} eq "somatic sv")
			{
				if ($first == 1)
				{
					$first = 0;
					unless ($noVariants == 1)
					{
						$html .= "<hr>\n";
					}
					push(@table, "Rearrangement Type", "Positions", "Gene(s) at Other Breakpoint");
				}

				$noVariants = 0;

				push(@table, "$genes->{$g}{variants}{$v}{mutation_class} $genes->{$g}{variants}{$v}{mutation_type}");
				push(@table, $genes->{$g}{variants}{$v}{position});
				push(@table, $genes->{$g}{variants}{$v}{fusion_genes});

			}
		}
		unless ($noVariants == 1)
		{
			$html .= print_table(3,1,"variant",\@table);
		}


		# cnv last
		$first = 1;
		@table = ();
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			if ($genes->{$g}{variants}{$v}{mutation_class} eq "somatic cnv")
			{
				if ($first == 1)
				{
					$first = 0;
					unless ($noVariants == 1)
					{
						$html .= "<hr>\n";
					}
					push(@table, "Copy Number Variation Type", "Position");
				}

				$noVariants = 0;

				push(@table, "$genes->{$g}{variants}{$v}{mutation_class} $genes->{$g}{variants}{$v}{mutation_type}");
				push(@table, $genes->{$g}{variants}{$v}{position});

			}
		}
		unless ($noVariants == 1)
		{
			$html .= print_table(2,1,"variant",\@table);
		}
	}


	if ($noVariants == 1)
	{
		$html .= "<p>No variants of interest detected.</p>\n";
	}


	my $text = "";

	# kb info
	if (exists $genes->{$g}{genekb})
	{
		for my $tag (sort keys %{ $genes->{$g}{genekb} })
		{
			for my $kb (sort keys %{ $genes->{$g}{genekb}{$tag} })
			{
				if (exists $genes->{$g}{genekb}{$tag}{$kb}{report_text})
				{
					$text .= $genes->{$g}{genekb}{$tag}{$kb}{report_text} . " ";
				}
			}
		}
	}
	unless ($text eq "")
	{
		$text .= "<br>\n";
	}

	if (exists $varText{1})
	{
		for my $i (sort keys %varText)
		{
			$text .= "<br>";
			$text .= "($i) $varText{$i}";
		}
	}

	unless ($text eq "")
	{
		$html .= print_text_box($text);
	}

	return $html;

}






sub full_gene_page
{
	my $g = shift;
	my $genes = shift;

	my $pathToGeneReport = "../../../genes/g/";

	my $html = "";

	# gene info (name, full, synonyms
	
	$html .= "<table id=\"$g\" class=\"full_gene_header\"><tr><td align=\"left\" style=\"padding:0px\"><i>$g</i></td><td align=\"center\" style=\"padding:0px\">$genes->{$g}{full_name}</td><td align=\"right\" style=\"padding:0px\">$genes->{$g}{gene_chr}</td></tr></table>\n";

	my @table;
	if (exists $genes->{$g}{biallelic_lof})
	{
		push(@table, "<b>Status:</b> <mark class=\"biallelic\">Inactivated</mark>");
	}
	elsif (exists $genes->{$g}{gof})
	{
		push(@table, "<b>Status:</b> <mark class=\"gof\">Activated</mark>");
	}
	else
	{
		push(@table, "<b>Status:</b> Intact");
	}

#	my $kbTags = "";
#	if (exists $genes->{$g}{genekb})
#	{
#		for my $kb (sort keys %{ $genes->{$g}{genekb} })
#		{
#			$kbTags .= "$kb, ";
#		}
#	}
#	$kbTags =~ s/, $//;

#	push(@table, "<b>KB Tags:</b> $kbTags");


	push(@table, "<b>Copy Number:</b> " . $genes->{$g}{copy_number});

	my $ab = $genes->{$g}{ab_counts};
	my $loh = "typical";

	my $rnaClass;
	my $cnClass;

	$ab =~ s/\./-/g;

	if ($ab =~ /^0/)
	{
		$loh = "loh";
	}
	push(@table, "<b>A-B Count:</b> <mark class=\"$loh\">$ab</mark>");

#	if (exists $genes->{$g}{exp_rank})
#	{
#		push(@table, "<b>Expression Rank:</b> " . $genes->{$g}{exp_rank});
#	}
#	else
#	{
#		push(@table, "<b>Expression Rank:</b> NA");
#	}

	if ((exists $genes->{$g}{rna}) and ($genes->{$g}{rna} ne ""))
	{
		if (exists $genes->{$g}{rna}{foldchange}) {
			$rnaClass = get_rna_mark_class($genes->{$g}{rna}{foldchange});
		} else {
			$genes->{$g}{rna}{foldchange} = 0;
			$rnaClass = 0;
		}
		push(@table, "<b>Expression Fold Change:</b> <mark class=\"$rnaClass\">" . sprintf("%.3f", $genes->{$g}{rna}{foldchange}) . "</mark>");
	}
	else
	{
		push(@table, "<b>RNA TPM Fold Change:</b> NA");
	}

	$html .= print_table(4,0,"regular",\@table);
	@table = ();

	my $report;
	my $first;
	my $ids;
	my ($isoform, $change);

	my %varText = ();
	my $varTextCounter = 1;
	my $varTextInsert;

	my $noVariants = 1;
	
	$html .= "<hr>\n";
	# variant list (if exists)
	if (exists $genes->{$g}{variants})
	{
		$first = 1;
		@table = ();
	
		# do all snps and indels first
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			if (($genes->{$g}{variants}{$v}{mutation_class} =~ /snp/) or ($genes->{$g}{variants}{$v}{mutation_class} =~ /indel/) or ($genes->{$g}{variants}{$v}{mutation_class} =~ /snv/))
			{
				$report = 0;
				if ($genes->{$g}{variants}{$v}{mutation_class} =~ /^germline/)
				{
					# only report germline if they're confirmed in the varkb
					if (exists $genes->{$g}{variants}{$v}{varkb})
					{
						for my $kb (sort keys %{ $genes->{$g}{variants}{$v}{varkb} })
						{
							if ($genes->{$g}{variants}{$v}{varkb}{$kb}{effect} =~ /function/)
							{
								$report = 1;

								if ((defined $genes->{$g}{variants}{$v}{varkb}{$kb}{report_text}) and ($genes->{$g}{variants}{$v}{varkb}{$kb}{report_text} ne ""))
								{
									$varText{$varTextCounter} .= $genes->{$g}{variants}{$v}{varkb}{$kb}{report_text} . "  ";
								}
							}
						}
					}
				}
				else
				{
					$report = 1;

					if (exists $genes->{$g}{variants}{$v}{varkb})
					{
						for my $kb (sort keys %{ $genes->{$g}{variants}{$v}{varkb} })
						{
							if ((defined $genes->{$g}{variants}{$v}{varkb}{$kb}{report_text}) and ($genes->{$g}{variants}{$v}{varkb}{$kb}{report_text} ne ""))
							{
								$varText{$varTextCounter} .= $genes->{$g}{variants}{$v}{varkb}{$kb}{report_text} . "  ";
							}
						}

					}
				}

				if ($report == 1)
				{
					$noVariants = 0;
					if ($first == 1)
					{
						$first = 0;
						# print header
						push(@table, "Variant Type", "Position", "Nuc Change", "AA Change", "T Freq (Depth)", "N Freq (Depth)");
					}

					$varTextInsert = "";
					if (exists $varText{$varTextCounter})
					{
						$varTextInsert = "($varTextCounter)";
						$varTextCounter++;
					}


					push (@table, "$varTextInsert $genes->{$g}{variants}{$v}{mutation_class} $genes->{$g}{variants}{$v}{mutation_type}");
					push (@table, $genes->{$g}{variants}{$v}{position});
#					push (@table, $genes->{$g}{variants}{$v}{base_change});
					($isoform,$change) = get_best_isoform($g,$v,$genes->{$g}{variants}{$v}{nuc_context},$genes);
					push (@table, $change);
					($isoform,$change) = get_best_isoform($g,$v,$genes->{$g}{variants}{$v}{aa_context},$genes);
					push (@table, $change);

					push (@table, "$genes->{$g}{variants}{$v}{tumour_freq}% ($genes->{$g}{variants}{$v}{tumour_depth})");
					push (@table, "$genes->{$g}{variants}{$v}{normal_freq}% ($genes->{$g}{variants}{$v}{normal_depth})");
					
					$ids = "$isoform, ";
					for my $i (qw/dbsnp cosmic/)
					{
						if ((defined $genes->{$g}{variants}{$v}{$i}) and ($genes->{$g}{variants}{$v}{$i} ne "") and ($genes->{$g}{variants}{$v}{$i} ne "."))
						{
							$ids .= "$genes->{$g}{variants}{$v}{$i}, "
						}
					}
					$ids =~ s/, $//;
#					push (@table, $ids);



				}
			}
		}
		unless ($noVariants == 1)
		{
			$html .= print_table(6,1,"variant",\@table);
		}

		# sv breakpoints next
		$first = 1;
		@table = ();
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			if ($genes->{$g}{variants}{$v}{mutation_class} eq "somatic sv")
			{
				if ($first == 1)
				{
					$first = 0;
					unless ($noVariants == 1)
					{
						$html .= "<hr>\n";
					}
					push(@table, "Rearrangement Type", "Positions", "Gene(s) at Other Breakpoint");
				}

				$noVariants = 0;

				push(@table, "$genes->{$g}{variants}{$v}{mutation_class} $genes->{$g}{variants}{$v}{mutation_type}");
				push(@table, $genes->{$g}{variants}{$v}{position});
				push(@table, $genes->{$g}{variants}{$v}{fusion_genes});

			}
		}
		unless ($noVariants == 1)
		{
			$html .= print_table(3,1,"variant",\@table);
		}


		# cnv last
		$first = 1;
		@table = ();
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			if ($genes->{$g}{variants}{$v}{mutation_class} eq "somatic cnv")
			{
				if ($first == 1)
				{
					$first = 0;
					unless ($noVariants == 1)
					{
						$html .= "<hr>\n";
					}
					push(@table, "Copy Number Variation Type", "Position");
				}

				$noVariants = 0;

				push(@table, "$genes->{$g}{variants}{$v}{mutation_class} $genes->{$g}{variants}{$v}{mutation_type}");
				push(@table, $genes->{$g}{variants}{$v}{position});

			}
		}
		unless ($noVariants == 1)
		{
			$html .= print_table(2,1,"variant",\@table);
		}
	}


	if ($noVariants == 1)
	{
		$html .= "<p>No variants of interest detected.</p>\n";
	}


	my $text = "";

	# kb info
	if (exists $genes->{$g}{genekb})
	{
		for my $tag (sort keys %{ $genes->{$g}{genekb} })
		{
			for my $kb (sort keys %{ $genes->{$g}{genekb}{$tag} })
			{
				if (exists $genes->{$g}{genekb}{$tag}{$kb}{report_text})
				{
					$text .= $genes->{$g}{genekb}{$tag}{$kb}{report_text} . " ";
				}
			}
		}
	}
	unless ($text eq "")
	{
		$text .= "<br>\n";
	}

	if (exists $varText{1})
	{
		for my $i (sort keys %varText)
		{
			$text .= "<br>";
			$text .= "($i) $varText{$i}";
		}
	}

	unless ($text eq "")
	{
		$html .= print_text_box($text);
	}

	return $html;

}






sub get_preferred_isoform
{
	my $g = shift;
	my $genes = shift;
	my $kbtag = shift;	# "any" to check all tags


	my $prefIsoform = "";

	my ($nm1, $nm2);

	if (exists $genes->{$g}{genekb})
	{
		for my $tag (sort keys %{ $genes->{$g}{genekb} })
		{
			if (($tag eq $kbtag) or ($kbtag eq "any"))
			{
				for my $kb (sort keys %{ $genes->{$g}{genekb}{$tag} })
				{
					if ((defined $genes->{$g}{genekb}{$tag}{$kb}{preferred_isoform}) and ($genes->{$g}{genekb}{$tag}{$kb}{preferred_isoform} ne ""))
					{
						if (($prefIsoform eq "") or ($prefIsoform eq $genes->{$g}{genekb}{$tag}{$kb}{preferred_isoform}))
						{
							$prefIsoform = $genes->{$g}{genekb}{$tag}{$kb}{preferred_isoform};
						}
						else
						{
							warn "Preferred isoform conflict for $g in kb $kbtag ($prefIsoform or $genes->{$g}{genekb}{$tag}{$kb}{preferred_isoform}).  Using lower ID.\n";

							$nm1 = 999999999999;
							$nm2 = 999999999999;

							if ($prefIsoform =~ /NM_(.*)$/)
							{
								$nm1 = $1;
							}
							if ($genes->{$g}{genekb}{$tag}{$kb}{preferred_isoform} =~ /NM_(.*)$/)
							{
								$nm2 = $1;
							}

							if ($nm2 < $nm1)
							{
								$prefIsoform = $genes->{$g}{genekb}{$tag}{$kb}{preferred_isoform};
							}
						}
					}
				}

			}

		}
	}

	return $prefIsoform;
}


sub get_best_isoform
{
	my $g = shift;
	my $v = shift;
	my $list = shift;
	my $genes = shift;

	# "nuc_context" : "NM_007298:c.A1525G|NM_007300:c.A4900G|NM_007299:c.A1525G|NM_007297:c.A4696G|NM_007294:c.A4837G",
	# "aa_context" : "NM_007298:p.S509G|NM_007300:p.S1634G|NM_007299:p.S509G|NM_007297:p.S1566G|NM_007294:p.S1613G",
	
	my $prefIsoform = "";
	my $tag;

	# check gene kb for preferred isoform
	$prefIsoform = get_preferred_isoform($g,$genes,"any");

	my ($iso, $context);
	my $bestContext = "";
	my $bestIso = "";

	for my $f (sort split (/\|/, $list))
	{
		($iso, $context) = split(/:/, $f);
		if ($iso eq $prefIsoform)
		{
			$bestContext = $context;
			$bestIso = $iso;
		}
		elsif ($bestContext eq "")
		{
			$bestContext = $context;
			$bestIso = $iso;
		}
	}

	return ($iso, $bestContext);
}



sub pass_icon
{
	my $common = set_common_path();
	return "<img src=\"$common/pass.png\" style=\"height:20px\">";
}

sub fail_icon
{
	my $common = set_common_path();
	return "<img src=\"$common/fail.png\" style=\"height:20px\">";
}

sub pass_icon_size
{
	my $size = shift;
	my $common = set_common_path();
	return "<img src=\"$common/pass.png\" style=\"height:${size}px\">";
}

sub fail_icon_size
{
	my $size = shift;
	my $common = set_common_path();
	return "<img src=\"$common/fail.png\" style=\"height:${size}px\">";
}

sub pass_or_fail_icon
{
	my $val = shift;
	my $ord = shift;
	my $cut = shift;

	my $common = set_common_path();

	if ($ord eq "min")
	{
		if ($val >= $cut)
		{
			return "<img src=\"$common/pass.png\" style=\"height:20px\">";
		}
		else
		{
			return "<img src=\"$common/fail.png\" style=\"height:20px\">";
		}
	}
	else
	{
		if ($val <= $cut)
		{
			return "<img src=\"$common/pass.png\" style=\"height:20px\">";
		}
		else
		{
			return "<img src=\"$common/fail.png\" style=\"height:20px\">";
		}

	}
}

sub pass_or_fail_icon_size
{
	my $val = shift;
	my $ord = shift;
	my $cut = shift;
	my $size = shift;

	my $common = set_common_path();

	if ($ord eq "min")
	{
		if ($val >= $cut)
		{
			return "<img src=\"$common/pass.png\" style=\"height:${size}px\">";
		}
		else
		{
			return "<img src=\"$common/fail.png\" style=\"height:${size}px\">";
		}
	}
	else
	{
		if ($val <= $cut)
		{
			return "<img src=\"$common/pass.png\" style=\"height:${size}px\">";
		}
		else
		{
			return "<img src=\"$common/fail.png\" style=\"height:${size}px\">";
		}

	}
}

sub print_table
{
	my $cols = shift @_;	# number of columns
	my $header = shift @_;	# 1 if the first row is a header
	my $class = shift @_;
	my $content = shift @_;	# ref to array of values to put in the table

	my $count = 0;
	my $inHeader = 0;

	if ($header == 1)
	{
		$inHeader = 1;
	}

	my $html = "<table class=\"$class\">\n";


	for my $con (@$content)
	{
		if ($count == 0)
		{
			$html .= "<tr>\n";
		}
		$count++;


		if ($inHeader == 1)
		{
			$html .= "<th>";
		}
		else
		{
			$html .= "<td>";
		}

		$html .= $con;

		if 	($inHeader == 1)
		{
			$html .= "</th>";
		}
		else
		{
			$html .= "</td>";
		}


		if ($count >= $cols)
		{
			$html .= "</tr>\n";
			$count = 0;
			$inHeader = 0;
		}
	}
	unless ($count == 0)
	{
		$html .= "</tr>";
	}

	$html .= "</table>\n";


	return $html;
}

sub print_text_box
{
	my $text = shift;
	my $html = "";

	$html .= "<div class=\"comment_box\">\n";
	$html .= "<p style=\"margin-top:0px;margin-bottom:0px\">$text</p>\n";
	$html .= "</div>\n";

	return $html;

}

sub print_gene_box
{
}


sub print_mof_gene_box
{
	my $html = "";
	my $genes = shift;


	$html .= "<table width=\"100%\"><tr><td style=\"text-align:center\" class=\"mofclassic\" width=\"50%\">Classic Expression Fold Change</td><td style=\"text-align:center\" class=\"mofbasal-like\" width=\"50%\">Basal-like Expression Fold Change</td></tr></table>\n";

	my @mofBGenes = get_gene_list("mof_basal", $genes);
	my @mofCGenes = get_gene_list("mof_classic", $genes);

	my $exp;

	my @table = ();
	for my $g (sort @mofCGenes)
	{
		$exp = get_rna_mark_class($genes->{$g}{rna}{foldchange});
		push (@table, "<i><mark class=\"$exp\">$g</mark></i>");
	}

	$html .= "<table width=\"50%\" style=\"float:left\">\n";

	my $pos;
	for (my $i = 0; $i < 5; $i++)
	{
		$html .= "<tr>\n";
		for (my $j = 0; $j < 5; $j++)
		{
			$pos = ($i * 5) + $j;
			$html .= "<td style=\"text-align:center\" width=\"20%\">$table[$pos]</td>\n";
		}
		$html .= "</tr>\n";
	}
	$html .= "</table>\n";



	@table = ();
	for my $g (sort @mofBGenes)
	{
		$exp = get_rna_mark_class($genes->{$g}{rna}{foldchange});
		push (@table, "<i><mark class=\"$exp\">$g</mark></i>");
	}

	$html .= "<table width=\"49%\">\n";

	for (my $i = 0; $i < 5; $i++)
	{
		$html .= "<tr>\n";
		for (my $j = 0; $j < 5; $j++)
		{
			$pos = ($i * 5) + $j;
			$html .= "<td style=\"text-align:center\" width=\"20%\">$table[$pos]</td>\n";
		}
		$html .= "</tr>\n";
	}
	$html .= "</table>\n";

	return $html;
}



sub get_gene_list
{
	my $class = shift @_;
	my $genes = shift @_;

	my %list;

	for my $g (keys %{ $genes })
	{
		if (exists $genes->{$g}{genekb}{$class})
		{
			$list{$g}++;
		}
	}

	return sort keys %list;

}


sub commify
{
	my $text = reverse $_[0];
	$text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text;
}

1;

