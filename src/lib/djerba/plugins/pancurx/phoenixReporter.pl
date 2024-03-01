#!/usr/bin/perl

# use strict;
# use warnings;
# use Time::Local;

# my $summaryFile = shift;
# my $variantsFile = shift;
# my $outPath = shift;

# unless (defined $outPath)
# {
# 	$outPath = ".";
# }

# my (%data);
# $data{out_path} = $outPath;

# my $commonLocation = "/.mounts/labs/PCSI/reports/phoenix_common";


# my $l;
# my (@f, @header);
# $data{version} = "v0.2";

# # read summary file
# open (FILE, $summaryFile) or die "Couldn't open $summaryFile\n";

# while ($l = <FILE>)
# {
#     chomp $l;
#     if ($l =~ /^donor/)
#     {
#         @header = split(/,/, $l);
#     }
#     else
#     {
#         @f = split(/,/, $l);

#         for (my $i = 0; $i < scalar(@f); $i++)
#         {
#             $data{$header[$i]} = $f[$i];
#         }
#     }
# }
# close FILE;


# # read variants file
# my (%genes, %row);
# my $i;
# open (FILE, $variantsFile) or die "Couldn't open $variantsFile\n";

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
        $i = 0;
        for my $val (split(/,/, $l))
        {
            $row{$header[$i]} = $val;
            $i++;
        }
		for my $h (qw/ploidy copy_number ab_counts gene_position maf_mean maf_p/)
		{
			$genes{$row{gene}}{$h} = $row{$h};
		}

		for my $h (qw/mutation_class mutation_type position fusion_genes base_change tumour_freq tumour_depth normal_freq normal_depth nuc_context aa_context dbsnp cosmic cadd_phred rarity 1000G_all ExAC_all ESP6500siv2_all clinvar cosmic_census_flag cosmic_census_data/)
		{
			$genes{$row{gene}}{variants}{"$row{mutation_type},$row{position},$row{base_change}"}{$h} = $row{$h};
		}
    }
}
close FILE;


# $data{plot_dir} = "../plots/"; ##AMY htmls in results/reports/ but want to embed common and plots in results/ dir as before
# $data{common_dir} = "$data{plot_dir}/common"; ##AMY

# #`mkdir -p $outPath/plots`;  ##AMY
# `ln -s $commonLocation $outPath/$data{common_dir}`; ##AMY

# my $currentTime = localtime;

my @germlineGeneOrder = qw/POLE POLD1 EPCAM MLH1 MSH2 MSH6 PMS1 PMS2 BRCA1 BRCA2 PALB2 ATM APC MUTYH CDKN2A STK11 TP53 SMAD4 RAD51C CHEK2 BRIP1 RNF43/;
my @driverGeneOrder = qw/KRAS TP53 CDKN2A SMAD4 MAP2K4 ARID1A RNF43 TGFBR2 KDM6A/;

my %germlineGenes;
for my $g (@germlineGeneOrder)
{
	$germlineGenes{$g} = 1;
}

my %driverGenes;
for my $g (@driverGeneOrder)
{
	$driverGenes{$g} = 1;
}

my %nonsilent =  (
	"frameshift" => 1,
	"nonframeshift" => 1,
	"nonsynonymous" => 1,
	"stoploss" => 1,
	"stopgain" => 1,
	"splicing" => 1,
	"strong amplification" => 1,
	"homozygous deletion" => 1,
	"deletion breakpoint" => 1,
	"duplication breakpoint" => 1,
	"inversion breakpoint" => 1,
	"translocation breakpoint" => 1,
	# leaving altered promoters and potential fusions for now
);


# my %limsTissue = (
# 	"Ad" => "Adipose",
# 	"Ap" => "Appendix",
# 	"Ag" => "Adrenal Gland",
# 	"As" => "Ascites",
# 	"Bm" => "Bone Marrow",
# 	"Bn" => "Brain",
# 	"Br" => "Breast",
# 	"Bu" => "Buccal cells",
# 	"Cb" => "Cored Blood",
# 	"Cn" => "Central Nervous System",
# 	"Du" => "Duodenum",
# 	"Ep" => "Esophagus",
# 	"Es" => "Esophagus",
# 	"Fs" => "Foreskin",
# 	"Gb" => "Gallbladder",
# 	"Hr" => "Heart",
# 	"Ki" => "Kidney",
# 	"Le" => "Leukocyte",
# 	"Li" => "Large Intestine",
# 	"Ln" => "Lymph Node",
# 	"Lu" => "Lung",
# 	"Lv" => "Liver",
# 	"Lx" => "Larynx",
# 	"Ly" => "Lymphocyte",
# 	"Md" => "Mediastinum",
# 	"Me" => "Mesenchyme",
# 	"Nk" => "Neck",
# 	"Oc" => "Oral Cavity",
# 	"Om" => "Omentum",
# 	"Ov" => "Ovary",
# 	"Pa" => "Pancreas",
# 	"Pb" => "peripheral blood",
# 	"Pm" => "Peritoneum",
# 	"Pr" => "Prostate",
# 	"Sa" => "Saliva",
# 	"Sg" => "Salivary Gland",
# 	"Si" => "Small Intestine",
# 	"Sk" => "Skin",
# 	"Sm" => "Skeletal Muscle",
# 	"Sp" => "Spleen",
# 	"St" => "Stomach",
# 	"Ta" => "Tail (Referring to Models)",
# 	"Tr" => "Trachea",
# 	"Mu" => "Muscle",
# 	"Wm" => "Worm (Referring to Models)",
# 	"Um" => "Umbilical Cord",
# 	"nn" => "Unknown"
# );





# ### Print website
# my $html;
# open ($html, ">$outPath/$data{tumour}_summary.html") or die "Couldn't open >$outPath/$data{tumour}_summary.html\n";


# printHtmlHead($html, \%data);

# printHeader($html, \%data, $currentTime);

# printGenomicsSummary($html, \%data, \%limsTissue);

# printSeqMetrics($html, \%data);

# printClassifications($html, \%data, \%genes);


# printDsbrScore($html, \%data, \%genes);
# printMmrScore($html, \%data, \%genes);
# print $html "<p style=\"page-break-after:always;\"></p>\n";

# printGermlineVariants($html, \%data, \%genes, \@germlineGeneOrder, \%nonsilent);
print $html "<h2>Somatic Mutations</h2>\n";
printWholeGenome($html, \%data);
print $html "<p style=\"page-break-after:always;\"></p>\n";

printMutationLoad($html, \%data);

printDriverGenes($html, \%genes, \@driverGeneOrder, \%nonsilent);

print $html "<p style=\"page-break-after:always;\"></p>\n";


printPloidyAndChromo($html, \%data);

printFullSomaticList($html, \%data, \%genes, \%nonsilent);





#print $html "<p style=\"page-break-after:always;\"></p>\n";
#printLegend($html);



print $html "<p style=\"page-break-after:always;\"></p>\n";

print $html "<h1>Appendix</h1>\n";

#printLaneQC($html, \%data);

printWholeChroms($html, \%data);

print $html "</div>\n";

printFooter($html, \%data);
print $html "</body>\n</html>\n";








# sub printHtmlHead
# {
# 	my $html = shift;
# 	my $data = shift;

# 	print $html "<!DOCTYPE html>
# <html>
# <head>
# <title>$data->{tumour} Genome Summary</title>

# <link rel=\"stylesheet\" type=\"text/css\" href=\"$data->{common_dir}/report_style.css\">

# <script src=\"$data->{common_dir}/sorttable.js\"></script>
# <style type=\"text/css\">
# th, td {
#   padding: 3px !important;
# }
# table
# {
# border-collapse:collapse;
# }
# /* Sortable tables */
# table.sortable thead {
#         background-color:#FFFFFF;
#         color:#000000;
#         font-weight: bold;
#         cursor: default;
# }
# </style>

# </head>\n";

# }


# sub printHeader
# {
# 	my $html = shift;
# 	my $data = shift;
# 	my $currentTime = shift;

# 	print $html "<body>\n";

# 	print $html "<div id=\"header\">\n";
# 	print $html "<table border=\"0\" style=\"width:100%\"><tr style=\"vertical-align:top;\">
# 	<td style=\"width:33.33%\"><h4>Genome Summary $data->{version}</h4></td>
# 	<td style=\"width:33.33%\"><h4 style=\"text-align:center;\">$data{tumour}</h4></td>
# 	<td style=\"width:33.33%\"><h4 style=\"text-align:right;\">$currentTime</h4></td></tr></table>\n";

# 	print $html "</div>\n";
	
# 	print $html "<div id=\"main\">\n";


# }


# sub printGenomicsSummary
# {
# 	my $html = shift;
# 	my $data = shift;
# 	my $limsTissue = shift;

# 	my $tumTissue = "NA";
# 	if ($data{tumour} =~ /^[A-Z0-9]+_...._(..)_/)
# 	{
# 		if (exists $limsTissue->{$1})
# 		{
# 			$tumTissue = $limsTissue->{$1};
# 		}
# 		else
# 		{
# 			$tumTissue = $1;
# 		}
# 	}

# 	my $norTissue = "NA";
# 	if ($data{normal} =~ /^[A-Z0-9]+_...._(..)_/)
# 	{
# 		if (exists $limsTissue->{$1})
# 		{
# 			$norTissue = $limsTissue->{$1};
# 		}
# 		else
# 		{
# 			$norTissue = $1;
# 		}
# 	}


# 	print $html "<h1>PanCuRx Genomics Summary</h1>\n";

# 	print $html "<table style=\"width:100%\">\n";
# 	print $html "<tr>\n";
# 	print $html "<td><b>Donor ID:</b> $data->{donor}</td>";
# 	print $html "<td><b>Tumour Sample:</b> $data->{tumour}</td>";
# 	print $html "<td><b>Normal Sample:</b> $data->{normal}</td>";
# 	print $html "</tr>\n<tr>\n";
# 	print $html "<td><b>External IDs:</b> $data->{external_id}</td>";
# 	print $html "<td><b>Tumour Tissue:</b> $tumTissue</td>";
# 	print $html "<td><b>Normal Tissue:</b> $norTissue</td>";
# 	print $html "</tr>\n";
# 	print $html "</table>\n";

# 	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}-onco_slice-1000_no_name.png\" style=\"width:100%\">\n";

# 	print $html "<p><b>This report is for research purposes only.</b></p>\n";

# }

# sub printSeqMetrics
# {
# 	my $html = shift;
# 	my $data = shift;

# 	print $html "<h3>DNA Sequencing Quality Metrics</h3>\n";
	
# 	print $html "<table style=\"width:100%\">\n";
# 	print $html "<tr>\n";
# 	print $html "<td><b>Tumour Coverage:</b> " . sprintf("%.1f",$data{tumour_coverage}) . "x (>30x)</td>";
# #	print $html "<td><b>Error Rate:</b> " . sprintf("%.2f",$data{tumour_error_rate} * 100) . "% (<1%)</td>";
# #	print $html "<td><b>Soft Clip Rate:</b> " . sprintf("%.2f",$data{tumour_soft_clip} * 100) . "% (<1%)</td>";
# #	print $html "</tr>\n<tr>\n";
# 	print $html "<td><b>Normal Coverage:</b> " . sprintf("%.1f",$data{normal_coverage}) . "x (>30x)</td>";
# #	print $html "<td><b>Error Rate:</b> " . sprintf("%.2f",$data{normal_error_rate} * 100) . "% (<1%)</td>";
# #	print $html "<td><b>Soft Clip Rate:</b> " . sprintf("%.2f",$data{normal_soft_clip} * 100) . "% (<1%)</td>";
# #	print $html "</tr>\n";
# #	print $html "<tr>\n";
# 	print $html "<td><b>Cellularity:</b> " . ($data{cellularity} * 100) . "%</td>";
# #	print $html "<td><b>Genotype Concordance:</b> " . sprintf("%.01f", $data{genotype_concordance} * 100) . "%</td>";
# #	print $html "<td></td>";
# 	print $html "</tr>\n";
# 	print $html "</table>\n";
# }


# sub printGermlineVariants
# {
# 	# my $html = shift;
# 	# my $data = shift;
# 	my $genes = shift;
# 	# my $germlineGeneOrder = shift;
# 	# my $nonsilent = shift;

# 	# my $germVars = $data->{germline_snv_count} + $data->{germline_indel_count};
# 	# my $germNonsil = $data->{germline_missense_count} + $data->{germline_nonsense_count};
# 	my $germNonsilGenes = 0;
# 	my $germNonsilGenesRare = 0;
# 	my $germPathogenic = 0;

# 	for my $g (@{ $germlineGeneOrder })
# 	{
# 		for my $v (sort keys %{ $genes->{$g}{variants} })
# 		{
# 			if ($genes->{$g}{variants}{$v}{mutation_class} =~ /germline/)
# 			{
# 				if (exists $nonsilent->{ $genes->{$g}{variants}{$v}{mutation_type} })
# 				{
# 					$germNonsilGenes++;
# 					if ($genes->{$g}{variants}{$v}{rarity} ne "common")
# 					{
# 						$germNonsilGenesRare++;

# 						if (($genes->{$g}{variants}{$v}{mutation_type} eq "frameshift") or ($genes->{$g}{variants}{$v}{mutation_type} eq "stopgain") or ($genes->{$g}{variants}{$v}{clinvar} =~ /^CLINSIG=pathogenic/))
# 						{
# 							if ($genes->{$g}{variants}{$v}{dbsnp} eq 'rs11571833') {
# 								next;
# 							}
# 							$germPathogenic++;
# 						}
# 					}

# 				}
# 			}
# 		}
# 	}


# 	# print $html "<h2>Germline Variants</h2>\n";

# 	# print $html "<p><b>Germline Variants:</b> $germVars<br>\n";
# 	# print $html "<b>Germline Non-silent Variants:</b> $germNonsil<br>\n";
# 	# print $html "<b>Germline Non-silent Variants in HPCS Genes:</b> $germNonsilGenes<br>\n";
# 	# print $html "<b>Germline Rare Non-silent Variants in HPCS Genes:</b> $germNonsilGenesRare<br>\n";
# 	# print $html "<b>Germline Pathogenic Variants:</b> $germPathogenic<br>\n";
	
# 	# print $html "<table style=\"width:100%\" class=\"sortable\">\n";
# 	# print $html "<tr>\n";
# 	# print $html "<th>Gene</th><th>Variant</th><th>Copy Number</th><th>AB Count</th><th>Additional Information</th>";
# 	# print $html "</tr>\n";
	
# 	for my $g (@{ $germlineGeneOrder })
# 	{
# 		for my $v (sort keys %{ $genes->{$g}{variants} })
# 		{
# 			if ($genes->{$g}{variants}{$v}{mutation_class} =~ /germline/)
# 			{
# 				if (exists $nonsilent->{ $genes->{$g}{variants}{$v}{mutation_type} })
# 				{
# 					if ($genes->{$g}{variants}{$v}{rarity} ne "common")
# 					{
# 						printSimpleGeneRow($g, $v, $genes, $html, ""); ##AMY
# 					}
# 				}
# 			}
# 		}
# 	}
	
# 	print $html "</table>\n";

# }


# sub printClassifications
# {
# 	my $html = shift;
# 	my $data = shift;
# 	my $genes = shift;

# 	my $ploidy = "diploid ($data->{ploidy})";
# 	if ($data->{ploidy} > 2.3)
# 	{
# 		$ploidy = "polyploid ($data->{ploidy})";
# 	}

# 	my $HLAtypes = $data->{hla_types};
# 	$HLAtypes =~ s/\|/, /g;

# 	print $html "<h2>Classifications</h2>\n";

# 	print $html "<table style=\"width:100%\">\n";
# 	print $html "<tr>\n";
	print $html "<td><b>Alexandrov Class:</b> $data->{alexandrov_class}</td>\n";
	print $html "<td><b>Collisson Class:</b> $data->{collisson}</td>\n";
# 	print $html "<td><b>Ploidy:</b> $ploidy</td>\n";
# 	print $html "<tr>\n";

# 	print $html "<tr>\n";
	print $html "<td><b>Waddell Class:</b> $data->{waddell}</td>\n";
	print $html "<td><b>Moffitt Class:</b> $data->{moffitt}</td>\n";
# 	print $html "<td><b>Inferred Sex:</b> $data->{inferred_sex}</td>\n";
# 	print $html "</tr>\n";

# 	print $html "<tr>\n";
	print $html "<td colspan=\"3\"><b>HLA Types:</b> $HLAtypes</td>\n";
# 	print $html "</tr>\n";

# 	print $html "</table>\n";

# }

# sub printDsbrScore
# {
# 	# my $html = shift;
# 	# my $data = shift;
# 	# my $genes = shift;

# 	# my %niceName = (
# 	# 	"dsbr_snv_load" => "SNV Load >",
# 	# 	"dsbr_ct_ratio" => "SNV C>T Ratio <",
# 	# 	"dsbr_del4_load" => "4bp+ Deletion Load >",
# 	# 	"dsbr_del4_ratio" => "4bp+ Deletion Ratio >",
# 	# 	"dsbr_delsv_load" => "100-10kbp Deletion Load >",
# 	# 	"dsbr_delsv_ratio" => "100-10kbp Deletion Ratio >",
# 	# 	"dsbr_dup_load" => "10k-1mbp Duplication Load >",
# 	# 	"dsbr_sv_load" => "Structural Variant Load >",
# 	# 	"dsbr_first_genes" => "First Gene Hit(s)",
# 	# 	"dsbr_second_genes" => "Second Gene Hit(s)"
# 	# );

# 	# my %tableCell;

# 	$data->{dsbr_ct_ratio} = sprintf("%.02f", $data->{dsbr_ct_ratio});
# 	$data->{dsbr_del4_ratio} = sprintf("%.02f", $data->{dsbr_del4_ratio});
# 	$data->{dsbr_delsv_ratio} = sprintf("%.02f", $data->{dsbr_delsv_ratio});


# 	my $yesCircle = "<img src=\"$data->{common_dir}/DSBR-circle.png\" style=\"height:40%\">";
# 	my $noCircle = "<img src=\"$data->{common_dir}/Empty-circle.png\" style=\"height:40%\">";

# 	for my $type (qw/dsbr_snv_load dsbr_del4_load dsbr_del4_ratio dsbr_delsv_load dsbr_delsv_ratio dsbr_dup_load dsbr_sv_load/)
# 	{
# 		if ($data->{$type} > $data->{"${type}_cut"})
# 		{
# 			$tableCell{$type} = "$yesCircle <b>$niceName{$type} " . $data->{"${type}_cut"} . ":</b> $data->{$type}\n";
# 		}
# 		else
# 		{
# 			$tableCell{$type} = "$noCircle <b>$niceName{$type} " . $data->{"${type}_cut"} . ":</b> $data->{$type}\n";
# 		}
# 	}

# 	for my $type (qw/dsbr_ct_ratio/)
# 	{
# 		if ($data->{$type} < $data->{"${type}_cut"})
# 		{
# 			$tableCell{$type} = "$yesCircle <b>$niceName{$type} " . $data->{"${type}_cut"} . ":</b> $data->{$type}\n";
# 		}
# 		else
# 		{
# 			$tableCell{$type} = "$noCircle <b>$niceName{$type} " . $data->{"${type}_cut"} . ":</b> $data->{$type}\n";
# 		}
# 	}

# 	for my $type (qw/dsbr_first_genes dsbr_second_genes/)
# 	{
# 		if ($data->{$type} eq "")
# 		{
# 			$tableCell{$type} = "$noCircle <b>$niceName{$type}:</b> $data->{$type}";
# 		}
# 		else
# 		{
# 			$data->{$type} =~ s/\|/, /g;
# 			$tableCell{$type} = "$yesCircle <b>$niceName{$type}:</b> $data->{$type}";
# 		}
# 	}


# 	print $html "<h3>HR Deficiency Hallmark Score</h3>\n"; ##changed from DSBR Deficiency Hallmark Score

# 	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}-DSBR_score_bar.png\" style=\"width:100%\"></td>";

# 	print $html "<table style=\"width:100%\">\n";
# 	print $html "<tr>\n";
# 	print $html "<td>$tableCell{dsbr_snv_load}</td>\n<td>$tableCell{dsbr_delsv_load}</td>\n";
# 	print $html "</tr>\n";
# 	print $html "<tr>\n";
# 	print $html "<td>$tableCell{dsbr_ct_ratio}</td>\n<td>$tableCell{dsbr_delsv_ratio}</td>\n";
# 	print $html "</tr>\n";
# 	print $html "<tr>\n";
# 	print $html "<td>$tableCell{dsbr_del4_load}</td>\n<td>$tableCell{dsbr_dup_load}</td>\n";
# 	print $html "</tr>\n";
# 	print $html "<tr>\n";
# 	print $html "<td>$tableCell{dsbr_del4_ratio}</td>\n<td>$tableCell{dsbr_first_genes}</td>\n";
# 	print $html "</tr>\n";
# 	print $html "<tr>\n";
# 	print $html "<td>$tableCell{dsbr_sv_load}</td>\n<td>$tableCell{dsbr_second_genes}</td>\n";
# 	print $html "</tr>\n";

# 	print $html "</table>\n";

# }

# sub printMmrScore
# {
# 	# my $html = shift;
# 	# my $data = shift;
# 	# my $genes = shift;

# 	# my %niceName = (
# 	# 	"mmr_snv_load" => "SNV Load >",
# 	# 	"mmr_indel_load" => "Indel Load >",
# 	# 	"mmr_first_genes" => "First Gene Hit(s)",
# 	# 	"mmr_second_genes" => "Second Gene Hit(s)"
# 	# );

# 	# my %tableCell;

# 	# my $yesCircle = "<img src=\"$data->{common_dir}/MMR-circle.png\" style=\"height:40%\">";
# 	# my $noCircle = "<img src=\"$data->{common_dir}/Empty-circle.png\" style=\"height:40%\">";

# 	for my $type (qw/mmr_snv_load mmr_indel_load/)
# 	{
# 		if ($data->{$type} > $data->{"${type}_cut"})
# 		{
# 			$tableCell{$type} = "$yesCircle <b>$niceName{$type} " . $data->{"${type}_cut"} . ":</b> $data->{$type}\n";
# 		}
# 		else
# 		{
# 			$tableCell{$type} = "$noCircle <b>$niceName{$type} " . $data->{"${type}_cut"} . ":</b> $data->{$type}\n";
# 		}
# 	}

# 	for my $type (qw/mmr_first_genes mmr_second_genes/)
# 	{
# 		if ($data->{$type} eq "")
# 		{
# 			$tableCell{$type} = "$noCircle <b>$niceName{$type}:</b> $data->{$type}";
# 		}
# 		else
# 		{
# 			$data->{$type} =~ s/\|/, /g;
# 			$tableCell{$type} = "$yesCircle <b>$niceName{$type}:</b> $data->{$type}";
# 		}
# 	}


# 	print $html "<h3>MMR Deficiency Hallmark Score</h3>\n";

# 	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}-MMR_score_bar.png\" style=\"width:100%\"></td>";

# 	print $html "<table style=\"width:100%\">\n";
# 	print $html "<tr>\n";
# 	print $html "<td>$tableCell{mmr_snv_load}</td>\n<td>$tableCell{mmr_first_genes}</td>\n";
# 	print $html "</tr>\n";
# 	print $html "<tr>\n";
# 	print $html "<td>$tableCell{mmr_indel_load}</td>\n<td>$tableCell{mmr_second_genes}</td>\n";
# 	print $html "</tr>\n";
# 	print $html "</table>\n";

}




sub printMutationLoad
{
	my $html = shift;
	my $data = shift;

	# terribly hacky!
	my $snvText = `grep main $data->{out_path}/$data->{plot_dir}/$data->{tumour}-histogram.Rcode | grep SNVs | grep bp | sed 's/.*main="//' | sed 's/",.*//'`;
	my $indelText = `grep main $data->{out_path}/$data->{plot_dir}/$data->{tumour}-histogram.Rcode | grep Indel | grep bp | sed 's/.*main="//' | sed 's/",.*//'`;
	my $svText = `grep main $data->{out_path}/$data->{plot_dir}/$data->{tumour}-histogram.Rcode | grep SVs | grep bp | sed 's/.*main="//' | sed 's/",.*//'`;
#	my $neoText = `grep main $data->{out_path}/$data->{plot_dir}/$data->{tumour}-histogram.Rcode | grep Neo | grep bp | sed 's/.*main="//' | sed 's/",.*//'`;


	print $html "<h3>\tMutational Load</h3>\n";
	
	print $html "<table style=\"width:100%\">\n";
	print $html "<tr>\n";
	print $html "<td style=\"text-align:center\"><b>$snvText</b></td>";
	print $html "<td style=\"text-align:center\"><b>$indelText</b></td>";
	print $html "<td style=\"text-align:center\"><b>$svText</b></td>";
#	print $html "<td style=\"text-align:center\"><b>$neoText</b></td>";
#	print $html "<td style=\"text-align:center\"><b>$data->{snv_count} SNVs (%PR)</b></td>";
#	print $html "<td style=\"text-align:center\"><b>$data->{indel_count} Indels (%PR)</b></td>";
#	print $html "<td style=\"text-align:center\"><b>$data->{sv_count} SVs (%PR)</b></td>";
#	print $html "<td style=\"text-align:center\"><b>$data->{neo_antigens} Neo-antigens (%PR)</b></td>";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-histbox-snv_count-240.png\" style=\"width:100%\"></td>";
	print $html "<td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-histbox-indel_count-240.png\" style=\"width:100%\"></td>";
	print $html "<td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-histbox-sv_count-240.png\" style=\"width:100%\"></td>";
#	print $html "<td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-histbox-neo_antigens-240.png\" style=\"width:100%\"></td>";
	print $html "</tr>\n";
	print $html "</table>\n";
	
	print $html "<table style=\"width:100%\">\n";
	print $html "<tr>\n";
	print $html "<td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-snv_context-400x300.png\"></td>";
	print $html "<td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-indel_sv_bins-400x300.png\"></td>";
	print $html "</tr>\n";
	print $html "</table>\n";

	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}-snv_vaf-no_sample.png\">\n";
	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}-indel_vaf-no_sample.png\">\n";

}


sub printDriverGenes
{
	my $html = shift;
	my $genes = shift;
	my $driverGeneOrder = shift;
	my $nonsilent = shift;

	print $html "<h3>Driver Gene Status</h3>\n";
	
	print $html "<table style=\"width:100%\" class=\"sortable\">\n";
	print $html "<tr>\n";
	print $html "<th>Gene</th><th>Variant</th><th>Copy Number</th><th>AB Count</th><th>Additional Information</th>";
	print $html "</tr>\n";
	
	my $printed;
	
	for my $gene (@{ $driverGeneOrder })
	{
		$printed = 0;
		for my $v (sort keys %{ $genes->{$gene}{variants} })
		{
			unless ($genes->{$gene}{variants}{$v}{mutation_class} =~ /germline/)
			{
				unless ($genes->{$gene}{variants}{$v}{mutation_class} =~ /NA/)
				{
					printSimpleGeneRow($gene, $v, $genes, $html, "");
					$printed = 1;
				}
			}
		}
	
		for my $v (sort keys %{ $genes->{$gene}{variants} })
		{
			if ($printed == 0)
			{
				printSimpleGeneRow($gene, $v, $genes, $html, "gene only");
				$printed = 1;
			}
		}
	
	}
	print $html "</table>\n";
}

sub printWholeGenome
{
	my $html = shift;
	my $data = shift;

	print $html "<h3>Whole Genome Summary</h3>\n";
	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}.whole_genome.png\" style=\"width:100%\">\n";
}

sub printPloidyAndChromo
{
	my $html = shift;
	my $data = shift;

	print $html "<h3>Ploidy and Chromothripsis</h3>\n";
	print $html "<img src=\"$data->{plot_dir}/$data{tumour}-celluloid_contour.png\" style=\"width:100%\">\n";
}

sub printFullSomaticList
{
	my $html = shift;
	my $data = shift;
	my $genes = shift;
	my $nonsilent = shift;

	print $html "<h3>Full Deleterious Somatic Variant List</h3>\n";
	
	print $html "<table style=\"width:100%\">\n";
	print $html "<tr>\n";
	print $html "<td><b>Nonsynonymous-SNVs:</b> $data->{nonsyn_count}</td>\n";
	print $html "<td><b>Frameshift-deletions:</b> $data->{del_frameshift_count}</td>\n";
	print $html "<td><b>Deletion-Gene-BPs:</b> $data->{sv_del_bp_gene_count}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td><b>Stopgain-SNVs:</b> $data->{stopgain_count}</td>\n";
	print $html "<td><b>Nonframeshift-deletions:</b> $data->{del_nonframeshift_count}</td>\n";
	print $html "<td><b>Duplication-Gene-BPs:</b> $data->{sv_dup_bp_gene_count}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td><b>Stoploss-SNVs:</b> $data->{stoploss_count}</td>\n";
	print $html "<td><b>Frameshift-Insertions:</b> $data->{ins_frameshift_count}</td>\n";
	print $html "<td><b>Inversion-Gene-BPs</b> $data->{sv_inv_bp_gene_count}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td><b>Splicing-SNVs:</b> $data->{splice_count}</td>\n";
	print $html "<td><b>Nonframeshift-Insertions:</b> $data->{ins_nonframeshift_count}</td>\n";
	print $html "<td><b>Translocation-Gene-BPs:</b> $data->{sv_tra_bp_gene_count}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td><b>Total Deleterious SNVs:</b> " . ($data->{nonsyn_count} + $data->{stopgain_count} + $data->{stoploss_count} + $data->{splice_count}) . "</td>\n";
	print $html "<td><b>Total Deleterious Indels:</b> " . ($data->{del_frameshift_count} + $data->{del_nonframeshift_count} + $data->{ins_frameshift_count} + $data->{ins_nonframeshift_count}) . "</td>\n";
	print $html "<td><b>Total Deleterious SV Breakpoints:</b> " . ($data->{sv_del_bp_gene_count} + $data->{sv_dup_bp_gene_count} + $data->{sv_inv_bp_gene_count} + $data->{sv_tra_bp_gene_count}) . "</td>\n";
	print $html "</tr>\n";
	print $html "</table>\n";
	
	print $html "<br>\n";
	
	print $html "<table style=\"width:100%\" class=\"sortable\">\n";
	print $html "<tr>\n";
	print $html "<th>Gene</th><th>Variant</th><th>Copy Number</th><th>AB Count</th><th>Additional Information</th>";
	print $html "</tr>\n";
	
	for my $gene (sort keys %{ $genes })
	{
		for my $v (sort keys %{ $genes->{$gene}{variants} })
		{
			if (exists $nonsilent->{ $genes->{$gene}{variants}{$v}{mutation_type} })
			{
				unless ($genes->{$gene}{variants}{$v}{mutation_class} =~ /germline/)
				{
					printSimpleGeneRow($gene, $v, $genes, $html, "");
				} else {
#					printSimpleGeneRow($gene, $v, $genes, $html, "");
				}
			}
		}
	}
	print $html "</table>\n";

}


sub printWholeChroms
{
	my $html = shift;
	my $data = shift;

	print $html "<h3>Chromosome Summary</h3>\n";


	for my $chr (qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7-8 chr9-10 chr11-12 chr13-14 chr15-17 chr18-22 chrX-Y/)
	{
		print $html "<img src=\"$data->{plot_dir}/$data->{tumour}.whole_$chr.png\"><hr>\n";
	}
}





sub printFooter
{
	my $html = shift;
	my $data = shift;

	print $html "<div id=\"footer\">\n";
	print $html "<td style=\"width=33%;text-align:center;vertical-align:middle\"><img src=\"$data->{common_dir}/OICR_logo.jpg\"style=\"height:125px;\"></td>";
	print $html "<td style=\"width=33%;text-align:center;vertical-align:middle\"><img src=\"$data->{common_dir}/UHN_logo.png\" style=\"height:50px;\"></td>";
	print $html "<td style=\"width=33%;text-align:center;vertical-align:middle\"><img src=\"$data->{common_dir}/COMPASS_logo.png\" style=\"height:150px;\"></td>";
	print $html "</tr></table>\n";
	print $html "</div>\n";

}



sub getPercentRank
{
	my $value = shift;
	my $arrayRef = shift;

	my $rank = "";
	my $lowerCount = 0;
	my $equalCount = 0;

	for (my $i = 0; $i < scalar(@$arrayRef); $i++)
	{
		if ($value > $arrayRef->[$i])
		{
			$lowerCount++;
		}
		elsif ($value == $arrayRef->[$i])
		{
			$equalCount++;
		}
	}

	$rank = ($lowerCount + (0.5 * $equalCount)) / (scalar(@$arrayRef) + 1);

	print $rank . "\n";
	return sprintf("%.0f", ($rank * 100));
}


#printSimpleGeneRow($gene, $variant, \%vars, $html);
sub printSimpleGeneRow
{
	my $gene = shift;
	my $v = shift;
	my $vars = shift;
	my $html = shift;
	my $geneOnlyFlag = shift;

	my $geneOnly = 0;
	if ($geneOnlyFlag eq "gene only")
	{
		$geneOnly = 1;
	}

	my ($nuc,$aa,$cons);

	print $html "<tr>";
	print $html "<td>$gene</td>";

	if ($geneOnly == 0)
	{
		if ($vars->{$gene}{variants}{$v}{mutation_class} eq "somatic sv")
		{
			print $html "<td>$vars->{$gene}{variants}{$v}{mutation_type} ($vars->{$gene}{variants}{$v}{position})</td>";
		}
		elsif (($vars->{$gene}{variants}{$v}{mutation_class} eq "somatic snv") or ($vars->{$gene}{variants}{$v}{mutation_class} eq "somatic indel") or ($vars->{$gene}{variants}{$v}{mutation_class} eq "germline snp") or ($vars->{$gene}{variants}{$v}{mutation_class} eq "germline indel"))
		{
			$nuc = $vars->{$gene}{variants}{$v}{nuc_context};
			$nuc =~ s/\|.*$//;		# first isoform only
			$nuc =~ s/NM.*?://;		# drop the isoform name
	
			$aa = $vars->{$gene}{variants}{$v}{aa_context};
			$aa =~ s/\|.*$//;		# first isoform only
			$aa =~ s/NM.*?://;		# drop the isoform name
	
			print $html "<td>$vars->{$gene}{variants}{$v}{mutation_type} ($nuc;$aa)</td>";
		}
		else
		{
			print $html "<td>$vars->{$gene}{variants}{$v}{mutation_type} ($vars->{$gene}{variants}{$v}{position})</td>";
		}
	}
	else
	{
		print $html "<td></td>";
	}

	my $cn = 99999;
	my $cnCount = 0;
	my $cnFlag = "";
	my $ab;

	my @abs = split(/\|/, $vars->{$gene}{ab_counts});

	for my $state (split(/\|/, $vars->{$gene}{copy_number}))
	{
		unless ($state eq "NA")
		{
			if ($state < $cn)
			{
				$cn = sprintf("%0.3f", $state);
				$ab = $abs[$cnCount];
			}
		}
		$cnCount++;
	}

	if ($cnCount > 1)
	{
		$cnFlag = "*";
	}

	if ($cn == 99999)
	{
		$cn = "NA";
		$ab = "NA";
	}

	if ($cn eq "NA")
	{
		print $html "<td>$cn$cnFlag</td>";
	}
	elsif ($cn < 0.5)
	{
		print $html "<td class=\"cn_loss\">$cn$cnFlag</td>";
	}
	elsif ($cn >= ($vars->{$gene}{ploidy}*4))
	{
		print $html "<td class=\"cn_very_gain\">$cn$cnFlag</td>";
	}
	elsif ($cn < 1.5)
	{
		print $html "<td class=\"cn_loh\">$cn$cnFlag</td>";
	}
	elsif ($cn > ($vars->{$gene}{ploidy} + 1))
	{
		print $html "<td class=\"cn_gain\">$cn$cnFlag</td>";
	}
	else
	{
		print $html "<td>$cn$cnFlag</td>";
	}

	if ($ab =~ /^0\./)
	{
		print $html "<td class=\"cn_loh\">$ab$cnFlag</td>";
	}
	else
	{
		print $html "<td>$ab$cnFlag</td>";
	}


	my $info = "";

	if ($geneOnly == 0)
	{
		if (($vars->{$gene}{variants}{$v}{dbsnp} ne ".") and ($vars->{$gene}{variants}{$v}{dbsnp} ne "NA"))
		{
			$info .= "$vars->{$gene}{variants}{$v}{dbsnp}; ";
		}
		if (($vars->{$gene}{variants}{$v}{cosmic} ne ".") and ($vars->{$gene}{variants}{$v}{cosmic} ne "NA"))
		{
			$info .= "$vars->{$gene}{variants}{$v}{cosmic}; ";
		}
	
		if (($vars->{$gene}{variants}{$v}{rarity} ne ".") and ($vars->{$gene}{variants}{$v}{rarity} ne "NA"))
		{
			$info .= "$vars->{$gene}{variants}{$v}{rarity}; ";
		}
	
		my $clinvar = $vars->{$gene}{variants}{$v}{clinvar};
		if (($clinvar ne "NA") and ($clinvar ne "."))
		{
			$clinvar =~ s/CLINSIG=//;
			$clinvar =~ s/;.*//;
			$clinvar =~ s/\|.*//;
	
			$info .= "$clinvar; ";
		}
	
		if ($vars->{$gene}{variants}{$v}{cosmic_census_flag} eq "cosmic_mutation")
		{
			$info .= "COSMIC census hit; ";
		}
		$info =~ s/; $//;
	}

	print $html "<td>$info</td>";

	print $html "</tr>";
	


}




