#!/usr/bin/perl
#Noor D. White
# Take 
# >CNGA1_Ggal_exon1
# ATGAAGATAGGAGTGATTGAGACCCATCACTCCCATACAGTTGTTCCCAGCGTGGTAGTGCATGACACCAGTGAGGACCCTGGACCGATGGACAAAGGGGGAAACAGATTTGCCAN
# >CNGA1_Ggal_exon2
# GCAATAGTATCTACCTGGTGCATTTGCACGCTACAATATTAACAACAATAGCAATAAAGATGA
#  
# to
# >CNGA1
# ATGAAGATAGGAGTGATTGAGACCCATCACTCCCATACAGTTGTTCCCAGCGTGGTAGTGCATGACACCAGTGAGGACCCTGGACCGATGGACAAAGGGGGAAACAGATTTGCCAN
# GCAATAGTATCTACCTGGTGCATTTGCACGCTACAATATTAACAACAATAGCAATAAAGATGA


open (FILE, $ARGV[0]) || die "where is the file?\n";
$gene = ();

$line1 = <FILE>;
chomp $line1;
$line1 =~ s/>//g;
$line2 = <FILE>;
chomp $line2;
@header = split ('_' , $line1);

print ">".$header[0]."\n".$line2."\n";

$gene = $header[0];

#Read in file
until (eof FILE) {
	$line1 = <FILE>;
	chomp $line1;
	$line1 =~ s/>//g;
	$line2 = <FILE>;
	chomp $line2;
	@header = split ('_' , $line1);
	
	if ($header[0] eq $gene){

		print $line2."\n";

		} else {

		print ">".$header[0]."\n".$line2."\n";
		$gene = $header[0];

	}

}

end;
