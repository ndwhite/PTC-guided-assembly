#!/usr/bin/perl
#Noor D. White Takes an interleaved fasta file and turns it into a regular one.

#Open file
open (FILE, $ARGV[0]) || die "where is the file?\n";

#Create outfile
unlink "Temp.txt";
open (TEMP, ">>Temp.txt") || die "Couldn't make tempfile\n";

#Process line, printing to temporary output file
until (eof FILE) {
	$line = <FILE>;
	chomp $line;

	if ($line =~ m/>/){
		print TEMP $sequence;
		$sequence = ();
		print TEMP "\n".$line."\n";
		} else {
		$sequence .= $line;
		}
	}
print TEMP $sequence;
close TEMP;

#Open tempfile, remove extra newline
open (FILE2, "Temp.txt");
$firstline = <FILE2>;
until (eof FILE2) {
	$line2 = <FILE2>;
	chomp $line2;
	print $line2."\n";
	}

unlink "Temp.txt";
