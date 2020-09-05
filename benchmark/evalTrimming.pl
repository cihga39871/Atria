#!/usr/bin/env perl

# This code is part of Skewer (https://sourceforge.net/projects/skewer/). The License:
#
# The MIT License (MIT)
#
# Copyright (c) 2013-2014 by Hongshan Jiang
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

use strict;

if(@ARGV < 3){
	print STDERR "Usage: $0 fullLen lengths.tab file1.fastq [file2.fastq [lengths2.tab]] > summary\n";
	exit(1);
}

my ($full_len, $tab_file, $file1, $file2, $tab_file2) = @ARGV;
our ($tp, $tn, $fp, $fp2, $fn, $fn2) = (0, 0, 0, 0, 0, 0);

open(TAB, "<$tab_file") or die("Can not open $tab_file for reading\n");

if(!open(IN, "<$file1")){
	close TAB;
	die("Can not open $file1 for reading\n");
}

my $id = &calcMetrics(\*TAB, \*IN);
die("read $id is in $file1 but not in $tab_file\n") if(defined $id);

close IN;
if(defined $file2){
	my $bUseTab2;
	if(defined $tab_file2){
		if(open(TAB2, "<$tab_file2")){
			$bUseTab2 = 1;
		}
		else{
			print STDERR "Warning: Can not $tab_file2 for reading, using $tab_file instead\n";
			$bUseTab2 = 0;
		}
	}
	else{
		$bUseTab2 = 0;
	}
	#
	if(open(IN, "<$file2")){
		my $fh = $bUseTab2 ? \*TAB2 : \*TAB;
		my $fname = $bUseTab2 ? $tab_file2 : $tab_file;
	    $id = &calcMetrics($fh, \*IN);
		die("read $id is in $file2 but not in $fname\n") if(defined $id);
	}
	else{
		print STDERR "Warning: Can not $file2 for reading, using information in $file1 only\n";
	}
	if($bUseTab2){
		close TAB2;
	}
}

close TAB;

my $ppv = ($tp+$fp+$fp2+$fn2) > 0 ? $tp/($tp+$fp+$fp2+$fn2) : 0;
my $sen = ($tp+$fn+$fp2+$fn2) > 0 ? $tp/($tp+$fn+$fp2+$fn2) : 0;
my $spec = $tn/($tn+$fp);
my $dom = sqrt(($tp+($fp+$fp2))*($tp+($fn+$fn2))*($tn+($fp+$fp2))*($tn+($fn+$fn2)));
my $cc = ($dom > 0) ? (($tp * $tn - ($fp+$fp2) * ($fn+$fn2)) / $dom) : 0;
my $fpr = (1 - $spec);
print "TP\tFP_ft\tFP_ot\tFN_fr\tFN_ut\tTN\tPPV\tSen.\tSpec.\tmCC\n";
print "$tp\t$fp\t$fp2\t$fn\t$fn2\t$tn\t$ppv\t$sen\t$spec\t$cc\n";
print "(FPR, TPR) = ($fpr, $sen)\n";

exit(0);

sub calcMetrics
{
    my ($fh_tab, $fh_file) = @_;
	our ($tp, $tn, $fp, $fp2, $fn, $fn2);

    my $line;
    my ($id, $len);
    my ($id2, $len2, $seq);
    while($line = <$fh_tab>){
        chomp($line);
       ($id, $len) = split(/\t/, $line);

       $id2 = <$fh_file>; chomp($id2);
       $seq = <$fh_file>; chomp($seq);
       <$fh_file>; <$fh_file>;
       ($id2) = split(/\//, substr($id2,1));
       while($id2 ne $id){
           if($len == 0){
               $tp++;
           }
           else{
			   if($len == $full_len){
	               $fp++;
			   }
			   else{
				   $fp2++;
			   }
           }
		   if(!($line=<$fh_tab>)){
			   return $id2;
	           #die("read $id2 is $file1 but not in $tab_file\n");
		   }
           chomp($line);
           ($id, $len) = split(/\t/, $line);
       }
       $len2 = length($seq);
       if($len == $len2){
           if($len == $full_len){
               $tn++;
           }
           else{
               $tp++;
           }
       }
       else{ # $len != $len2
		   if($len < $len2){
			   if($len2 == $full_len){
				   $fn++;
			   }
			   else{
				   $fn2++;
			   }
		   }
		   else{ # $len > $len2
			   if($len == $full_len){
				   $fp++;
			   }
			   else{
				   $fp2++;
			   }
		   }
       }
    }
	return undef;
}
