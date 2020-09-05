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

if(@ARGV != 1){
    print STDERR "Usage: $0 file.aln > lengths.tab\n";
    exit(1);
}
my ($aln_file) = @ARGV;

my ($line, $no);
my @columns;
my ($id, $len, $len2);
my @chars;
$no = -1;
open(ALN, "<$aln_file") or die("Can not open $aln_file for reading\n");
while($line = <ALN>){
    chomp($line);
    if($line =~ /^>/){
        @columns = split(/\t/, $line);
        @columns = split(/\//, $columns[1]);
        $id = $columns[0];
        $no = 0;
        next;
    }
	next if($no < 0);
    $no++;
    if($no == 1){ # first sequence
        $len = length($line);
		next;
    }
	if($no == 2){ # second sequence
		@chars = split(//, substr($line,0,$len));
		my $del=0;
		for(my $i=$#chars; $i>=0; $i--){
			if($chars[$i] eq '-'){
				$del++;
			}
		}
		$len -= $del;
        print "$id\t$len\n";
	}
}
close ALN;

exit(0);
