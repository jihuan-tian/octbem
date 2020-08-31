#!/usr/bin/perl

use strict;
use warnings;

my ($octfile, $matfile) = @ARGV;

die "Please provide an Octave file name: $!" if (!defined($octfile)) || (!(-f $octfile));

if (!defined($matfile)) {
    $matfile = $octfile . ".mlb";
    print "Output Matlab file has been created as $matfile\n";
}

if (! open OCTAVE, '<', $octfile) {
    die "Cannot open the Octave file $octfile: $!";
}

if (! open MATLAB, '>', $matfile) {
    die "Cannot open the Matlab file $matfile: $!";
}

# Process each line of the Octave script.
while(<OCTAVE>) {
    # Replace double comment symbols.
    s/##/%/g;
    
    # Replace single comment symbols.
    s/#/%/g;
    
    # Replace double quotes.
    s/"/'/g;

    if (!/^\s.%/) {
	# Ignore functions.
	if (/^\s*graphics_toolkit/) {
	    next;
	}

	if (/^\s*pkg load/) {
	    next;
	}

	# Replace logical negation operator ! with ~.
	s/!/~/g;
	
	# Replace graythresh function.
	s/(graythresh\s*\([a-zA-Z][a-zA-Z0-9_]*),\s*'.*'\s*\)/$1\)/g;

	# Replace GenVideos.sh script with GenVideos.bat.
	s/GenVideos\.sh/GenVideos\.bat/g;

	# Swap arguments of save function.
	s/(^\s*save\s*\(\s*)('-\w*')\s*,(\s*.*\s*)\)\s*;/$1$3, $2\);/g;

	# Replace ending keywords
	s/endfunction/end/g;
	s/endif/end/g;
	s/endswitch/end/g;
	s/endfor/end/g;
	s/endwhile/end/g;
    }

    print MATLAB;
}


close OCTAVE;
close MATLAB;
