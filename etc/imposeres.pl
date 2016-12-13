#!/usr/bin/perl

BEGIN {
	die "Usage: $0 OPENRES SHUTRES\n" . $! if scalar @ARGV < 2;
	$last = 0;
	$state = "";
	$hdr = 1;
	@opens = ();
	@shuts = ();
	$openres = shift;
	$shutres = shift;
}

if ($hdr) {
	$hdr = 0;
	next;
}
($amp, $dur, $prop) = split(/\t/, $_, -1);
if ($amp != 0) {
	if ($dur > $openres) {
		$state = "open";
		push @opens, [ ($amp, $dur + $last) ];
		$last = 0;
	} else {
		if ($state ne "") {
			$last += $dur;
		}
	}
} else {
	if ($dur > $shutres) {
		$state = "shut";
		push @shuts, [ ($amp, $dur) ];
	}
}

END {
	print join("\t", @{ $_ }) foreach (@opens);
	print join("\t", @{ $_ }) foreach (@shuts);
}
