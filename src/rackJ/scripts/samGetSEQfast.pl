#!/usr/bin/perl -w
use strict;

my $usage   = "Usage: samGetSEQfast.pl [-unmap] <inSAM> [<FASTA>]+\n";

#Retrieve optional parameter
my $unmap=0;
my @arg_idx=(0..@ARGV-1);
for my $i (0..@ARGV-1) {
	if ($ARGV[$i] eq '-unmap') {
		$unmap=1;
		delete $arg_idx[$i];
	}
}
my @new_arg;
for (@arg_idx) { push(@new_arg,$ARGV[$_]) if (defined($_)); }
@ARGV=@new_arg;


# SAM
my $samfile = shift or die $usage;
die $usage if @ARGV < 1;

$samfile = "-" if $samfile eq "STDIN";
my $samFILEHANDLE;
open($samFILEHANDLE, $samfile) or die("Error: cannot open sam file '$samfile'\n");

# get first line in inSAM
my $currentSamRec=getCurrentSAMRec($samFILEHANDLE);
if (@$currentSamRec){
	# skip header lines
	while($$currentSamRec[0] =~ /^@/){
		print join("\t",@$currentSamRec)."\n";
		$currentSamRec=getCurrentSAMRec($samFILEHANDLE);
	}
}else{
	exit;
}

# parse FASTA
for my $i ( 0 ... @ARGV - 1 ) {
	my $SEQId="";
	my $SEQ="";
	open(my $FASTAstream, $ARGV[$i]) or die("Error: cannot open FASTA file '$ARGV[$i]'\n");
	while (my $line = <$FASTAstream>) {
		if ($line =~ /^\>/){
			if($SEQId ne ""){
				if (@$currentSamRec && $SEQId eq $$currentSamRec[0]) {
					samFill($SEQId, $SEQ, $samFILEHANDLE, \$currentSamRec);
				}else{
					outSeqUnalign($SEQId, $SEQ) if ($unmap);
				}
			}
			
			$SEQId = trim(substr($line,1));
			$SEQ = "";
		}else{    
			$SEQ .= trim($line);
		}
	}

	# handle the last sequence
	if (@$currentSamRec && $SEQId eq $$currentSamRec[0]) {
		samFill($SEQId, $SEQ, $samFILEHANDLE, \$currentSamRec);
	}else{
		outSeqUnalign($SEQId, $SEQ) if ($unmap);
	}

	close $FASTAstream;
}

close $samFILEHANDLE;


#-------------------------------
# subroution
#-------------------------------
sub getCurrentSAMRec{
	my $fh = shift;
	
	my @ss;
	if (defined (my $sline = <$fh>)){
		@ss = split("\t",trim($sline));
	}
	
	return \@ss;
}

sub samFill{	
	my ($seqId,$sequence,$fh,$sam) = @_;
	
	while(@{$$sam} && $seqId eq ${$$sam}[0]){
		${$$sam}[9] = $sequence;
		${$$sam}[9] =~ tr/ATCG/TAGC/ if ${$$sam}[1] & 16;
		${$$sam}[9] = reverse ${$$sam}[9] if ${$$sam}[1] & 16;
		${$$sam}[5] =~ tr/H/S/;
		print join( "\t", @{$$sam} )."\n";

		$$sam = getCurrentSAMRec($fh);
	}
}

sub outSeqUnalign{	
	my ($seqId,$sequence) = @_;
	
	my @unalign = ($seqId,4,"*",0,255,"*","*",0,0,$sequence,"*");
	print join( "\t", @unalign )."\n";
}


sub trim {
	my $str = shift;
	$str =~ s/\s+$//g;
	$str =~ s/^\s+//g;
	return $str;
}