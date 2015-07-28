#!/usr/bin/perl
use strict;
use lib '/home/sinhas/sinhas/lib/';
use Bio::SeqIO;

die "Usage: perl $0 enhancer_file.fasta Genomename_and_release(E.G. Dere2) output_prefix\n" unless @ARGV==3;

my $input_seqfile = $ARGV[0];
#my $orthmapdir = $ARGV[1];
#my $orthgenomedir = $ARGV[2];
my $genome_name = $ARGV[1];
my $genomes_base = "/shared-mounts/sinhas/wholegenome/Genomes/";
my $orthmapdir = "${genomes_base}/${genome_name}/maps";
my $orthgenomedir = "${genomes_base}/${genome_name}/masked";

my $out_prefix =$ARGV[2];



## reading the maps between mel and ortholog species
my %map=();
for my $chr ("2L","2R","3L","3R","X"){
	print STDERR "Reading $chr Map file\n";
	open (IN,"$orthmapdir/$chr.map") or die "Could not open map file";
	while (my $line=<IN>){
		chomp($line);
		my @tmp= split(/\s+/, $line);
		my $loc = shift(@tmp);
		$map{$chr}{$loc}=join(",",@tmp);
	}
	close(IN);
}


my $crmfile = $input_seqfile;
$crmfile =~ /([^\/]+)$/;
my $all_seqs = Bio::SeqIO->new(-file=>$crmfile, -format=>"Fasta");

open (OUT, ">${out_prefix}.fna");

while ( my $in_seq = $all_seqs->next_seq() ) { 
	
	my ($id, $seq)=GetOrthCRM($in_seq);
	my $baseid = (split / /, $id)[0];
	print OUT ">" . $id . "\n";
	print OUT "$seq\n";
}
close(OUT);

## return the ortholog crms
sub GetOrthCRM {
	my $f1 = shift;	

	my ($chr, $start, $stop)= split(/\s+/,$f1->desc());
	my $crmname = $f1->id();

	my $newstart = int($start/50)*50;
	my $newstop = int($stop/50)*50;
	## this case should not happen
	if(!defined($map{$chr}{$newstart}) || !defined($map{$chr}{$newstop})){ 
		print STDERR "Ignoring $crmname : No Mapping with Mel\n";
		return ($crmname,"");
        }
	## finding the first position that begin and end of ref crm has ortholog
	while ($map{$chr}{$newstart} =~ m/\?/ && ($newstart+50)<=$newstop) {
		$newstart+=50;
	}
        while ($map{$chr}{$newstop} =~ m/\?/ && $newstart<=($newstop-50)) {
                $newstop-=50;
        }
	if ($newstart eq $newstop){
		print STDERR "Ignoring $crmname: No ortholog crm found\n";
		return ($crmname,"");
	}
	## at this stage newstart and newstop have both none "?" value
	my $id=$crmname . " " . $map{$chr}{$newstart} . " " . $map{$chr}{$newstop};
	my $seq="";
	while($newstart+100<=$newstop){
		my ($bscafold, $bstart, $bdir)=split(/,/,$map{$chr}{$newstart});
		my ($escafold, $eend, $edir)=split(/,/,$map{$chr}{$newstart+100});
		if(($bscafold ne $escafold) or ($bdir ne $edir)){
			$newstart+=100;
			next;
		}
		if($bscafold =~ m/\?/){
			$newstart+=100;
			next;	
		}
		if(($bdir =~ m/\+/ and $bstart>=$eend) or ($bdir =~ m/\-/ and $bstart<=$eend)){
			$newstart+=100;
			next;
		}
		## if the length of the ortholog region is more than twice the length of region then ignore
		if(abs($bstart-$eend)>200){ 
			$newstart+=100;
			next;
		}
        #my $fo2 = Bio::SeqIO->new(-file=> "$orthgenomedir/$bscafold.fa", -format=>'Fasta');
        my $fo2 = Bio::SeqIO->new(-file=> "$orthgenomedir/$bscafold.TRMASKED", -format=>'Fasta');
		my $f2 = $fo2->next_seq();
		
		my $oseq="";
		if($bdir=~ /\+/){
			$oseq=substr($f2->seq(), $bstart, $eend-$bstart);
		}
		else{
			$oseq=substr($f2->seq(), $eend, $bstart-$eend);
			## do reverse complement here
			$oseq = reverse $oseq;
			$oseq =~ tr/[ACGTacgt]/[TGCAtgca]/;
		}
		$oseq =~ tr/[acgt]/[ACGT]/;
		$oseq =~ s/[^ACGT]/N/g;
		$seq.=$oseq;
		$newstart+=100;
	}
	return ($id, $seq);
}
