#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use File::Temp qw/ tempfile tempdir /;
use Getopt::Std;


#my %mrna2spe=("^AT"=>"A.thaliana", '^\d\d'=>"A.lyrata",  "^Ca"=>"C.rubella", "^CA"=>"C.hirsuta","^ae"=>"A.arabicum", "^Br"=>"B.rapa", "^Th"=>"E.salsugineum", "^Tp"=>"S.parvula" );
#my @spes=("^AT", '^\d\d', "^Ca", "^CA", "^ae", "^Br", "^Th", "^Tp"); #keep the original order
####my @spes= keys %mrna2spe;   


my %opts = (e=>undef, s=>0, i=>undef);
getopts('e:s:i', \%opts);
if (@ARGV<2) {
	die "usage: pg_assotable.pl [-s 2] mcl.table info\n";
}

my $excl=$opts{e};
my $thr_size=$opts{s};

my $MCL=$ARGV[0];
my $info=$ARGV[1];
open my $hmcl, "< $MCL" or die "Can not find file $MCL\n";


open my $hinfo, "<$info" or die "Can not find find $info\n";
my %mrna2spe;
my @spes;
my %mrna2file;
while(<$hinfo>) {
   next if /^#/;
   chomp;
   my @t=split(/\t/);
   push(@spes, $t[0] );
   $mrna2spe{ $t[0] }   = $t[1];
   if ( @t > 2 )  {
      $mrna2file{ $t[0] }  = $t[2];
   } else {
      $mrna2file{ $t[0] }  = $t[1];
   }
}


my $tmp=0;
my %maps=map { $_ => $tmp++} @spes;

my $k=scalar( @spes );
#my @spes_fullname=map{ $mrna2spe{$_} } @spes;
#print join("\t", "id", @spes_fullname), "\n"; 
my $m=0;
while (<$hmcl>) {
    chomp;
    my @genefs;
    my @s=split(/[[:space:]]/);
    $m++;
    next if ($thr_size >1 && scalar(@s)< $thr_size ) ;

    my $num=0;
    my %taxa=();
    my @t=( 0 ) x $k ;
    foreach my $gene (@s) {
       foreach my $pm (keys %mrna2spe) {
           if ($gene =~/$pm/) {
	      push(@genefs,"$gene");
	      $taxa{$pm}=1;
	      $num++;
	      last;
           }
       }
    }
    if ($num>0) {
       my $ntaxa=scalar(keys %taxa);
       foreach my $ip (@genefs) {
           print "cluster_".$m, "\t",  $ip,  "\n";
       }
    }
}

print STDERR "$k species and $m clusters have been processed\n"; 

#foreach my $ky (keys %types ) {
#   print STDERR "$ky:",join(",", @{$types{$ky}}), " ";   
#}
#print STDERR "\n";
exit 0;

