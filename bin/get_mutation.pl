#!/usr/bin/perl -w
use strict;

my ( $blast3 ) = @ARGV;

my %reverse = ('A'=>'T', 'T'=>'A', 'G'=>'C', 'C'=>'G', 'a' => 't', 't' => 'a', 'g' => 'c', 'c' => 'g');
open IN, "less $blast3 |awk '(\$1 ~ /Query/ || \$1==0 ) && \$4' |" or die $!;
#open IN, "less contig3.blast |awk '(\$1 ~ /Query/ || \$1==0 ) && \$4' |" or die $!;

my %snp; my %inses; my @dels = ();
while (<IN>) {
    chomp;
    my $que = $_;
    my $db = <IN>;
    chomp $db;
    my @dbs = split /\s+/, $db;
    my @ques = split /\s+/, $que;
    if ( abs($dbs[1] - $dbs[3]) < 5 and length($ques[2]) > 10 ) {
        print "$db\n";
        $db = <IN>;
        chomp $db;
        @dbs = split /\s+/, $db;
    }

    my @refbase = split //, $dbs[2];
    my @quebase = split //, $ques[2];

    my $format = "";
    my $mutpos;
    for my $i ( 0..$#refbase ) {
        if ( $refbase[$i] ne "." ) {
            if ( $dbs[1] < $dbs[3] ) {
                my $sub = substr($dbs[2], 0, $i);
                $sub =~ s/-//g;
                $mutpos = $dbs[1] + length($sub);
            }
            else {
                my $sub = substr($dbs[2], $i, length($dbs[2]) - $i);
                $sub =~ s/-//g;
                $mutpos = $dbs[3] + length($sub) - 1;
            }
            
            if ( $refbase[$i] =~ /[ATGCatgc]/ and $quebase[$i] =~ /[ATGCatgc]/ ) {
                if ( $dbs[1] > $dbs[3] ) {
                    $refbase[$i] = $reverse{$refbase[$i]};
                    $quebase[$i] = $reverse{$quebase[$i]};
                }
                $format = "NC_012920.1:m.$mutpos"."$refbase[$i]".">$quebase[$i]";
                $snp{$format} = 1;
            }
            elsif ( $refbase[$i] eq "-" and $quebase[$i] =~ /[ATGCatgc]/ ) {
                my $start = $mutpos - 1;
                my $end = $mutpos + 1;
                my $pformat = "NC_012920.1:m.$start"."_"."$mutpos"."ins";
                if ( $dbs[1] > $dbs[3] ) {
                   $pformat = "NC_012920.1:m.$mutpos"."_"."$end"."ins";
                }
                if ( !$inses{$pformat} ) {
                    $quebase[$i] = $reverse{$quebase[$i]} if ( $dbs[1] > $dbs[3] );
                    $inses{$pformat} = "$quebase[$i]";
                }
                else {
                    if ( $dbs[1] > $dbs[3] ) { $inses{$pformat} = "$reverse{$quebase[$i]}"."$inses{$pformat}"; }
                    else { $inses{$pformat} = "$inses{$pformat}"."$quebase[$i]"; }
                }
            }
            elsif ( $refbase[$i] =~ /[ATGC]/ and $quebase[$i] eq "-" ) {
            #    $format = "NC_012920.1:m.$mutpos"."_"."$mutpos"."del";
                push @dels, $mutpos;
            }
            
        }
    }    
}

for my $s ( sort keys %snp ) {
    print "$s\n";
}
for my $i ( sort keys %inses ) {
    my $format = "$i"."$inses{$i}";
    print "$format\n";
}

if ( $#dels >= 0 ) {
@dels = sort {$a<=>$b} @dels;
#print "@dels\n";
my ( $start, $end ) = ( $dels[0], $dels[0] );

for my $i ( 1..$#dels ) {
    if ( ($dels[$i] - $end) == 1 ) {
        $end = $dels[$i];
    }
    elsif ( ($dels[$i] - $end) > 1 ) {
        print "NC_012920.1:m.$start"."_"."$end"."del\n";
        $start = $dels[$i];
        $end = $dels[$i];
    }
}
print "NC_012920.1:m.$start"."_"."$end"."del\n";
}
        
close IN;
