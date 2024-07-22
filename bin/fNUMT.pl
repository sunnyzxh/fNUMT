#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use File::Basename;

my ( $Help, $bam, $sample, $ref, $outdir );

GetOptions(
	'help|?'     => \$Help,
	'bam=s'      => \$bam,
	'sample=s'   => \$sample,
        'ref=s'      => \$ref,
	'outdir=s'   => \$outdir
);

die `pod2text $0` if ( $Help or !$bam or !$sample or !$outdir );

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Configuration
my $path_curf = File::Spec->rel2abs(__FILE__);
my ( $vol, $dirs, $file ) = File::Spec->splitpath($path_curf);
my $home = dirname($dirs);

my $samtools = "$home/bin/samtools";
my $cap3 = "$home/bin/cap3";
my $blastn = "$home/bin/blastn";
my $ref_chrm = "$home/lib/chrM.fa";
my $annovar = "$home/bin/annotate_variation.pl";
$ref ||= "hg19";

my $dir = "$outdir/$sample";
`mkdir -p $dir`;
`chmod u=rwx,g=xr,o=x $outdir/$sample`;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open IN, "$samtools view -F 3072 $bam |awk '\$3 != \"chrM\" && \$3 !~ /_/' |grep SA:Z:chrM |sort -k3,3 -k4,4n |" or die $!;

open OUT, ">$dir/tmp.reads";
open OUT1, ">$dir/$sample.numts.txt";
print OUT1 "Sample\tChr_n\tPos1\tPos2\tDistance\tn_region\tn_gene\tChr_mt\tStart\tEnd\tLength\tmt_region\tmt_gene\tSM_support\tMS_support\tSupported_readpairs\tMatched_breakpoints\t#Blastn\tSoft-clip_info\n";
open OUTM, ">$dir/$sample.nonref_NUMT_mutation.txt";


&showLog("Extracting soft-clip reads for $sample...");

my %hash; my %screads;
while (<IN>) {
    chomp;
    my @info = split;
    my $line = $_;
    my @sainfo = split /,/, $info[11];
    if ( $info[5] =~ /(\d+)S(\d+)M/ ) {
        if ( $sainfo[3] =~ /(\d+)S(\d+)M/ ) {
            my ( $pos1, $pos2 ) = ( $info[5] =~ /(\d+)S(\d+)M/ )[0, 1];
            my ( $sapos1, $sapos2 ) = ( $sainfo[3] =~ /(\d+)S(\d+)M/ )[0, 1];
            if ( abs ( $pos1 - $sapos2 ) <= 10 and abs ( $pos2 - $sapos1 ) <= 10 ) {
                my $softseq = substr($info[9], 0, $pos1);
                print OUT "$line\n";
                my $format = "SM:$info[2]:$info[3],chrM:$sainfo[1]";
                $hash{$format} ++;
                if ( !$screads{$format} ) {
                    $screads{$format} = ">$info[0]:SC\n$softseq";
                }
                else {
                    $screads{$format} = "$screads{$format}\n>$info[0]:SC\n$softseq";
                }
               
            }
        }
        else {
            print OUT "SMdiscordant\t$line\n";
        }
    }
    if ( $info[5] =~ /(\d+)M(\d+)S/ ) {
        if ( $sainfo[3] =~ /(\d+)M(\d+)S/ ) {
            my ( $pos1, $pos2 ) = ( $info[5] =~ /(\d+)M(\d+)S/ )[0, 1];
            my ( $sapos1, $sapos2 ) = ( $sainfo[3] =~ /(\d+)M(\d+)S/ )[0, 1];
            if ( abs ( $pos1 - $sapos2 ) <= 10 and abs ( $pos2 - $sapos1 ) <= 10 ) {
                my $softseq = substr($info[9], 150-$pos2, $pos2);
                print OUT "$line\n";
                $info[3] += $pos1 - 1;
                $sainfo[1] += $sapos1 - 1;
                my $format = "MS:$info[2]:$info[3],chrM:$sainfo[1]";
                $hash{$format} ++;
                if ( !$screads{$format} ) {
                    $screads{$format} = ">$info[0]:SC\n$softseq";
                }
                else {
                    $screads{$format} = "$screads{$format}\n>$info[0]:SC\n$softseq";
                }
            }
        }
        else {
            print OUT "MSdiscordant\t$line\n";
        }
    }
}

&showLog("Clustering, Assembling and NUMT-FPs detection...");
for my $k1 ( keys %hash ) {
    for my $k2 ( keys %hash ) {
        next if ( $k1 eq $k2 );
        #MS:chr12:132024304,chrM:5345
        my @info1 = split /[:,]/, $k1;
        my @info2 = split /[:,]/, $k2;
        next if ( !$hash{$k1} or !$hash{$k2} );
        if ( $info1[0] eq $info2[0] and $info1[1] eq $info2[1] and abs($info1[2] - $info2[2]) < 10 and abs($info1[4] - $info2[4]) < 10 ) {
            if ( $hash{$k1} >= $hash{$k2} ) {
                $hash{$k1} += $hash{$k2};
                $screads{$k1} = "$screads{$k1}\n$screads{$k2}";
                delete($hash{$k2});
            }
            else {
                $hash{$k2} += $hash{$k1};
                $screads{$k2} = "$screads{$k2}\n$screads{$k1}";
                delete($hash{$k1});
            }
        }
    }
}

my %dedup;
for my $k1 ( keys %hash ) {
    for my $k2 ( keys %hash ) {
        next if ( $k1 eq $k2 );
        my @tmp = ( $k1, $k2 );
        @tmp = sort(@tmp);
        my $seglen;
        #e.g., MS:chr11:49883569,chrM:61:3,SM:chr11:49883572,chrM:16089:3
        #e.g., MS:chr11:123964734,chrM:1388:5,SM:chr11:123964736,chrM:1321:4
        my @info1 = split /[:,]/, $tmp[0];
        my @info2 = split /[:,]/, $tmp[1];        
        if ( $info1[0] ne $info2[0] and $info1[1] eq $info2[1] and abs($info1[2] - $info2[2]) < 10 ) {
            $seglen = $info1[4] - $info2[4];
            if ( $seglen < 0 ) {
                $seglen = $info1[4] + 16569 - $info2[4];
            }
            
            if ( !$dedup{"$tmp[0]\t$tmp[1]"} ) {
                #print OUT1 "$tmp[0]\t$hash{$tmp[0]}\t$tmp[1]\t$hash{$tmp[1]}\t$seglen\n";
                my @nposes = ( $info1[2], $info2[2] ); @nposes = sort(@nposes);
                my @mtposes = ( $info2[4], $info1[4] );
                my $dis = $nposes[1] - $nposes[0];
                open OUT2, ">$dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq";
                print OUT2 "$screads{$tmp[0]}\n$screads{$tmp[1]}\n";
                #print "$tmp[0]\t$tmp[1]\n";                
                my $support_readpairs = 0;
                open IND, "$samtools view -F 3072 $bam |awk '\$3 == \"chrM\" && \$7 eq \"$info1[1]\" && \$8 > $nposes[0]-1000 && \$8 < $nposes[1]+1000' |" or die $!;
                while (<IND>) {
                    chomp;
                    $support_readpairs ++;
                    my @info = split;
                    my $mtseq = "";
                    if ( $info[5] eq "150M" ) {
                        $mtseq = $info[9];
                    }
                    elsif ( $info[5] =~ /(\d+)S(\d+)M/ ) {
                        my ($p1, $p2) = ( $info[5] =~ /(\d+)S(\d+)M/ )[0, 1];
                        $mtseq = substr($info[9], $p1, $p2);
                    }
                    elsif ( $info[5] =~ /(\d+)M(\d+)S/ ) {
                        my $p = ( $info[5] =~ /(\d+)M(\d+)S/ )[0];
                        $mtseq = substr($info[9], 0, $p);
                    }
                    print OUT2 ">$info[0]:SC\n$mtseq\n" if ( length($mtseq) > 2 );
                }
                close IND;
                close OUT2;
                if ( $support_readpairs > 200 ) {
                    `rm -f $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq`;
                    next;
                }
                #print OUT1 "$sample\t$info1[1]\t$nposes[0]\t$nposes[1]\t$dis\tchrM\t$mtposes[0]\t$mtposes[1]\t$seglen\t$hash{$tmp[0]}\t$hash{$tmp[1]}\t$support_readpairs\n";

                `$cap3 $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq > $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq.contig`;
                my $contig = `less $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq.cap.contigs |wc -l`;
                chomp $contig;
                #print "$contig\t$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1]\n";
                if ( $contig == 0 ) {
                    `rm -f $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq*`;
                    next;
                }
                `$blastn -query $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq.cap.contigs -db $ref_chrm -out $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq.cap.contigs.blast.mt7 -outfmt 7`;
                `$blastn -query $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq.cap.contigs -db $ref_chrm -out $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq.cap.contigs.blast.mt3 -outfmt 3`;
                
                open INBLA, "$dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq.cap.contigs.blast.mt7";
                my @blasts = (); my $flag1 = 0; my $flag2 = 0;
                while (<INBLA>) {
                    chomp;
                    next if /\#/;
                    my @blastinfo = split;
                    push @blasts, "$blastinfo[8]-$blastinfo[9]-len$blastinfo[3]-mis$blastinfo[4]-gap$blastinfo[5]";
                    if ( abs($blastinfo[8] - $mtposes[0]) < 5 or abs($blastinfo[9] - $mtposes[0]) < 5 ) { $flag1 = 1; }
                    if ( abs($blastinfo[8] - $mtposes[1]) < 5 or abs($blastinfo[9] - $mtposes[1]) < 5 ) { $flag2 = 1; }
                }
                my $flag = $flag1 + $flag2;
                my $mapping = join ",", @blasts;
                if ( ($hash{$tmp[0]} <= 1 and $hash{$tmp[1]} <= 1) or ($flag == 0) or ($#blasts >= 5) ) {
                    `rm -f $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq*`;
                    next;
                }
                `echo "$info1[1]\t$nposes[0]\t$nposes[1]\t0\t0\t0" > $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].n.anno`;
                if ( $mtposes[0] < $mtposes[1] ) { `echo "MT\t$mtposes[0]\t$mtposes[1]\t0\t0\t0" > $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].mt.anno`; }
                else { `echo "MT\t1\t$mtposes[1]\t0\t0\t0" > $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].mt.anno`; 
                       `echo "MT\t$mtposes[0]\t16569\t0\t0\t0" >> $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].mt.anno`; }

                `perl $annovar --geneanno --buildver $ref $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].n.anno $home/lib`;
                my $nanno = `less $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].n.anno.variant_function |tail -1`; chomp $nanno;
                my ( $n_region, $n_gene ) = ( split /\s+/, $nanno )[0, 1];
                `rm -f $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].n.anno*`;

                `perl $annovar -buildver GRCh37_MT -dbtype ensGene $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].mt.anno $home/lib`;
                my ( $mt_region, $mt_gene );
                if ( $mtposes[0] < $mtposes[1] ) { 
                    my $mtanno = `less $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].mt.anno.variant_function |tail -1`; 
                    chomp $mtanno;
                    ( $mt_region, $mt_gene ) = ( split /\s+/, $mtanno )[0, 1];
                }
                else {
                    my $mtanno = `less $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].mt.anno.variant_function |tail -2`;
                    chomp $mtanno;
                    my @mtannos = split /\n/, $mtanno;
                    my ( $mt_region1, $mt_gene1 ) = ( split /\s+/, $mtannos[0] )[0, 1];
                    my ( $mt_region2, $mt_gene2 ) = ( split /\s+/, $mtannos[1] )[0, 1];
                    $mt_region = "$mt_region1-$mt_region2";
                    $mt_gene = "$mt_gene1-$mt_gene2";
                }
                `rm -f $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].mt.anno*`;
                
                print OUT1 "$sample\t$info1[1]\t$nposes[0]\t$nposes[1]\t$dis\t$n_region\t$n_gene\tchrM\t$mtposes[0]\t$mtposes[1]\t$seglen\t$mt_region\t$mt_gene\t$hash{$tmp[0]}\t$hash{$tmp[1]}\t$support_readpairs\t$flag\t$mapping\t$tmp[0]:$hash{$tmp[0]},$tmp[1]:$hash{$tmp[1]}\n";
                close INBLA;
                    
                my $mutations = `perl $home/bin/get_mutation.pl $dir/$info1[1]-$nposes[0]-$nposes[1]-chrM-$mtposes[0]-$mtposes[1].seq.for.cap3.fastq.cap.contigs.blast.mt3`;
                chomp $mutations;
                print OUTM "$mutations\n" if ( length($mutations)>1 ); 
            }
            $dedup{"$tmp[0]\t$tmp[1]"} = 1;
        }
    }
}
close IN;

&showLog("DONE SUCCESSFULLY!");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub showLog {
    my @t = localtime();
    printf STDERR "[%04d-%02d-%02d %02d:%02d:%02d]\t%s\n", $t[5] + 1900,
     $t[4] + 1, @t[ 3, 2, 1, 0 ], $_[0];
}


=head1 Function
    Detecting non-ref NUMTs and the derived NUMT-FPs from short-read WGS data

=head1 Usage
    perl fNUMT.pl -b example.bam -s sample.name -o outdir

=head1 Options
    -h|-help            help
    -b|-bam      [s]    bam file
    -s|-sample   [s]    sample name
    -r|-ref      [s]    reference version [hg19/hg38]
    -o|-outdir   [s]    output dir

=head1 Author
    sunnyzxh@connect.hku.hk

=head1 Version
    v1.1; 2022.10.28

=cut
