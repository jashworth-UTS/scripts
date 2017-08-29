#!/usr/bin/perl -w

use strict;

if(@ARGV<1){die "give nucleotide sequence no spaces\n";}

my$inseq=$ARGV[0];
$inseq=~tr/A-Z/a-z/;
my@seq=split(//,$inseq);

my%cod = (ttt => "F", tct => "S", tat => "Y", tgt => "C", ttc => "F",
        tcc => "S", tac => "Y", tgc => "C", tta => "L", tca => "S",
        taa => ".", tga => ".", ttg => "L", tcg => "S", tag => ".",
        tgg => "W", ctt => "L", cct => "P", cat => "H", cgt => "R",
        ctc => "L", ccc => "P", cac => "H", cgc => "R", cta => "L",
        cca => "P", caa => "Q", cga => "R", ctg => "L", ccg => "P",
        cag => "Q", cgg => "R", att => "I", act => "T", aat => "N",
        agt => "S", atc => "I", acc => "T", aac => "N", agc => "S",
        ata => "I", aca => "T", aaa => "K", aga => "R", atg => "M",
        acg => "T", aag => "K", agg => "R", gtt => "V", gct => "A",
        gat => "D", ggt => "G", gtc => "V", gcc => "A", gac => "D",
        ggc => "G", gta => "V", gca => "A", gaa => "E", gga => "G",
        gtg => "V", gcg => "A", gag => "E", ggg => "G");

my@a=();
my@b=();
my@c=();

for(my$i=0;$i<@seq/3;++$i){
  $a[$i]=&codon_lookup(substr_array(\@seq,$i*3,3));
  unless($i> (@seq-3)/3-1 ){
    $b[$i]=&codon_lookup(substr_array(\@seq,$i*3+1,3));
    $c[$i]=&codon_lookup(substr_array(\@seq,$i*3+2,3));
  }
}

$"="";

print "\n\n@a\n\n@b\n\n@c\n\n";

sub substr_array{
 my($array,$start,$length)=@_;
 my$string="";
 for(my$i=$start;$i<$start+$length;++$i){
  $string=$string.$array->[$i];
 }
 return $string;
}

sub codon_lookup{
	my($codon)=@_;
	if($codon=~/n/){return " X ";}
	else{return $cod{$codon};}
}
