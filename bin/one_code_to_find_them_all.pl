#!/usr/bin/perl

#     One code to find them all  -- perl utility to extract information from RepeatMasker output files
#     Copyright (C) 2014  Bailly-Bechet Marc
#
#     Edit: Jul-06-2022 by Clément Goubert for GraffiTE pipeline (https://github.com/cgroza/GraffiTE)
#           - Remove "LTR/" added when Unknown (or "Unspecified") TEs are in the RM output
#
#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see <http://www.gnu.org/licenses/>.


### This script takes an output from Repeat Masker
### and try to find transposons inside
### It reports the list of transposons found
### and some quantitative information about them
### It can also retrieve fasta sequences from a local file

use FileHandle;
use Getopt::Long;
use File::Find;
use File::Basename;

### Default value of choice and strict
$choice='';
$strict='';
$dry="";
$length_file="";
$unknown="";
$fasta_file="no_value";
$flanking=0;
$insert=-1;
$help="";
%sequence=();

GetOptions('rm=s' => \$rm_file,
        'ltr=s' => \$ltr_file,
        'length=s' =>\$length_file,
        'fasta:s' =>\$fasta_file,
        'flanking=i' =>\$flanking,
        'insert=i'=>\$insert,
        'strict' => \$strict,
        'choice' => \$choice,
        'unknown'=>\$unknown,
        'dry-run' =>\$dry,
        'help'=>\$help);

if($help==1){
  warn "#################\n";
  warn "Usage:\none_code_to_find_them_all.pl --rm <file> --ltr <file2> [--length <file3>] [--fasta <file4> [--flanking <num>]] [--insert <num2>] [--strict] [--choice] [--unknown] [--dry-run]\n";
  warn "Option --rm requires a file or directory to be passed  in order to work correctly;\n";
  warn "Options  --ltr and  --length require a file to be passed in order to work correctly;\n";
  warn "Option --fasta can work without any file given; if a file is given it must have the .fa extension\n";
  warn "Options --flanking and --insert require a positive integer value to work correctly;\n";
  warn "#################\n";
  warn "See README file for more details concerning the different options usage\n";
  exit(0);
}

### Sanity checks
if($rm_file eq ""){die("Dying! Option --rm must be passed\n")}
if($ltr_file eq ""){die("Dying! Option --ltr must be passed\n")}
if($dry==1){$choice=1;$global_choice_counter=0};
if($flanking <0){
  warn "Given flanking value was negative; flanking value replaced by 0\n";
  $flanking=0;
}

### Checking if we work on a directory or a single file either
if(-d $rm_file){$file_ok=0}else{$file_ok=1}

### Subroutines min and max
sub min{
  my $v1=$_[0];
  my $v2=$_[1];
  my $min=$v1;
  if($v2<$min){$min=$v2;}
  return($min);
}

sub max{
  my $v1=$_[0];
  my $v2=$_[1];
  my $max=$v1;
  if($v2>$max){$max=$v2;}
  return($max);
}

#### The subroutine combining two lines into one element
sub combine_lines{
  my @now=@{$_[0]};
  my @second=@{$_[1]};
  my $match_length_now=$now[6]-$now[5]+1;
  my $match_length_second=$second[6]-$second[5]+1;
  my $combined=$now[0]."/".$second[0]."\t";
  my $id=($now[1]*$match_length_now+$second[1]*$match_length_second)/($match_length_now+$match_length_second);
  $id=sprintf("%.3f",$id);
  my $idins=($now[2]*$match_length_now+$second[2]*$match_length_second)/($match_length_now+$match_length_second);
  $idins=sprintf("%.3f",$idins);
  my $iddel=($now[3]*$match_length_now+$second[3]*$match_length_second)/($match_length_now+$match_length_second);
  $iddel=sprintf("%.3f",$iddel);
  $combined.=$id."\t".$idins."\t".$iddel."\t".$now[4]."\t".$now[5]."\t".$second[6]."\t";
  my $cover=$now[7]+$second[7];
  if(($now[8] eq "+") & ($now[12] >= $second[11]) ){
    $cover-=($now[12]-$second[11]+1);
  }
  if(($now[8] eq "C") & ($second[12] >= $now[13]) ){
    $cover-=($second[12]-$now[13]+1);
  }
  $combined.=$cover."\t";
  $combined.=$now[8]."\t".$now[9]."\t".$now[10]."\t";
  $combined.=$now[11]."\t";
  if($now[8] eq "+"){
    $combined.=$second[12];
  }else{
    $combined.=$now[12];
  }
  $combined.="\t".$second[13]."\t".$now[14]."/".$second[14]."\t";
  $sum=$now[15]+$second[15];
  $combined.=$sum;
  return($combined);
}


### Subroutine associating LTR and internal parts into one copy
sub build_ltr{
  my @now=@{$_[0]};
  my @second=@{$_[1]};
  my $match_length_now=$now[6]-$now[5]+1;
  my $match_length_second=$second[6]-$second[5]+1;
  my $combined=$now[0]."/".$second[0]."\t";
  my $id=($now[1]*$match_length_now+$second[1]*$match_length_second)/($match_length_now+$match_length_second);
  $id=sprintf("%.3f",$id);
  my $idins=($now[2]*$match_length_now+$second[2]*$match_length_second)/($match_length_now+$match_length_second);
  $idins=sprintf("%.3f",$idins);
  my $iddel=($now[3]*$match_length_now+$second[3]*$match_length_second)/($match_length_now+$match_length_second);
  $iddel=sprintf("%.3f",$iddel);
  $combined.=$id."\t".$idins."\t".$iddel."\t".$now[4]."\t".$now[5]."\t".$second[6]."\t";
  my $cover=$now[7]+$second[7];
  $combined.=$cover."\t";
  $combined.=$now[8]."\t".$now[9]."\t".$now[10]."\t";
  $combined.="NA\tNA\tNA\t";
  $combined.=$now[14]."/".$second[14]."\t";
  $sum=$now[15]+$second[15];
  $combined.=$sum;
  return($combined);
}


### Subroutine for position-based sorting
sub by_position{
  @c=split /\s+/,$a;
  @d=split /\s+/,$b;
  {$c[5] <=>$d[5]};
}

### Subroutine retrieving fasta sequences
sub get_fasta_seq{
  my $t=$_[0];
  my $mem=$_[1];
  my @now=split/\t/,$t;
  my $contig=$now[4];
  my $sense=$now[8];
  my $fam=$now[10];
  my @lines= split /\n/,$mem;
  my @lims=();
  @F=split /\t/,$lines[0];
  $lims[0]=&max($F[5]-$flanking,1);
  my $memelem=$F[9];
  my $local_end=&min($F[6]+$flanking,$contig_length{$contig});
  my $elem = $memelem;
  my $pos="";
  my $seq="";
  my $tmpseq="";

  shift @lines;
  $i=0;
  foreach $l (@lines){
    @F=split /\t/,$l;
    if($F[9] ne $memelem){
      $elem.=":".$F[9];
      $memelem=$F[9];
    }
    if($F[5]-$flanking>$local_end+1){
      push @lims,$local_end;
      push @lims,&max($F[5]-$flanking,1);
    }
    $local_end=&min($F[6]+$flanking,$contig_length{$contig});
  }
  push @lims,$local_end;

### BLOCK ERASED DUE TO ERROR IN REPORTING REVERSE COMPLEMENT FASTA SEQUENCES
  #   if($sense eq "C"){
#     $pos.=")|";
#     for($i=$#lims;$i>0;$i=$i-2){
#       if($i<$#lims){$pos=",".$pos;}
#       $pos=$lims[$i-1]."..".$lims[$i].$pos;
#       $tmpseq.=substr($sequence{$contig},$lims[$i-1]-1,$lims[$i]-$lims[$i-1]+1);
#     }
#     $tmpseq=~tr/("a","c","g","t","A","C","G","T")/("t","g","c","a","T","G","C","A")/;
#     $seq=reverse($tmpseq);
#     $pos="|c(".$pos;
#   }else{
#     $pos.="|(";
#     for($i=0;$i<$#lims;$i=$i+2){
#       if($i>0){$pos.=","}
#       $pos.=$lims[$i]."..".$lims[$i+1];
#       $seq.=substr($sequence{$contig},$lims[$i]-1,$lims[$i+1]-$lims[$i]+1);
#     }
#     $pos.=")|";
#   }
### END OF BLOCK ERASED

### CORRECTION TO GET CORRECT REVERSE COMPLEMENT FASTA SEQUENCES
for($i=0;$i<$#lims;$i=$i+2){
      if($i>0){$pos.=","}
      $pos.=$lims[$i]."..".$lims[$i+1];
      $tmpseq.=substr($sequence{$contig},$lims[$i]-1,$lims[$i+1]-$lims[$i]+1);
    }
$pos.=")|";
if($sense eq "C"){
  $pos = "|c(".$pos;
  $tmpseq=~tr/("a","c","g","t","A","C","G","T")/("t","g","c","a","T","G","C","A")/;
        $seq=reverse($tmpseq);
}else{
  $pos="|(".$pos;
  $seq=$tmpseq;
}
###END OF CORRECTION
  $chain=">".$contig.$pos.$fam."|".$elem."\n".$seq."\n";
  return($chain);
}

### Setting up the research of fasta file
### by choosing on which file/directory it will takes place
if($fasta_file ne "no_value"){
  if($fasta_file eq ""){
    if($file_ok==0){
      $fasta_file=$rm_file;
      warn "Searching for fasta files in $fasta_file\n";
    }
  }
  if($fasta_file eq ""){
    if($file_ok==1){
      my @tmp=fileparse($rm_file);
      $fasta_file=$tmp[1];
      warn "Searching for fasta files in $fasta_file\n";
    }
  }
  find(\&Wanted_Fasta, $fasta_file);
  foreach $local_contig (keys %sequence){
    $contig_length{$local_contig}=length($sequence{$local_contig});
  }
}

### Subroutine handling the research of all available fasta sequences
sub Wanted_Fasta{
  if(/.*\.fa$/){
    warn "Looking for fasta sequences  in $File::Find::name\n";
    ### Checking existence of the fasta file if necessary
    open FASTA,"<$_" or die("Cannot open file $fasta_file; dying");
    ### Then  putting all contigs in memory
    while(<FASTA>){
      if(/^\>(.*)/){
  @getname=split /\s/, $1;
  $local_contig=$getname[0];
  chomp $local_contig;
  next;
      }
      chomp;
      $sequence{$local_contig}.=$_;
    }
  }
}



### Subroutine handling the construction of the reference length file, if necessary
### Here we only associate elements and lengths, and do not decide yet which length is the good one
%elements=();
%elem_length=();

sub Wanted_length{
  #Looking for .fa.out files in a directory, or taking directly the input RM file
  if(/.*\.out$/  | ($file_ok==1)){
    warn "Looking for element lengths  in $File::Find::name\n";
    open FILE, "<$_" or die("Cannot open file $File::Find::name\n");
    while(<FILE>){
      @F=split;
      ### Accounting RC Helitron as DNA transposons
      if($F[10] eq "RC/Helitron"){$F[10]="DNA/RC";}
      ### Filtering everything excepted transposable elements
      next unless(($F[10] =~ /^DNA/) | ($F[10] =~/^SINE/) | ($F[10] =~ /^LINE/) | ($F[10] =~ /^LTR/) | (($F[10] =~/^Unknown|^Unspecified/) & ($unknown == 1)));
      ### Filtering RM first lines
      next if(($F[0] eq "SW") | ($F[0] eq "score") | /^\s+$/);
      ### computing reference length for the line under study
      if($F[8] eq "+"){
  $remain=$F[13];
  $remain=~s/[\(\)]//g;
      }else{
  $remain=$F[11];
  $remain=~s/[\(\)]//g;
      }
      $lt=$F[12]+$remain;
      $elements{$F[9]}{$lt}++;
    }
  }
}

if($length_file eq ""){
  ### No file_length was given as input so we compute lengths
  find(\&Wanted_length, $rm_file);
  warn "\n###############################\n";
  warn "Computing elements lengths \n";

  ### Choosing the most common length as the good ones
  $length_file=$rm_file;
  if($length_file=~/(.*)\/$/){$length_file=$1};
  $length_file.=".length";
  unless($dry==1){
    open LENGTH_FILE, ">$length_file" or die("Cannot write/create file $length_file\n");
  }
  foreach $elem (sort keys %elements){
    @ll=sort {$elements{$elem}{$b} <=> $elements{$elem}{$a}} keys %{$elements{$elem}};
    unless($dry ==1){print LENGTH_FILE "$elem\t$ll[0]\n";}
    $elem_length{$elem}=$ll[0];
  }
}else{
  ### A file length was given, no computation done, reading only
  warn "\n###############################\n";
  warn "Reading elements lengths in file $length_file \n";
  ### We read the length file
  open LENGTH, "<$length_file" or die("Cannot open file $length_file");
  while(<LENGTH>){
      chomp;
      @F=split;
      $elem_length{$F[0]}=$F[1];
  }
}

### Summary of the elements that will be studied in this run
$elem_number = scalar(keys %elem_length);
warn "Element lengths found for $elem_number elements\n";
warn "###############################\n";



### Parsing the dictionary file to know about LTR/internal elements associations
%ltr=();
%internal=();

open LTR,"<$ltr_file" or die("Cannot open file $ltr_file");
warn "\n###############################\n";
warn "Reading LTR dictionary $ltr_file\n";

while(<LTR>){
  chomp;
  @F=split /\s+/,$_;
  if($#F>=1){
    @G=split /:/,$F[0];
    foreach $g (@G){
      $ltr{$g}=$F[1];
    }
    @H=split /:/,$F[1];
    foreach $h (@H){
      $internal{$h}=$F[0];
    }
  }else{
    ### All elements with no known structure are considered as internal -- with no consequence
    @G=split /:/,$F[0];
    foreach $g (@G){
      $ltr{$g}="NA";
    }
  }
}
$l=scalar(keys %ltr);

### Summary of the dictionary contents
warn "$l entries read\n";
warn "###############################\n";


### Here is the main code, which is one big find subroutine, in order to be recursive
find(\&Wanted, $rm_file);

sub Wanted{
  if(/.*\.out$/ | ($file_ok==1)){
    warn "Working on  $File::Find::name\n";
    open RM,"<$_" or die("Cannot open file $File::Find::name");
    $file=$_;
    warn "File treated is $file\n";

    ### Initialization
    %h=();
    %full_list=();
    %mem=();
    %chrom_size=();

    ### Setting up log file
    $outlog=$file.".log.txt";
    open OUTLOG,">$outlog" or die("Cannot open (create?) file $outlog");

    ### First we read once all RM files to put everything needed in memory
    while(<RM>){
      next if(/^\s+?SW/);
      next if(/^score/);
      next if(/^\s+$/);
      @F=split;
      ### Computing scaffold size
      if($chrom_size{$F[4]}==0){
  $F[7]=~s/[\(\)]//g;
  $chrom_size{$F[4]}=$F[6]+$F[7];
      }
      ### Filtering for other repeats and recording simple statistics about them
      if(($F[10]=~/^Low_complexity/) | ($F[10] =~/^Simple_repeat/) | ($F[10] =~/^Satellite/)){
  $useless_class = $F[10] ;
  $length = $F[6]-$F[5]+1;
  $non_TE_length{$F[4]}{$useless_class}+=$length;
  $non_TE_count{$F[4]}{$useless_class}++;
  next;
      }
      ### If no --unknown passed, then "Unknown" repeats are not considered, except for summary statistics
      if(($F[10] =~/^Unknown|^Unspecified/) & ($unknown != 1)){
  $useless_class = $F[10] ;
  $length = $F[6]-$F[5]+1;
  $non_TE_length{$F[4]}{$useless_class}+=$length;
  $non_TE_count{$F[4]}{$useless_class}++;
  next;
      }
      next unless(($F[10] =~ /^DNA/) | ($F[10] =~/^SINE/) | ($F[10] =~ /^LINE/) | ($F[10] =~ /^LTR/) | ($F[10] =~ /^RC/) | (($F[10] =~/^Unknown|^Unspecified/) & ($unknown == 1)));
      ### We do not treat elements not present in the length file.
      next unless(exists($elem_length{$F[9]}));
      ### Some modifications on the line to get useful info such as the real length of the element
      chomp;
      $_=~s/^\s+//g;
      @now=split /\s+/,$_;
      if($now[8] eq "+"){
  $lt=$now[12]-$now[11]+1;
      }else{
  $lt=$now[12]-$now[13]+1;
      }
      splice(@now,7,1,$lt);
      ### If there is something in the 15th column (typically a * symbol) we remove it
      if($#now >=15){
  splice(@now,15);
      }
      ### Adding one column which will contain the fragment number inside each copy
      $now[15]=1;
      ### Accounting RC Helitron as DNA transposons
      if($now[10] eq "RC/Helitron"){$now[10]="DNA/RC";}
      ### Accounting for Unknown LTR elements
      if((($now[10] eq "Unknown") | ($now[10] eq "Unspecified")) & (exists($ltr{$now[9]}) | exists($internal{$now[9]})) ){$now[10]="Unknown"}
      $tmp=join "\t",@now;
      ### %h is the main hash with everything recorded by element and family
      push @{$h{$F[4]}{$now[10]}{$now[9]}},$tmp;
      ### %mem is the hash associating each copy to its composition in fragments
      $mem{$tmp}=$tmp;
      ### %full_list is a cruder hash containing all RM lines with elements, to be displayed
      ### if choice is asked.
      if(($choice==1) & ($dry!=1)){push @{$full_list{$F[4]}}, $tmp;}
    }

    warn "RepeatMasker File contains ", scalar(keys %h)," contig(s)\n";
    print OUTLOG "RepeatMasker File contains ", scalar(keys %h)," contig(s)\n";

    ### Checking that the contig names in the fasta files
    ### correspond with those given in the RM files
    if($fasta_file ne "no_value"){
      @temp_cont=sort keys %h;
      @temp_cont_2=sort keys %sequence;
      for($i=0;$i<=$#temp_cont;$i++){
  $flag_cont=0;
  for($j=0;$j<=$#temp_cont_2;$j++){
    if($temp_cont[$i] eq $temp_cont_2[$j]){
      $flag_cont=1;
      last;
    }
  }
  if($flag_cont==0){
    ### If at least one contig in the RM file has no fasta counterpart
    ### disactivation of the fasta mode
    warn "FASTA MODE DISACTIVATED\n";
    warn "Contig $temp_cont[$i] not found in the provided fasta files\n";
    warn "Contigs found in RM file : @temp_cont\n";
    warn "Contigs found in Fasta file: @temp_cont_2\n";
    print OUTLOG "FASTA MODE DISACTIVATED\n";
    print OUTLOG "Contig $temp_cont[$i] not found in the provided fasta files\n";
    print OUTLOG "Contigs found in RM file : @temp_cont\n";
    print OUTLOG "Contigs found in Fasta file: @temp_cont_2\n";
    $fasta_file="no_value";
    $flanking=0;
    last;
       }else{
    for($i=0;$i<=$#temp_cont_2;$i++){
      warn "Fasta contig $temp_cont_2[$i] has length $contig_length{$temp_cont_2[$i]}\n";
      print OUTLOG "Fasta contig $temp_cont_2[$i] has length $contig_length{$temp_cont_2[$i]}\n";
    }
       }
      }
    }

    #### Writing the title in output files
    @tit=("Score","%_Div","%_Del","%_Ins","Query","Beg.","End.", "Length","Sense","Element","Family","Pos_Repeat_Beg","Pos_Repeat_End","Pos_Repeat_Left","ID","Num_Assembled", "%_of_Ref");
    $title=join "\t",@tit;
    $title.="\n";

    ### Fragment assembly:
    ### all fragments are assembled as copies
    ### As this is the only step for non-LTR transposons,
    ### they are detected and printed here.
    warn "\n###############################\n";
    warn "Finding non LTR elements\n";
    print OUTLOG "\n###############################\n";
    print OUTLOG "Finding non LTR elements\n";

    if(($choice==1) & ($dry!=1)){warn "\nTo exit choice mode type q when asked for a solution\n";}
    warn "###############################\n\n";

    if(($choice==1) & ($dry!=1)){print OUTLOG "\nTo exit choice mode type q when asked for a solution\n";}
    print OUTLOG "###############################\n\n";

    ### Then we work on each scaffold independently
    foreach $contig (sort keys %h){
      warn "Starting to work on contig $contig\n";
      print OUTLOG "Starting to work on contig $contig\n";

      ### Initialization
      $counterchoice=0;
      $counternonltr=0;
      %mem_stat=();
      $mem_solo_ltr=();
      %family_list=();
      %mem_number=();
      @full_list_by_position=();
      @positions_list=();
      %all_copies=();

      ### Making list of all contig elements by position
      @full_list_by_position=sort by_position @{$full_list{$contig}};
      foreach $p (@full_list_by_position){
  @F=split /\t/,$p;
  push @positions_list,$F[5];
      }

      ### Here  we define all scaffold dependent output files were resultats will be written
      ### Nothing is created except the log file in --dry-run mode
      unless($dry==1){
  $out=$file."_".$contig.".transposons.csv";
  open OUT,">$out" or die("Cannot open (create?) file $out");

  $outltr=$file."_".$contig.".ltr.csv";
  open OUTLTR,">$outltr" or die("Cannot open (create?) file $outltr");

  $outorder=$file."_".$contig.".elem_sorted.csv";
  open OUTORDER,">$outorder" or die("Cannot open (create?) file $outorder");

  $outstat=$file."_".$contig.".copynumber.csv";
  open OUTSTAT,">$outstat" or die("Cannot open (create?) file $outstat");

  if($fasta_file ne "no_value"){
    $outfasta=$file."_".$contig;
    if($flanking>0){$outfasta.="_flank_".$flanking;}
    $outfasta.=".fasta";
    open OUTFASTA,">$outfasta" or die("Cannot open (create?) file $outfasta");
  }

  print OUT $title;
  print OUTLTR $title;
  print OUTORDER $title;
  print OUTSTAT "Family\tElement\tLength\tFragments\tCopies\tSolo_LTR\tTotal_Bp\tCover\n";
      }

      ### Then we work on each class and family.
      foreach $family (sort keys %{$h{$contig}}){
  foreach $elem (sort keys %{$h{$contig}{$family}}){
    $counterelem=0;
    $totbase_elem=0;

    ### @trans is the array containing all relevant RM lines for the element under study
    @trans=@{$h{$contig}{$family}{$elem}};
    $num=$#trans;
    $t=$trans[0];

    ### We compute the length $lt of the element in the database
    ### This is done line-wise to prevent some trouble when the reference length
    ### does not match the reference computed on the line
    @now=split /\s+/,$t;
    if($now[8] eq "+"){
      $remain=$now[13];
      $remain=~s/[\(\)]//g;
      $lt=$now[12]+$remain;
    }else{
      $remain=$now[11];
      $remain=~s/[\(\)]//g;
      $lt=$now[12]+$remain;
    }

    ### Sanity check, to be sure that the counts of LTR and internal subparts are well cumulated
    $numbis=$num+1;
    $mem_number{$family}{$elem}+=$numbis;
    if(exists($internal{$elem})){$mem_number{$family}{$internal{$elem}}+=$numbis;}
    if(exists($ltr{$elem})){$mem_number{$family}{$ltr{$elem}}+=$numbis;}

    unless($dry==1){
      warn "Assembling element $elem of family $family; length $lt; there are $numbis fragments\n";
    }
    print OUTLOG "Assembling element $elem of family $family; length $lt; there are $numbis fragments\n";
    @newtrans=();

    ### Taking in turn each fragment to be the beginning of a combination of fragments
    while($#trans>=0){
      if($choice==1){
        @combined=();
        @posflag=();
        %tmpmem=();
      }
      $t=shift @trans;
      @now=split /\s+/,$t;
      ### $flag will be 0 if we find no other part of this element
      ### and will pass to 1 if we find one.
      $flag=0;
      ### Checking if two lines of the same element
      ### Could be a single element
      ### Even if they do not follow each other
      ### Conditions: sense should be the same, positions should be close enough
      ### All sub-parts should be ordered in the same way in the element found and the real element (ie no recombination)
      $j=-1;
      ### Here we parse all other fragments for the same element, and check each one in turn
      ### for an assembly with the element considered
      foreach $tbis (@trans){
        $j++;
        ### Could $t and $tbis be two parts of the same element?
        @second=split /\s+/,$tbis;
        if($insert>-1){
    last unless($second[5]-$now[6]-1 <= $insert);
        }else{
    last unless($second[6]-$now[5] < 2*$lt);
        }
        next unless($second[8] eq $now[8]);

        ### Here we check that the two elements we want to assemble are not
        ### separated by an internal part if they are LTR
        ### or by a LTR part if they are internal
        $man_in_the_middle=0;
        if(exists($ltr{$elem})){
    ### Then $elem is an internal elements
    ### We take every possible LTR corresponding element and test for its presence
    @possible_ltr=split /:/,$ltr{$elem};
    foreach $poss (@possible_ltr){
      @middle=@{$h{$contig}{$family}{$poss}};
      foreach $m (@middle){
        @mid=split /\s+/,$m;
        if(($mid[5] > $now[5]) & ($mid[5]<$second[5]) & ($mid[8] eq $now[8])){
          $man_in_the_middle=1;
          $middleguy=$m;
          last;
        }
      }
    }
        }
        if(exists($internal{$elem})){
    ### Then $elem is an LTR elements
    ### Then we take every possible internal corresponding element and test for its presence
    @possible_int=split /:/,$internal{$elem};
    foreach $poss (@possible_int){
      @middle=@{$h{$contig}{$family}{$poss}};
      foreach $m (@middle){
        @mid=split /\s+/,$m;
        if(($mid[5] > $now[5]) & ($mid[5]<$second[5]) & ($mid[8] eq $now[8])){
          $man_in_the_middle=1;
          $middleguy-$m;
          last;
        }
      }
    }
        }
        ### If there is a "man in the middle", we do not go further
        ### LTR and internal part assembly will be made later on
        last if($man_in_the_middle==1);

        ### Check if the genomic positions follow each other at least partially (useful only if long-range choices have been previously made)
        next unless(($second[5]>$now[5]) & ($second[6]>$now[6]));

        ### Check of the positions on the reference element follow each other at least partially
        if($now[8] eq "+"){
    $end_delta=$second[12]-$now[12];
    $ori_delta=$second[11]-$now[11];
        }else{
    $ori_delta=$now[12]-$second[12];
    $end_delta=$now[13]-$second[13];
        }
        if(($ori_delta>0) & ($end_delta>0)){
    ### We found another line combining well with the one under study
    ### We will replace the first transposon in the table
    ### by the combination of both lines, in order to find for
    ### fragmented transposons
    $flag++;
    ### Here we do different things depending on the fact that we
    ### want to chose or not
    if($choice!=1){
      ### Erasing the secondary element
      splice(@trans,$j,1);
      ### Computing characteristics of the combined element
      $tmpcombined=&combine_lines(\@now,\@second);
      ### Inserting combined element in the array
      unshift @trans,$tmpcombined;
      ### Memorizing the lines used to compute this element
      $mem{$tmpcombined}=$mem{$t}."\n".$tbis;
      last;
      ### If we have found a combination, we do not check for other independent ones
      ### but research will restart starting with the "combined" element
    }else{
      ## Where we memoriz, in order to choose
      $posflag[$flag-1]=$j;
      $combined[$flag-1]=&combine_lines(\@now,\@second);
      $tmpmem{$combined[$flag-1]}=$mem{$t}."\n".$tbis;
    }
        }
      }
      if($flag==0){
        ### Nothing combines well with the line under study
        unless(exists($ltr{$elem}) | exists($internal{$elem})){
    ### Then, unless it is a LTR subpart,  print it
    ### without strict all the time,
    ### if strict if has to follow the 80-80 rule
    ### 80 pb minimum  and 80% identity to real elements
    if(($strict eq "") |  (($now[7]>=80) & ($now[1])<20) ){
      ### Adjusting global counters
      $counternonltr++;
      $counterelem++;
      $totbase_elem+=$now[7];
      ### Computing length relative to the reference, if any
      if(exists($elem_length{$now[9]})){
        $relative_length=$now[7]/$elem_length{$now[9]};
        $relative_length_print=sprintf("%.3f",$relative_length);
      }else{
        $relative_length_print="No_ref_available";
      }
      unless($dry==1){
        print OUT "\n###$t\t",$relative_length_print,"\n";
        ### %all_copies is an hash containing all assembled fragments,
        ### that will be used to sort them and print them in genome order
        ### at the end of the run
        $all_copies{$t}="\n###$t\t$relative_length_print\n";
        ### If the element was a combination of fragments,
        ### print the fragments under it
        if(($mem{$t}) ne $t){
          print  OUT "$mem{$t}\n";
          $all_copies{$t}.="$mem{$t}\n";
        }
        if($fasta_file ne "no_value"){
          print OUTFASTA &get_fasta_seq($t,$mem{$t});
        }
      }
    }
        }else{
    ### LTR element : we just put it back in memory to be assembled with other in next part.
    $counterelem++;
    push @newtrans, $t;
        }
        ### Going directly to another start for fragment assembly.
        next;
      }
      if($choice == 1){
        ### Then we want to choose, and there may be ambiguous cases
        if($flag==1){
    ### Only one combination was found
    ### it is accepted without user choice
    ### Then, exactly as if it was accepted by user
    ### put it back into the list
    splice(@trans,$posflag[0],1);
    unshift @trans,$combined[0];
    $mem{$combined[0]}=$tmpmem{$combined[0]};
        }
        if($flag>1){
    ### Many solutions found
    ### First we check if RM Block ID can guide the user in assembling subparts
    $check_rm_block_id=0;
    $sol="";
    for($k=0;$k<$flag;$k++){
      @F=split /\t/,$combined[$k];
      @G=split /\//,$F[14];
      for($l=0;$l<$#G;$l++){
        if($G[$l] == $G[$#G]){
          ### RM Block IDs match : the choice is bypassed and the assembly automatic
          unless($dry==1){warn "\nChoice bypassed -- Assembling fragments based on RM Block ID\n";}
          unless($dry==1){warn "###$combined[$k]\n$tmpmem{$combined[$k]}\n";}
          print OUTLOG "\nChoice bypassed -- Assembling fragments based on RM Block ID\n";
          print OUTLOG "###$combined[$k]\n$tmpmem{$combined[$k]}\n";
          $check_rm_block_id=1;
          $sol=$k;
          last;
        }
      }
      last if($check_rm_block_id==1);
    }
    ### End check RM Block ID
    ### %sol_list is the hash that will contains the multiple solutions from which to choose
    %sol_list=();
    if($check_rm_block_id==0){
      ### RM Block Id did not match, so we choose manually

      ### Size of the context to be printed
      $maxcountlines=11;

      ### Printing the possible solutions with fragment assembly
      unless($dry==1){warn "\nCannot decide between \n";}
      print OUTLOG "\nCannot decide between \n";
      for($k=0;$k<$flag;$k++){
        $sol_list{$k}=1;
        unless($dry==1){warn "\nSolution $k\n###$combined[$k]\n$tmpmem{$combined[$k]}\n";}
        print OUTLOG "\nSolution $k\n###$combined[$k]\n$tmpmem{$combined[$k]}\n";
      }
      ###Printing solution "alone"
      $k=$flag;
      $sol_list{$k}=1;
      $toprint=join("\t",@now);
      unless($dry==1){warn "\nSolution $k\n###$toprint\n";}
      print OUTLOG "\nSolution $k\n###$toprint\n";

      ### Printing the context to help choice
      unless($dry==1){
        warn "\n";
        warn "Context (type 'm' to see a broader context if necessary):\n";
        warn "$t\n";
      }
      print OUTLOG "\n";
      print OUTLOG "Context (type 'm' to see a broader context if necessary):\n";
      print OUTLOG "$t\n";
      @G=split /\t/,$trans[0];
      for($k=0;$k<=5;$k++){
        ### Check context to not display too much
        @F=split /\t/,$trans[$k];
        if(($F[5]-$G[5])>100000){
          last;
        }
        unless($dry==1){warn "$trans[$k]\n";}
        print OUTLOG "$trans[$k]\n";
      }
      unless($dry==1){
        warn "\n";
        warn "Which solution do you prefer?\n";
      }
      print OUTLOG "\n";
      print OUTLOG "Which solution do you prefer?\n";

      ### Reading user input and going on accordingly
      while(! exists($sol_list{$sol})){
        unless($dry==1){warn "Please type corresponding number\n";}
        print OUTLOG "Please type corresponding number\n";
        if($dry==1){$sol=0;$counterchoice++;$global_choice_counter++;}else{$sol=<>;}
        print OUTLOG $sol,"\n";
        if($sol eq "q\n"){
          ### Setting 0 for the current choice and exiting choice mode
          $sol=0;
          $choice=0;
          warn "User typed q; all further choices will be automatically 0\n";
          print OUTLOG "User typed q; all further choices will be automatically 0\n";
        }
        if($sol eq "m\n"){
          ### Printing new, broader context
          warn "User typed m; showing enlarged context\n";
          print OUTLOG "User typed m; showing enlarged context\n";
          $countlines=0;
          for($localpos=0;$localpos<=$#positions_list;$localpos++){
      next if($positions_list[$localpos]<$now[5]);
      $countlines++;
      if(($countlines>=$maxcountlines) | ($positions_list[$localpos]>$now[5]+100000)){
        ### Increasing maximum context size  to anwer to another m command
        $maxcountlines+=5;
        last;
      }
      warn $full_list_by_position[$localpos],"\n";
      print OUTLOG $full_list_by_position[$localpos],"\n";
          }
        }
        chomp $sol;
      }
      unless($dry==1){
        warn "Choice recorded\n\n";
      }
      print OUTLOG "Choice recorded\n\n";
    }

    ### If a numeric solution was given, then select if and go on for fragment assembly
    if($sol<=$flag){
      splice(@trans,$posflag[$sol],1);
      unshift @trans,$combined[$sol];
      $mem{$combined[$sol]}=$tmpmem{$combined[$sol]};
    }else{
      ### The "alone" solution  was chosen:
      ### the fragment is supposed to be a copy by itself
      ### and has to be printed.
      unless(exists($ltr{$elem}) | exists($internal{$elem})){
        ### Then, unless it is a LTR subpart, print it
        ### without strict all the time,
        ### if strict if has to follow the 80-80 rule
        ### 80 pb minimum  and 80% identity to real elements
        if(($strict eq "") |  (($now[7]>=80) & ($now[1])<20) ){
          ### Adjusting global counters
          $counternonltr++;
          $counterelem++;
          $totbase_elem+=$now[7];
          ### Computing length relative to the reference, if any
          if(exists($elem_length{$now[9]})){
      $relative_length=$now[7]/$elem_length{$now[9]};
      $relative_length_print=sprintf("%.3f",$relative_length);
          }else{
      $relative_length_print="No_ref_available";
          }
          unless($dry==1){
      print OUT "\n###$t\t",$relative_length_print,"\n";
      $all_copies{$t}="\n###$t\t$relative_length_print\n";
      if(($mem{$t}) ne $t){
        print  OUT "$mem{$t}\n";
        $all_copies{$t}.="$mem{$t}\n";
      }
      if($fasta_file ne "no_value"){
        print OUTFASTA &get_fasta_seq($t,$mem{$t});
      }
          }
        }
      }else{
        ### LTR element : we just put it bak in memory to be assembled with other in next part.
        $counterelem++;
        push @newtrans, $t;
      }
    }
        }
      }
    }
    ### Printing some summary of the full procedure for the element
    unless($dry==1){
      warn "$numbis fragments were assembled into $counterelem copie(s)\n";
    }
    print OUTLOG "$numbis fragments were assembled into $counterelem copie(s)\n";
    unless(exists($ltr{$elem}) | exists($internal{$elem})){
     ### Memorizing elements data for later computations.
      $percent=sprintf("%.2f",100*$totbase_elem/$chrom_size{$contig});
      $mem_stat{$family}{$elem}=[$lt,$numbis,$counterelem,$totbase_elem,$percent]
    }else{
      $h{$contig}{$family}{$elem}=[@newtrans];
    }
  }
      }

      warn "\n###############################\n";
      warn "$counternonltr copies of non-LTR elements were found\n";
      warn "###############################\n";
      print OUTLOG "\n###############################\n";
      print OUTLOG "$counternonltr copies of non-LTR elements were found\n";
      print OUTLOG "###############################\n";


      ### Finally we go through all assembled LTR elements
      ### And check if we find a structure LTR-I-LTR
      warn "\n###############################\n";
      warn "Finding complete LTR elements\n";
      warn "Results written in file $outltr\n";
      warn "###############################\n\n";
      print OUTLOG "\n###############################\n";
      print OUTLOG "Finding complete LTR elements\n";
      print OUTLOG "Results written in file $outltr\n";
      print OUTLOG "###############################\n\n";

      ### Initialization
      $counterltr=0;
      $counterltrpartial=0;
      $solo_ltr=0;
      %solo_ltr_elem=();

      ### We take every family in turn
      foreach $family (sort keys %{$h{$contig}}){
  foreach $elem (sort keys %{$h{$contig}{$family}}){
    $totbase_elem=0;
    $counterelem=0;
    ### Selection of a family of LTR and I elements
    next unless ((exists $ltr{$elem}) || ((exists($internal{$elem})) && (! exists($h{$contig}{$family}{$internal{$elem}}))));
    ### Here we make a list of all elements that could be assembled
    ### taking both LTR and internal subparts
    @possible_ltr = split /:/,$ltr{$elem};
    @tableau=@{$h{$contig}{$family}{$elem}};
    foreach $poss (@possible_ltr){
      @tableau=(@tableau,@{$h{$contig}{$family}{$poss}});
    }
    @tableau_sorted=sort by_position @tableau;
    $elem_family_size=scalar(@tableau_sorted);
    unless($dry==1){
      warn "Building LTR-ET with internal part $elem of family $family -- there are $elem_family_size such elements\n";
    }
    print OUTLOG "Building LTR-ET with internal part $elem of family $family -- there are $elem_family_size such elements\n";

    ### Now we go through tableau_sorted, looking for the pattern LTR-I-LTR
    ### We allow half the LTR length after a LTR before the I
    ### If not we dismiss the LTR and go directly to the next one
    @tableau_sorted_copy=@tableau_sorted;
    ### $pos is the current position in the array of elements
    $pos=-1;
    while($#tableau_sorted>=0){
      $pos++;
      $t=shift @tableau_sorted;
      @now=split /\s+/,$t;
      ### A priori we only choose LTR copies to start with
      next unless(exists($internal{$now[9]}));
      ### Computing reference length
      if($now[8] eq "+"){
        $remain=$now[13];
        $remain=~s/[\(\)]//g;
        $lt=$now[12]+$remain;
      }else{
        $remain=$now[11];
        $remain=~s/[\(\)]//g;
        $lt=$now[12]+$remain;
      }
      $tbis=$tableau_sorted[0];
      @second=split /\s+/,$tbis;
      ### Conditions on the internal part
      if((exists($ltr{$second[9]})) & ($second[5] - $now[6] < 0.5*$lt) & ($now[8] eq $second[8])){
        $tter= $tableau_sorted[1];
        @third=split /\s+/,$tter;
        ### Conditions on the final LTR
        if((exists($internal{$third[9]})) & ($third[5] - $second[6] < 0.5*$lt) & ($third[8] eq $second[8])){
    ### We found a complete element and report it
    $counterltr++;
    $tmpcombined=&build_ltr(\@now,\@second);
    @comb=split /\s+/,$tmpcombined;
    $tmpcombined=&build_ltr(\@comb,\@third);
    @checklength = split /\t/,$tmpcombined;
    ### Computation of the element relative size
    $total_elem_length=$elem_length{$second[9]}+$elem_length{$now[9]}+$elem_length{$third[9]};
    $relative_length_print=sprintf("%.3f",$checklength[7]/$total_elem_length);
    unless($dry==1){
      print OUTLTR "\n###$tmpcombined\t",$relative_length_print,"\n";
      $all_copies{$tmpcombined}="\n###$tmpcombined\t$relative_length_print\n";
      ### For each subpart (LTR, I, or LTR2),
      ### we check if it was built from fragments, and if so
      ### we print all fragments
      if(($mem{$t}) ne $t){
        print  OUTLTR "$mem{$t}\n";
        $all_copies{$tmpcombined}.="$mem{$t}\n";
      }else{
        print OUTLTR "$t\n";
        $all_copies{$tmpcombined}.="$t\n";
      }
      if(($mem{$tbis}) ne $tbis){
        print  OUTLTR "$mem{$tbis}\n";
        $all_copies{$tmpcombined}.="$mem{$tbis}\n";
      }else{
        print OUTLTR "$tbis\n";
        $all_copies{$tmpcombined}.="$tbis\n";
      }
      if(($mem{$tter}) ne $tter){
        print  OUTLTR "$mem{$tter}\n";
        $all_copies{$tmpcombined}.="$mem{$tter}\n";
      }else{
        print OUTLTR "$tter\n";
        $all_copies{$tmpcombined}.="$tter\n";
      }
      if($fasta_file ne "no_value"){
        print OUTFASTA &get_fasta_seq($t,$mem{$t}."\n".$mem{$tbis}."\n".$mem{$tter});
      }
    }
    ### Counter adjustment
    $totbase_elem+=$checklength[7];
    $counterelem++;
    ### And we erase the two following elements in the table,
    ### as they have been used to make a complete LTR retrotransposon
    splice(@tableau_sorted_copy,$pos,3);
    splice(@tableau_sorted,0,2);
    $pos=$pos-1;
        }
      }
    }
    ### Now we go looking for partial LTR elements
    ### We took away all copies, and remaining parts in tableau_sorted must be partial
    ### We just say that if an I element and a LTR element are closer than half of the LTR length, we are good
    while($#tableau_sorted_copy>=0){
      if($#tableau_sorted_copy==0){
        ### There is only one remaining element: we print it.
        $t=shift @tableau_sorted_copy;
        @now=split /\s+/,$t;
        if(($strict eq "") |  (($now[7]>=80) & ($now[1])<20) ){
          if(exists($ltr{$now[9]})){
      $counterltrpartial++;
      @ref_ltr_elem_local=split /:/,$ltr{$now[9]};
      if(exists($elem_length{$now[9]}) & exists($elem_length{$ref_ltr_elem_local[0]})){
        $total_elem_length=$elem_length{$now[9]}+2*$elem_length{$ref_ltr_elem_local[0]};
      }else{$total_elem_length="No_ref_available"}
          }else{
      ### Solo LTR special accounting and reference length
      if(exists($elem_length{$now[9]})){
        $total_elem_length=$elem_length{$now[9]};
      }else{
        $total_elem_length="No_ref_available";
      }
      $solo_ltr++;
      $solo_ltr_elem{$family}{$elem}++;
          }
          $totbase_elem+=$now[7];
          $counterelem++;
          if($total_elem_length>0){
      $relative_length_print=sprintf("%.3f",$now[7]/$total_elem_length);
          }else{
      $relative_length_print="No_ref_available";
          }
          unless($dry==1){
      print OUTLTR "\n###$t\t",$relative_length_print,"\n";
      $all_copies{$t}="\n###$t\t$relative_length_print\n";
      if(($mem{$t}) ne $t){
        print  OUTLTR "$mem{$t}\n";
        $all_copies{$t}.="$mem{$t}\n";
      }
      if($fasta_file ne "no_value"){
        print OUTFASTA &get_fasta_seq($t,$mem{$t});
      }
          }
        }
        last;
      }
      $t=shift @tableau_sorted_copy;
      @now=split /\s+/,$t;
      $tbis=$tableau_sorted_copy[0];
      @second=split /\s+/,$tbis;
      ### Check if the 2 first elements combine well
      if(($second[9] ne $now[9]) & ($second[5] - $now[6] < 0.5*$lt) & ($now[8] eq $second[8])){
        $tmpcombined=&build_ltr(\@now,\@second);
        @checklength=split /\t/,$tmpcombined;
        if(($strict eq "") |  (($checklength[7]>=80) & ($checklength[1])<20) ){
          if(exists($elem_length{$now[9]}) & exists($elem_length{$second[9]})){
      if(exists($ltr{$now[9]})){
        $total_elem_length=$elem_length{$now[9]}+2*$elem_length{$second[9]};
      }else{
        $total_elem_length=2*$elem_length{$now[9]}+$elem_length{$second[9]};
      }
          }else{$total_elem_length="No_ref_available"}
          $counterltrpartial++;
          $totbase_elem+=$checklength[7];
          $counterelem++;
          if($total_elem_length>0){
      $relative_length_print=sprintf("%.3f",$checklength[7]/$total_elem_length);
          }else{
      $relative_length_print="No_ref_available";
          }
          unless($dry==1){
      print OUTLTR "\n###$tmpcombined\t",$relative_length_print,"\n";
      $all_copies{$tmpcombined}="\n###$tmpcombined\t$relative_length_print\n";
      if(($mem{$t}) ne $t){
        print  OUTLTR "$mem{$t}\n";
        $all_copies{$tmpcombined}.="$mem{$t}\n";
      }else{
        print OUTLTR "$t\n";
        $all_copies{$tmpcombined}.="$t\n";
      };
      if(($mem{$tbis}) ne $tbis){
        print  OUTLTR "$mem{$tbis}\n";
        $all_copies{$tmpcombined}.="$mem{$tbis}\n";
      }else{
        print OUTLTR "$tbis\n";
        $all_copies{$tmpcombined}.="$tbis\n";
      };
      if($fasta_file ne "no_value"){
        print OUTFASTA &get_fasta_seq($t,$mem{$t}."\n".$mem{$tbis});
      }
          }
          splice(@tableau_sorted_copy,0,1);
        }
      }else{
        ### We return the first element, which is isolated, and go to the next one
        if(($strict eq "") |  (($now[7]>=80) & ($now[1])<20) ){
          if(exists($ltr{$now[9]})){
      $counterltrpartial++;
      @ref_ltr_elem_local=split /:/,$ltr{$now[9]};
      if(exists($elem_length{$now[9]}) & exists($elem_length{$ref_ltr_elem_local[0]})){
        $total_elem_length=$elem_length{$now[9]}+2*$elem_length{$ref_ltr_elem_local[0]};
      }else{$total_elem_length="No_ref_available"}
          }else{
      ### Solo LTR special accounting
      $solo_ltr++;
      $solo_ltr_elem{$family}{$elem}++;
      if(exists($elem_length{$now[9]})){
        $total_elem_length=$elem_length{$now[9]};
      }else{
        $total_elem_length="No_ref_available";
      }
          }
          $totbase_elem+=$now[7];
          $counterelem++;
          if($total_elem_length>0){
      $relative_length_print=sprintf("%.3f",$now[7]/$total_elem_length);
          }else{
      $relative_length_print="No_ref_available";
          }
          unless($dry==1){
      print OUTLTR "\n###$t\t",$relative_length_print,"\n";
      $all_copies{$t}="\n###$t\t$relative_length_print\n";
      if(($mem{$t}) ne $t){
        print  OUTLTR "$mem{$t}\n";
        $all_copies{$t}.="$mem{$t}\n";
      };
      if($fasta_file ne "no_value"){
        print OUTFASTA &get_fasta_seq($t,$mem{$t});
      }
          }
        }
      }
    }
    ### Memorizing elements data for later computations.
    $percent=sprintf("%.2f",100*$totbase_elem/$chrom_size{$contig});
    $mem_stat{$family}{$elem}=[$total_elem_length,$mem_number{$family}{$elem},$counterelem,$totbase_elem,$percent];
    if(!exists $solo_ltr_elem{$family}{$elem}){$solo_ltr_elem{$family}{$elem}=0};
        }
      }
      warn "\n###############################\n";
      warn "$counterltr complete LTR elements found\n";
      warn "$counterltrpartial partial LTR elements found\n";
      warn "$solo_ltr solo LTR elements found\n";
      warn "###############################\n";

      print OUTLOG "\n###############################\n";
      print OUTLOG "$counterltr complete LTR elements found\n";
      print OUTLOG "$counterltrpartial partial LTR elements found\n";
      print OUTLOG "$solo_ltr solo LTR elements found\n";
      print OUTLOG "###############################\n";

      ### If this was a dry run, then just report how many choices have to be done on this contig
      if($dry==1){
        warn "\nIn this dry run $counterchoice choices have been made in the file $File::Find::name\n";
        print OUTLOG "\nIn this dry run $counterchoice choices have been made in the file $File::Find::name\n";
        warn("Dry run finished with no apparent error\n\n\n");
      }

      ### Quantitative display
      @families=("Type:DNA","Type:LINE","Type:SINE","Type:LTR");
      if($unknown==1){
        push @families,"Type:Unknown";
      }
      unless($dry==1){
        foreach $bigfam (@families){
    @{$stat_summary{$bigfam}}=();
        }
        @{$stat_summary{"EVERYTHING"}}=();

        $stat_solo_ltr{"Type:LTR"}=0;
        foreach $fam (sort keys %mem_stat){

    if($fam=~/^DNA/){
      $bigfam="Type:DNA";
    }elsif($fam=~/^LINE/){
      $bigfam="Type:LINE";
    }elsif($fam=~/^SINE/){
      $bigfam="Type:SINE";
    }elsif($fam=~/^LTR/){
      $bigfam="Type:LTR";
      $stat_solo_ltr{$fam}=0;
    }elsif($fam=~/^Unknown|^Unspecified/){
      $bigfam="Type:Unknown";
    }
    @{$stat_summary{$fam}}=();
    push @{$family_list{$bigfam}},$fam;
    ### Making sums inside each category (class, but also DNA, LINE, SINE...)
    foreach $elem (sort keys %{$mem_stat{$fam}}){
      for($i=1;$i<=3;$i++){
        $stat_summary{$bigfam}[$i-1]+=$mem_stat{$fam}{$elem}[$i];
        $stat_summary{$fam}[$i-1]+=$mem_stat{$fam}{$elem}[$i];
        $stat_summary{"EVERYTHING"}[$i-1]+=$mem_stat{$fam}{$elem}[$i];
      }
      if($fam =~ /^LTR/){
        $stat_solo_ltr{$bigfam}+=$solo_ltr_elem{$fam}{$elem};
        $stat_solo_ltr{$fam}+=$solo_ltr_elem{$fam}{$elem};
      }
    }
        }

        ### Printing all quantiative information
        foreach $bigfam (@families){
    foreach $fam (@{$family_list{$bigfam}}){
      foreach $elem (sort keys %{$mem_stat{$fam}}){
        print OUTSTAT $fam,"\t",$elem;
        for($i=0;$i<=2;$i++){
          print OUTSTAT "\t${$mem_stat{$fam}{$elem}}[$i]";
        }
        if($bigfam eq "Type:LTR"){
          print OUTSTAT "\t$solo_ltr_elem{$fam}{$elem}";
        }else{
          print OUTSTAT "\tNA";
        }
        for($i=3;$i<=4;$i++){
          print OUTSTAT "\t${$mem_stat{$fam}{$elem}}[$i]";
        }
        print OUTSTAT "\n";
      }
      print OUTSTAT "###",$fam,"\tAll_elements\tNA";
      for($i=0;$i<=1;$i++){
        print OUTSTAT "\t",$stat_summary{$fam}[$i];
      }
      if($bigfam eq "Type:LTR"){
          print OUTSTAT "\t$stat_solo_ltr{$fam}";
      }else{
          print OUTSTAT "\tNA";
      }
      print OUTSTAT "\t",$stat_summary{$fam}[2],"\t";
      $percent=sprintf("%.4f",100*$stat_summary{$fam}[2]/$chrom_size{$contig});
      print OUTSTAT "$percent\n\n";
    }
    print OUTSTAT "######",$bigfam,"\tAll_elements\tNA";
    for($i=0;$i<=1;$i++){
      print OUTSTAT "\t",$stat_summary{$bigfam}[$i];
    }
    if($bigfam eq "Type:LTR"){
          print OUTSTAT "\t$stat_solo_ltr{$bigfam}";
    }else{
          print OUTSTAT "\tNA";
    }
    print OUTSTAT "\t",$stat_summary{$bigfam}[2],"\t";
    $percent=sprintf("%.4f",100*$stat_summary{$bigfam}[2]/$chrom_size{$contig});
    print OUTSTAT "$percent\n\n";
        }

        print OUTSTAT "#########Type:EVERYTHING_TE\tAll_elements\tNA";
        for($i=0;$i<=1;$i++){
    print OUTSTAT "\t",$stat_summary{"EVERYTHING"}[$i];
        }
        print OUTSTAT "\tNA\t",$stat_summary{"EVERYTHING"}[2],"\t";
        $percent=sprintf("%.4f",100*$stat_summary{"EVERYTHING"}[2]/$chrom_size{$contig});
        print OUTSTAT "$percent\n\n";

        #### Printing quantitative info for non TE elements
        foreach $nonTE (sort keys %{$non_TE_length{$contig}}){
    print OUTSTAT "######Type:",$nonTE,"\tAll_elements\tNA\tNA\t$non_TE_count{$contig}{$nonTE}\tNA\t$non_TE_length{$contig}{$nonTE}\t";
    $percent=sprintf("%.4f",100*$non_TE_length{$contig}{$nonTE}/$chrom_size{$contig});
    print OUTSTAT "$percent\n\n";
        }
      }

      ### Printing the sorted file with all copies in genomic order
      foreach $k (sort by_position keys %all_copies){
        print OUTORDER $all_copies{$k};
      }

    }### end of contig under study
  }### end of file under study
}### end of Wanted Subroutine

### If this was a dry run, then just report how many choices have to be done globally
if($dry==1){warn "\nIn this dry run $global_choice_counter choices have been made in total on all files\n";}
