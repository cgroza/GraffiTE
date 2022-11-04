#!/usr/bin/perl

#     One code to find them all  -- perl utility ot extract information from RepeatMasker output files
#     Copyright (C) 2014 Bailly-Bechet Marc
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


### This script takes either an output from Repeat Masker
### or a directory containing them (the script is recursive).
### It reports all combination of internal and LTR elements found
### Its guesses about the elements are based on syntax (i.e. element names) only


use FileHandle;
use Getopt::Long;
use File::Find;

$fuzzy=0;
$unknown=0;
$help="";

GetOptions('rm=s' => \$rm_file,
        'unknown'=>\$unknown,
        'fuzzy' => \$fuzzy,
        'help' =>\$help);


if($rm_file eq ""){die("Dying! Option --rm must be passed\n")}


if($help==1){
  warn "#################\n";
  warn "Usage:\nbuild_dictionary.pl --rm file [--fuzzy] [--unknown] > output\n";
  warn "Options --rm requires a file or directory to be passed in order to work correctly;\n";
  warn "#################\n";
  warn "See README file for more details concerning the different options usage\n";
  exit(0);
}


warn "#################################\n\n";

# Initialization
%ltr_elem=();
%internal_elem=();
%matches=();

### Checking if we work on a directory or a single file either
unless(-d $rm_file){$file_ok=1}

find(\&Wanted, $rm_file);

sub Wanted{
  if(/.*\.out$/  | ($file_ok==1)){
    warn "Searching in $File::Find::name\n";
    open FILE, "<$_" or die("Cannot open file $File::Find::name\n");
    while(<FILE>){
      @F=split;
      next unless(($F[10] =~ /LTR/) | (($F[10] =~/^Unknown|^Unspecified/) & ($unknown == 1)));
      ### Memorizing the element;
      if($F[9]=~/[-_](in|int|i)/i){
  $internal_elem{$F[9]}=1;
  next;
      }
      if($F[9]=~/LTR/){### If element name contains LTR, it is supposed to be an LTR part
  $ltr_elem{$F[9]}=1;
      }else{### Else we stay undecided
  $ltr_elem{$F[9]}=1;
  $internal_elem{$F[9]}=1;
      }
    }
  }
}

warn "\n#################################\n";
warn "Finding matching elements\n";
warn "#################################\n";

@INT= sort keys %internal_elem;
@LTR= sort keys %ltr_elem;

%inside=();
%outside=();
@uniques=();
@uniquebis=();
%found=();

warn "\n#################################\n";
warn "Phase 1 : exact matches\n";
warn "#################################\n";

$matchcount=0;

### Elements containing "LTR" should have an intenral part with a corresponding I,IN, INT or int
foreach $ltr (@LTR){
  foreach $sub ("I","IN", "INT","int"){
    $tmp=$ltr;
    $tmp=~s/LTR/$sub/;
    if(($internal_elem{$tmp}==1) & ($ltr ne $tmp)){
      $inside{$tmp}=$ltr;
      $outside{$ltr}=$tmp;
      $matchcount++;
      next;
    }
  }
}

### Elements not containing LTR could be matched with a following "-int", "_int", "_I", "_IN", etc...
foreach $ltr (@LTR){
  foreach $sub ("_I","_IN", "_INT","_int","-I","-IN","-INT","-int"){
    $tmp=$ltr;
    $tmp.=$sub;
    if(($internal_elem{$tmp}==1) & ($ltr ne $tmp)){
      $inside{$tmp}=$ltr;
      $outside{$ltr}=$tmp;
      $matchcount++;
      next;
    }
  }
}

### Combination of the two previous steps
foreach $ltr (@LTR){
  foreach $sub ("_I","_IN", "_INT","_int","-I","-IN","-INT","-int"){
    foreach $subbis ("I","IN", "INT","int"){
      $tmp=$ltr;
      $tmp=~s/LTR/$subbis/;
      $tmp.=$sub;
      if(($internal_elem{$tmp}==1) & ($ltr ne $tmp)){
  $inside{$tmp}=$ltr;
  $outside{$ltr}=$tmp;
  next;
      }
    }
  }
}

### Next pass : we simplify element names by erasing LTR, I, IN, INT, or symbols like -_, etc...

%elem_liste_ltr=();
%elem_liste_int=();

foreach $ltr (@LTR){
  chomp $ltr;
  $tmp=$ltr;
  unless($tmp =~ /^LTR/){
    $tmp=~s/LTR//;
  }
  $tmp=~s/[-_]//g;
  $elem_liste_ltr{$ltr}=$tmp;
}

foreach $int (@INT){
  chomp $int;
  $tmp=$int;
  $tmp=~s/int//i;
  $tmp=~s/IN//i;
  $tmp=~s/[-_]I//;
  $tmp=~s/[-_]//g;
  $elem_liste_int{$int}=$tmp;
}

### And we try to match the new names
foreach $el1 (sort keys %elem_liste_int){
    ### An element previously recognized as external should not be internal
    next if(exists($outside{$el1}));
    foreach $el2 (sort keys %elem_liste_ltr){
      ### An element already recognized should not appear with another internal
      next if(exists($outside{$el2}) | exists($inside{$el2}));
      next if($el1 eq $el2);
      if($elem_liste_ltr{$el2} eq $elem_liste_int{$el1}){
  $matchcount++;
  unless(exists($inside{$el1})){
    $inside{$el1}=$el2;
  }else{
    $inside{$el1}.=":".$el2;
  }
  $outside{$el2}.=":".$el1;
      }
    }
}

$matchcount=0;
foreach $int (keys %inside){
  @G= split /:/,$inside{$int};
  $matchcount+=scalar(@G);
}

warn "$matchcount matches found in non-fuzzy phase\n";



if($fuzzy==1){
### In this case we go further in element name simplification, by
### allowing 1 or 2 letter differences between LTR and internal parts
  @patterns=("[a-zA-Z]","[0-9]","[a-zA-Z0-9][a-zA-Z0-9]");
  warn "\n#################################\n";
  warn "Phase 2 : Fuzzy matches\n";
  warn "#################################\n";

  foreach $pat (@patterns){
    foreach $el1 (sort keys %elem_liste_int){
  ### An element previously recognized as external should not be internal
  next if(exists($outside{$el1}));
  foreach $el2 (sort keys %elem_liste_ltr){
    ### An element already recognized  should not appear with another internal
    next if(exists($outside{$el2}) | exists($inside{$el2}));
    next if($el1 eq $el2);
    ### Is $el1 included in $el2, with a few mismatches in the end?
    if(($elem_liste_ltr{$el2}=~/$elem_liste_int{$el1}$pat$/) | ($elem_liste_int{$el1}=~/$elem_liste_ltr{$el2}$pat$/)){
      unless(exists($inside{$el1})){
        $inside{$el1}=$el2;
      }else{
        $inside{$el1}.=":".$el2;
      }
      $outside{$el2}.=":".$el1;
    }
  }
    }
  }
  $matchcount=0;
  foreach $int (keys %inside){
    @G= split /:/,$inside{$int};
    $matchcount+=scalar(@G);
  }
  warn "$matchcount total matches found after fuzzy phase\n";
}


### Sorting all elements to write each one only once
%uniques=();
foreach $int (@INT){
  unless(exists($inside{$int}) | exists($outside{$int})){
    $uniques{$int}=1;
  }
}

foreach $ltr (@LTR){
  unless(exists($outside{$ltr}) | exists($inside{$ltr})){
    $uniques{$ltr}=1;
  }
}

@uniques_sorted=sort keys %uniques;

### Writing results
foreach $a (sort keys %inside){
  print "$a\t$inside{$a}\n";
}
foreach $a (@uniques_sorted){
  print "$a\n";
}

warn "\n#################################\n";
warn scalar(@uniques_sorted)," elements found without match\n";
warn "#################################\n";


warn "\n\n#################################\n";
warn "Output file should be manually edited to take into account all specificities of the considered organism!\n";
warn "#################################\n";
