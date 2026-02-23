#!/usr/bin/perl -w

#======================================================#
#--- This Perl script automatically alters start    ---#
#--- and end dates in a input.geos file             ---#
#--- Called via:                                    ---#    
#--- changedate $template $year $startmonth         ---#
#--- changedate input.geos-template 2015 01         ---#
#---                                                ---#
#---                              r.r.b. 03/05/2012 ---#
#======================================================#

#------------------------------------------------#
#--- Read in template file, years and months                 
#------------------------------------------------#
$template = shift || "";
$YYYY = shift || "";
$m1 = shift || "";
 if ($template eq "" || $YYYY eq "" || $m1 eq ""){
   print "Something is missing: template file, year, month, or day field empty\n";
   exit 0;
 }
$d1 = sprintf("%02d", 1); 
$m1 = sprintf("%02d", $m1); 
if ($m1 < 12 ){
   $m2 = sprintf("%02d", $m1 + 1);
   $YYYYII = $YYYY;
}
else{
   $m2 = sprintf("%02d",  1);
   $YYYYII = $YYYY + 1;
}

$sdate = $YYYY.$m1.$d1;  
$edate = $YYYYII.$m2.$d1;

#------------------------------------------------#
#--- Read in template file, change date,
#--- write new file
#------------------------------------------------#
open (IN, "<$template");
chomp(@lines = <IN>);
close (IN);
$outfile = $template;
$outfile =~ s/-template//g;
#print "$outfile\n";
open (OUT, ">$outfile");

for (@lines){
  if ($_ =~ /^  start_date: /){
    print OUT "  start_date: [$sdate, 000000]\n";
  }
  elsif ($_ =~ /^  end_date: /){
    print OUT "  end_date: [$edate, 000000]\n";
  }
  else{
    print OUT "$_\n";
  }
}


#for (@lines){
#  if ($_ =~ /^Start YYYYMMDD/){
#    print OUT "Start YYYYMMDD, HHMMSS  : $sdate 000000\n";
#  }
#  elsif ($_ =~ /^End   YYYYMMDD/){
#    print OUT "End   YYYYMMDD, HHMMSS  : $edate 000000\n";
#  }
#  else{
#    print OUT "$_\n";
#  }
#}

close (OUT);
