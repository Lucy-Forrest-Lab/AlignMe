#! /usr/bin/perl
use warnings;

print "Enter the name of the input file\n";
$infi =  <STDIN>;
chomp ($infi);
open (INFILE, "$infi") or die "no such file in directory\n";
print "INPUT file $infi has been successfully loaded\n";
print("\n");

print "enter the name of the output file\n";
$oufi = <STDIN>;
chomp ($oufi);
open (OUTFILE, ">$oufi") or die "cannot create $oufi\n";
print "$oufi has been created\n";

print OUTFILE "\@    yaxis  label \"averaged hydrophobicity\"\n";
print OUTFILE "\@    xaxis  label \"position\"\n";
print OUTFILE "\@    autoscale onread none\n";
print OUTFILE "\@    xaxis  tick on\n";
print OUTFILE "\@    xaxis  tick major 100\n";
print OUTFILE "\@    xaxis  tick minor ticks 1\n";
print OUTFILE "\@    yaxis  tick on\n";
print OUTFILE "\@    yaxis  tick major 1\n";
print OUTFILE "\@    yaxis  tick minor ticks 4\n";

$prof1 = 1;
$prof2 = 1;

while ($line = <INFILE>)
{
    if( substr($line, 0, 1) ne "#")
    {
	if ($line =~/(\d+)(\s+)(\S+)(\s+)(\S+)/)
        {
	    if ($3 ne 0)
	    {
		$prof1 = 0;
	    }
	    if ($5 ne 0)
	    {
		$prof2 = 0;
	    }
	    last;
	}
    }
}



############################
$set = 0;
$min = 0;

$check = $prof1;

seek(INFILE,0,0);

while ($line = <INFILE>)
{
    if( substr($line, 0, 1) ne "#")
    {
	if ($line =~/(\d+)(\s+)(\S+)(\s+)(\S+)/)
	{
	    if ($3<$min)
	    {
		$min = $3;
	    }
	    if ($5<$min)
	    {
		$min = $5;
	    }
	}
    }
}
###########################
seek(INFILE,0,0);

$max  = 0;
$xmin = 0;
$xmax = 0;

while ($line = <INFILE>)
{
    if( substr($line, 0, 1) ne "#")
    {
	
	if ($line =~/(\d+)(\s+)(\S+)(\s+)(\S+)/)
	{
	    if ($xmax < $1)
	    {
		$xmax = $1;
	    }
	    if ($max < $3)
	    {
		$max = $3;
	    }
	    if ($max < $5)
	    {
		$max = $5;
	    }
	}
    }
}
$difference = $max - $min;
$difference = $difference * 0.1;
$min = $min - $difference;
$max = $max + $difference;

#############################

seek(INFILE,0,0);

while ($line = <INFILE>) 
{
    if( substr($line, 0, 1) ne "#" && $line =~/[0-9]/)
    {
	if ($line =~/(\d+)(\s+)(\S+)/)
	{
	    $first_hydr = $3;
	}
	
	if ($first_hydr eq 0)
	{
	    if ($check eq 1)
	    {
		print OUTFILE "\&\n\@target G0.S$set\n";
		print OUTFILE "\@    s$set line linewidth 7\n";
		$set++;
		print OUTFILE "\@type xy\n";
		print OUTFILE "$1$2$min\n";
		$check = 0;
	    }
	    else
	    {
		print OUTFILE "$1$2$min\n";
	    }
	    
	}
	else 
	{
	    if ($check eq 0)
	    {
		print OUTFILE "\&\n\@target G0.S$set\n";
		print OUTFILE "\@    s$set line linewidth 1.5\n";
		$set++;
		print OUTFILE "\@type xy\n";
		print OUTFILE "$1$2$3\n";
		$check = 1
		}
	    else 
	    {
		print OUTFILE "$1$2$3\n";
	    }
	    
	}
    }
}

$remember_set = $set;

#################################

seek(INFILE,0,0);

print OUTFILE "&\n\n";
$check = $prof2;

while ($line = <INFILE>)
{
    if( substr($line, 0, 1) ne "#" && $line =~/[0-9]/)
    {
	if ($line =~/(\d+)(\s+)(\S+)(\s+)(\S+)/)
	{
	    $second_hydr = $5;
	}
	
	if ($second_hydr eq 0)
	{
	    if ($check eq 1)
	    {
		print OUTFILE "\&\n\@target G0.S$set\n";
		print OUTFILE "\@    s$set line linewidth 7\n";
		$set++;
		print OUTFILE "\@type xy\n";
		print OUTFILE "$1$2$min\n";
		$check = 0;
	    }
	    else
	    {
		print OUTFILE "$1$2$min\n";
	    }
	}
	else
	{
	    if ($check eq 0)
	    {
		print OUTFILE "\&\n\@target G0.S$set\n";
		print OUTFILE "\@    s$set line linewidth 1.5\n";
		$set++;
		print OUTFILE "\@type xy\n";
		print OUTFILE "$1$2$5\n";
		$check = 1
		}
	    else
	    {
		print OUTFILE "$1$2$5\n";
	    }
	}
    }
}

##########

for ($i = 0; $i<=$remember_set;$i++)
{
    print OUTFILE "\@    s$i line color  1\n";
}

for ($i = $remember_set; $i<$set;$i++)
{
    print OUTFILE "\@    s$i line color  2\n";
}
##########

$min = $min - $difference;

print OUTFILE "\@    world $xmin, $min, $xmax, $max\n";

close OUTFILE;
close INFILE;
