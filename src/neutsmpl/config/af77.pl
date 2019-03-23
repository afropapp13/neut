eval 'exec perl -w -S $0 ${1+"$@"}'
    if 0;
#
# The above is a trick to find perl wherever it is installed.
#
# Go around Absoft f77 bug which prevents argumets like -DVER=4.09.
#
# Bugs:
#  - -D NAME=value (with space after D) is not processed correctly.
#    Same thing for -I option if space is allowed there too.
#  - Absoft f77, cpp and rm are hard-coded. Should be customizable,
#    maybe by environmant variables.
#  - Only one Fortran file per compilation call is allowed,
#    maybe Absoft allows to compile more files at once?
#  - Options with space (like -o file) are not supported.
#
# $Id: af77.pl,v 1.1.1.1 2016/11/04 20:29:56 mdunkman Exp $
#
# $Log: af77.pl,v $
# Revision 1.1.1.1  2016/11/04 20:29:56  mdunkman
# Import of NEUT 5.3.3 v1r27 for use with T2KReWeight v1r27 (Winter OA 2017)
#
# Revision 1.2  1999/02/19 19:22:41  habig
# If the file is a .f or .for, only pass @f77_args.  cpp wont be run
# anyway, and some of the cpp args would cause f77 to die.
#
# Revision 1.1  1998/11/14 04:31:03  tomba
# Added a wrapper for Absoft f77 compiler which works around a bug
# that prevents use of most options like -DNAME=VALUE.
# The wrapper af77.pl takes out all cpp options,
# calls cpp and then calls f77 without cpp options.
#


#
# Helps catch some mistakes in the code:
#
use strict


#
# Those can be adjusted on a per-site basis.
#
$f77_cmd = "af77";
$cpp_cmd = "/lib/cpp";
$rm_cmd  = "/bin/rm";


#
# No more adjustable things below ---------------------------------------------
# 

@f77_args = ();
@cpp_args = ();
$file     = "";

# 
# Split command line arguments into those for cpp, those for f77
# and file name.
#
foreach $argument (@ARGV)
{
    # Make sure there is only one file and only as a last argument.
    if ($file ne "")
    {
	die "Argument ($argument) after file ($file)";
    }
    if ($argument =~ /^-/)
    {
	if ($argument =~ /^-[DI]/)
	{
	    @cpp_args = (@cpp_args, $argument);
	}
	else
	{
	    @f77_args = (@f77_args, $argument);
	}
    }
    else
    {
	$file = $argument
    }
}


#
# *.F and *.FOR are parsed by cpp, others are not.
#
if ($file =~ /[.]F$/ || $file =~ /[.]FOR$/)
{
    $base_name = substr($file, 0, rindex($file, "."));
    # The file name has to be *.f or *.for, 
    # otherwise Absoft f77 would not work right.
    # *.cpp in particular would not work.
    $cpp_file = $base_name . ".f";

    # This is Absoft f77 calls cpp.
    # Options -E and -o are ignored by cpp as far as I know.
    @cpp_call = ($cpp_cmd, "-traditional", @cpp_args, "-E", 
		 $file, "-o", $cpp_file);
    if (system @cpp_call) 
    {
	die "$cpp_cmd failed: $?";
    }

    @f77_call = ($f77_cmd, @f77_args, $cpp_file);
    if (system @f77_call)
    {
	die "$f77_cmd failed: $?";
	# Bug: temporary *.f file is not removed in such case.
    }

    @rm_call = ($rm_cmd, $cpp_file);
    if (system @rm_call)
    {
	die "$rm_cmd failed: $?";
    }
}
else
{
    # Just call f77 as if there was no this script.
    # but use only f77 args, since cpp won't be run anyway, and cpp
    # args can flake out the f77 shell script.
    @f77_call = ($f77_cmd, @f77_args, $file);
    if (system @f77_call)
    {
	die "raw $f77_cmd failed: $?";
    }
    # Could add removing of the annoying temporary .cpp file,
    # but none should be generated for *.f and *.for files.
    # Also could remove potentially fatal cpp options.
}
