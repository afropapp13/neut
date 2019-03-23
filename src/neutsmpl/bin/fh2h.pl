eval 'exec perl -w -S $0 ${1+"$@"}'
    if 0;
#
# $Id: fh2h.pl,v 1.1 2008-04-13 16:17:19 skrep Exp $
#
# Author: Tomasz Barszczak
#         tomba@hepxvt.ps.uci.edu
#         http://www.ps.uic.edu/~tomba/
#
# Converts fortran header files 
# which contain definitions of parameters and common blocks
# into C header files.
# Does not work with all valid fortran code, only with subset of it.
# It was tested on SK fortran include files.
#
# Why not use f2c?
#  - it lost information about PARAMETER definitions
#  - COMMON block structures were anonymous so that 
#    one could not make pinters to them
#  - multi-dimensional arrays were made flat one-dimensional arrays
#  - C preprocessor directives caused it to exit with error
#  - more ... (I don't remember)
#
# This is my first perl program, so please bear with it.
# 
# Bugs: (many, some are documented in the code)
#  - Fortran is insensitive to spaces, should be same with this script
#    but there could be bugs with it
#  - Careless about alignment of variables in common
#  - should accept command line options for tunable parameters
#  - accepts more than one file name in the command line
#  - Does not do complex parameters (others? logical?)
#  - does not check for declared variables which are not
#    in common block or PARAMETERs (i.e. does not warn about them).
#  - Ignores types of PARAMETERs other than those initialized with strings
#  - power (**) in constants is not translated
#  - chould combine char size with dimensions, and make dimensions
#    list of sizes instead of text
#  - real parameter (somename=10) will be integer in C
#  - does not translate differences in floating point constant syntax
#    (double and quadruple precision for example)
#  - ...
#
# Check http://www.ps.uci.edu/~tomba/fh2h/ for the newest version.
# (SK group members should check the CVS repository instead.)
#
#
# $Log: not supported by cvs2svn $
# Revision 1.2  1999/03/16 05:39:06  mcgrew
# Merge the fh2h.pl bugfixes from the ichikilo branch.
#
# Revision 1.1.2.1  1999/02/16 01:05:16  mcgrew
# Extra CPP directives must be passed through to the C.h file.  There
# may be a subtle problem with "defines".  If that rears it's ugly head,
# then it may be necessary to ignore define, but that is likely to cause
# lots of problems.
#
# Revision 1.1  1998/09/15 21:15:20  mcgrew
# This is taken from atmpd/src/programs/chgen/fh2h.pl.  It isn't really
# part of chgen, so it's OK to move it out of that directory.  This is
# exactly Tomasz's program, except that I added a line at the begin that
# is suppose to make it more portable.  We'll see.
#
# Revision 1.10  1998/06/18 18:07:12  tomba
# Made for and foreach loops more portable by moving 'my' in front of them.
#
# Revision 1.9  1998/06/01 04:09:05  tomba
# Added DATA statements and more.
#
# Revision 1.8  1998/05/28 09:18:52  tomba
# Protected contents of C comments.
#
# Revision 1.7  1998/05/28 06:41:11  tomba
# Forgot to remove some debug output.
#
# Revision 1.6  1998/05/28 06:37:36  tomba
# Fixed "TYPE[^*]" bug.
#
# Revision 1.5  1998/05/28 05:37:56  skrep
# Added support for more kinds of fortran code.
#
# Revision 1.4  1998/05/27 05:46:56  tomba
# Oops... forgot to change back default behavior (-w and $die_on_error).
#
# Revision 1.3  1998/05/27 05:44:48  tomba
# Added C++ support and pointers to common blocks.
#
# Revision 1.2  1998/05/27 03:56:21  skrep
# Fixed warnings, added #ifndef, error and info routines.
#
# Revision 1.1  1998/05/26 06:33:01  tomba
# First SK release.
#


#------ Tunable parameters ----------------------------------------------------


#
# Helps catch some mistakes in the code:
#
use strict;


#
# Whether to produce debug comments
#
my $debug_comments = 0;


#
# Whether to produce other information comments
#
my $info_comments = 1;


#
# Whether to preserve original fortran comments and empty lines
#
my $save_comments = 1;


#
# Whether to terminate the program on error or continue.
#
my $die_on_error = 1;


# How to do this?
# STDOUT->autoflush(1);


#------ Static data -----------------------------------------------------------


#
# Translation of fortran types to C types
#
my %types = 
    (
     "INTEGER[*]4"                     => "int",
     "INTEGER[*]2"                     => "short",
     "INTEGER"                         => "int",
     "REAL[*]16"                       => "long double",
     "REAL[*]8"                        => "double",
     "REAL[*]4"                        => "float",
     "REAL"                            => "float",
     "CHARACTER[*][(][a-zA-Z0-9_]+[)]" => "char_star",
     "CHARACTER[*][0-9]+"              => "char_star",
     "CHARACTER[*][(][*][)]"           => "char_err",
     "CHARACTER"                       => "char",
     "BYTE"                            => "char",
     "LOGICAL[*]4"                     => "int",
     "LOGICAL[*]2"                     => "short",
     "LOGICAL[*]1"                     => "char",
     "LOGICAL"                         => "int",
     "DOUBLEPRECISION"                 => "double",
     "COMPLEX[*]8"                     => "struct {float r,i;}",
     "COMPLEX[*]16"                    => "struct {double r,i;}",
     "COMPLEX[*]32"                    => "struct {long double r,i;}",
     "DOUBLECOMPLEX[*]8"               => "struct {double r,i;}",
     "COMPLEX"                         => "struct {float r,i;}",
     );


#
# What to do about simple fortran statements
# Bug: e.g. END can match before ENDIF for ENDIF statement.
#   This is because hash tables are unordered
#
my %statements =
    (
     "PROGRAM"                  => "err",
     "SUBROUTINE"               => "err",
     "FUNCTION"                 => "err",
     "ENTRY"                    => "err",
     "IMPLICITNONE"             => "warn",
     "IMPLICIT"                 => "err",
     "INCLUDE"                  => "warn",
     "BLOCKDATA"                => "unimp",
     "SAVE"                     => "ignore",
     "POINTER"                  => "unimp",
     "STRUCTURE"                => "unimp",
     "ENDSTRUCTURE"             => "unimp",
     "UNION"                    => "unimp",
     "ENDUNION"                 => "unimp",
     "RECORD"                   => "unimp", # Should be "err"?
     "EQUIVALENCE"              => "unimp", # Could use union?
     "INTRINSIC"                => "err",
     "IF"                       => "err",
     "ELSEIF"                   => "err",
     "ELSE"                     => "err",
     "ENDIF"                    => "err",
     "DOWHILE"                  => "err",
     "DO"                       => "err",
     "ENDDO"                    => "err",
     "GOTO"                     => "err",
     "CONTINUE"                 => "err",
     "RETURN"                   => "err",
     "CALL"                     => "err",
     "OPEN"                     => "err",
     "CLOSE"                    => "err",
     "REWIND"                   => "err",
     "READ"                     => "err",
     "WRITE"                    => "err",
     "PRINT"                    => "err",
     "TYPE"                     => "err",
     "FORMAT"                   => "warn",
     "PAUSE"                    => "err",
     "STOP"                     => "err",
     "END"                      => "err",
     # I don't know those:
     "ACCEPT"                   => "unimp",
     "ASSIGN"                   => "unimp",
     "AUTOMATIC"                => "unimp",
     "BACKSPACE"                => "unimp",
     "DECODE"                   => "unimp",
     "ENCODE"                   => "unimp",
     "MAP"                      => "unimp",
     "ENDMAP"                   => "unimp",
     "ENDFILE"                  => "unimp",
     "EXTERNAL"                 => "unimp",
     "INQUIRE"                  => "unimp",
     "NAMELIST"                 => "unimp",
     "OPTIONS"                  => "unimp",
     "PRAGMA"                   => "unimp",
     "STATIC"                   => "unimp",
     "VIRTUAL"                  => "unimp",
     "VOLATILE"                 => "unimp",
     );


# 
# C and C++ reserved words.
# Complete list for C (from ANSI C book),
# and selected frequent words for C++ (from memory).
#
my @reserved = 
    ("auto", "break", "case", "char", "const", "continue", "deault", 
     "do", "double", "else", "enum", "extern", "float", "for", 
     "goto", "if", "int", "long", "register", "return", "short", 
     "signed", "sizeof", "static", "struct", "switch", "typedef", 
     "union", "unsigned", "void", "volatile", "while",
     "new", "class", "private", "public", "protected", "delete",
     "virtual", "friend");


#------ Global variables ------------------------------------------------------


my $line;             # current line
my $next_line;        # next line
my %varhash;          # variable types
my %dimhash;          # variable dimensions
my %charhash;         # size of string variables
my %commonhash;       # common blocks
my %datahash;         # data statements
my @common_sequence;  # preserve sequence of common blocks
my @data_sequence;    # preserve sequence of data statements
my $ifdef_name;       # converted name of the input file


#------ Utility subroutines ---------------------------------------------------


#
# Make a comment with protected "/*" and "*/" in the comment
#
sub tcmt
{
    my $comment = $_[0];
    $comment =~ s(^/)( /);       # protect first slash
    $comment =~ s(/$)(/ );       # protect last slash
    $comment =~ s(/[*])(/ *)g;   # protect "/*" -> "/ *"
    $comment =~ s([*]/)(* /)g;   # protect "*/" -> "* /"
    print "/*$comment*/\n";
}


#
# Print error and exit conditionally
#
sub terror
{
    my $message = $_[0];
    tcmt("ERROR: $message");
    tcmt("Line:$line");
    if ($die_on_error)
    {
	die "ERROR: $message\nLine: $line\n";
    }
    else
    {
	print STDERR "ERROR: $message\nLine:$line\n";
    }
}


#
# Print warning
#
sub twarning
{
    my $message = $_[0];
    tcmt("WARNING: $message");
    print STDERR "WARNING: $message\n";
    if ($line ne "")
    {
	tcmt("Line:$line");
	print STDERR "Line:$line\n";
    }
}


#
# Make information comment
#
sub tinfo
{
    if ($info_comments)
    {
	my $message = $_[0];
	tcmt("$message");
    }
}


#
# Make debug comment
#
sub tdebug
{
    if ($debug_comments)
    {
	my $message = $_[0];
	tcmt("$message");
    }
}


#
# remove white space from variable name or other name
# I should have done it with one simple regular expression
#
sub remove_whitespace
{
    my $varname = "";
    my $i;
    for ($i=0; $i<length($_[0]); $i++)
    {
	my $character = substr($_[0], $i, 1);
	unless ($character =~ /\s/)
	{
	    $varname = $varname . $character;
	}
    }
    return $varname;
}


#
# split line of variable names on commas, 
# but don't touch commas in dimension declarations
#
sub split_varline
{
    my $typeline = "";
    my $inside_paren = 0;
    my $i;
    for ($i=0; $i<length($_[0]); $i++)
    {
	my $character = substr($_[0], $i, 1);
	if ($character eq "(")
	{
	    $inside_paren = 1;
	}
	if ($character eq ")")
	{
	    $inside_paren = 0;
	}
	if ($character eq ",")
	{
	    if (!$inside_paren)
	    {
		$character = ";";
	    }
	}
	$typeline = $typeline . $character;
    }
    # Split into individual variable definitions
    return split(";", $typeline);
}


#
# Translate "(...)*N" into "(N,...)"
#
sub star_to_dim
{
    my $dimf = $_[0];
    # remove "("
    $dimf = substr($dimf, 1);
    (my $reg_dim, my $star_dim) = split("[)]", $dimf);
    # remove "*"
    $star_dim = substr($star_dim, 1);
    return "($star_dim,$reg_dim)";
}


#
# Convert Fortran dimensions into C dimensions:
# "(dim1,lo:hi,dim3)*size" -> "[dim3][(hi)-(lo)+1][dim1][size]"
#
sub dim_f2c
{
    my $dimf = $_[0];
    if (substr($dimf, length($dimf)-1) ne ")")
    {
	$dimf = star_to_dim($dimf);
    }
    # remove parens (first and last character)
    $dimf = substr($dimf, 1, length($dimf)-2);
    my @dimf_list = split(",", $dimf);
    my $dimc = "";
    my $dimf_elem;
    foreach $dimf_elem (@dimf_list)
    {
	# Split low:high
	my $dimc_elem;
	(my $dim_lo, my $dim_hi) = split(":", $dimf_elem);
	if ($dim_hi)
	{
	    $dimc_elem = "($dim_hi)-($dim_lo)+1";
	}
	else
	{
	    $dimc_elem = $dimf_elem;
	}
	$dimc_elem = uc($dimc_elem);
	$dimc = "[$dimc_elem]" . $dimc;
    }
    return "$dimc";
}


#
# Split " varname (dim1, dim2)" into ("varname", "[dim2][dim1]")
#
sub split_var_dim
{
    # Extract array dimensions
    (my $varname2, my $varpar) = split("[(]", $_[0]);
    if ($varpar)
    {
	$varpar = "(" . $varpar;
	$varpar = remove_whitespace($varpar);
	$varpar = dim_f2c($varpar);
    }
    $varname2 = remove_whitespace($varname2);
    return ($varname2, $varpar)
}


#
# Add variable to dimension hash
#
sub add_to_dim_hash
{
    my $index = lc($_[0]);
    my $hash_dimvalue = $_[1];
    if ($hash_dimvalue)
    {
	if ((exists $dimhash{$index}) &&
	    ($dimhash{$index} ne $hash_dimvalue))
	{
	    terror("Redefined dimension: " . 
		   "$index was '$dimhash{$index}' " .
		   "now '$hash_dimvalue'");
	}
	$dimhash{$index} = $hash_dimvalue;
	tdebug("Added dim: $index: $dimhash{$index}");
    }
}


#
# Align most type lengths
#
sub nice_type
{
    my $type = $_[0];
    my $length_diff = length("double") - length($type);
    if ($length_diff > 0)
    {
	$type .= " " x $length_diff;
    }
    return $type;
}


#
# Returns type of the variable from the hash or if not found
# deduces default implicit type and prints warning.
#
sub get_type
{
    my $varname = $_[0];
    my $type = $varhash{$varname};
    if (! exists $varhash{$varname})
    {
	my $implicit;
	twarning("Untyped (undeclared) variable $varname");
	if (substr($varname, 0, 1) =~ /[i-n]/)
	{
	    $implicit = "int";
	}
	else
	{
	    $implicit = "float";
	}
	$type = "IMPLICIT $implicit";
    }
    return $type;
}


#
# Apend underscore and print warning if word is reserved
#
sub check_reserved
{
    my $varname = $_[0];
    my $resword;
    foreach $resword (@reserved)
    {
	if ($varname eq $resword)
	{
	    twarning("$varname is reserved word in C/C++");
	    $varname = $varname . "_";
	}
    }
    return $varname;
}


#
# Converts fortran constants into C constants
#   cons tant -> CONSTANT
#   'Some "\text"' -> "Some \"\\text\""
# Bug: does not process special fortran sequences (e.g. '', \t) 
# Bug: should look at the constant to find it is string.
# Bug: converts all ', not just the first and last ones.
#
sub const_f2c
{
    my $type = $_[0];
    my $fconst = $_[1];
    if ($type eq "char")
    {
	my $cconst = $fconst;
	$cconst =~ s/[\\]/\\\\/g;  # Double backslashes
	$cconst =~ s/"/\\"/g;      # Escape double quotes
        $cconst =~ s/'/"/g;        # Change single quotes to double quotes
        # Does nothing except fixes emacs highliting confused by the quotes
        $cconst =~ s/'/''/g;
	return $cconst;
    }
    else
    {
	my $cconst = uc(remove_whitespace($fconst)); 
	return $cconst;
    }
}


#------ The main code ---------------------------------------------------------


#
# Make the preamble
#

print "/*\n";
print " * Generated automatically by fh2h.pl\n";
print " * !!! DO NOT EDIT !!!\n";
print " * Edit the original fortran header file instead\n";
print " * or fix fh2h.pl if there is a translation bug.\n";
print " */\n";
print "\n";
print "\n";


if ($#ARGV+1)
{
    $ifdef_name = uc($ARGV[0]);
    # Translate dots, slashes, dashed into underscores
    # Bug: should translate more characters, 
    # i.e. all all non [0-9][a-z][A-Z][_]
    $ifdef_name =~ tr(./-)(_);
    $ifdef_name = "FH2H_" . $ifdef_name;
    print "#ifndef $ifdef_name\n";
    print "#define $ifdef_name\n";
    print "\n";
    print "\n";
}


print "#ifdef __cplusplus\n";
print "extern \"C\" {\n";
print "#endif\n";
print "\n";
print "\n";


print "#ifndef IMPLICIT\n";
print "#define IMPLICIT  /* Only to point out implicit types */\n";
print "#endif\n";
print "\n";
print "\n";


print 
"/*------ fortran header (without commons and data statements) ----------*/\n";
print "\n";


#
# The main loop
#
$next_line = <>;
INPUT_LINE: while (1)
{
    if ($next_line eq "")
    {
	last INPUT_LINE;
    }

    $line = $next_line;
    $next_line = <>;
    unless (defined $next_line)
    {
	$next_line = "";  # to avoid perl warning
    };

    # Remove the newline character
    chomp($line);

    # Ignore #include CPP directives
    if (substr($line, 0, 8) eq "#include" )
    {
	twarning("Ignoring include: $line");
	print "/* Include the corresponding *C.h file in your program */\n";
	next INPUT_LINE;
    }

    # Pass through other CPP directives
    # No, better reject them too
    # Absolutely NOT!!!!! The directives must be passed.
    if (substr($line, 0, 1) eq "#" )
    {
	print "$line\n";
	# tinfo("Ignoring cpp directive: $line");
	next INPUT_LINE;
    }

    # Convert full line comments
    if ($line =~ /^[CcDd*!]/)
    {
	if ($save_comments)
	{
	    # Remove comment character
	    my $comment = substr($line, 1);
	    tcmt("$comment");
	}
	next INPUT_LINE;
    }

    # Convert inline comments 
    if ($line =~ /!/)
    {
	# remove everything after "!"
	($line, my @rest) = split("!", $line, 2);
	# Bug: sensitive to "!" in quotes too, but this does not
        # happen usually in header files.
	if ($save_comments)
	{
	    tcmt("@rest");
	}
    }

    # Join continuation lines
    # Bug: Comments and blank lines before contiuation lines break this
    #
    # 5 spaces followed by a nonspace, nonzero
    if ($next_line =~ /^[ ][ ][ ][ ][ ][^ 0]/)
    {
	$next_line = $line . " " . substr($next_line, 6);
	next INPUT_LINE;
    }
    # tab followed by nonzero digit
    if ($next_line =~ /^[\t][1-9]/)
    {
	$next_line = $line . " " . substr($next_line, 2);
	next INPUT_LINE;
    }
    
    # Replace tabs with 8 spaces (incorrect but OK here)
    # Bug: will replace tabs in strings, but it is OK here
    #      should make a new variable instead and keep $line intact
    $line =~ s/[\t]/        /g;

    my $line_nolabel;
    # Remove line labels
    # This must have some bug. It even does not check if the 'label' is numeric
    if ((length($line) > 0) && (substr($line, 0, 1) ne "	"))
    {
	$line_nolabel = substr($line, 6);
    }
    else
    {
	$line_nolabel = $line;
    }

    # Remove spaces (also from strings, but it is OK here)
    my $line_nospaces = $line_nolabel;
    $line_nospaces =~ s/[ ]//g;

    # Skip empty lines
    if ((length($line_nolabel) == 0) ||
	(length(remove_whitespace($line_nolabel)) == 0))
    {
	if ($save_comments)
	{
	    print "$line\n";
	}
	next INPUT_LINE;
    }

    # 
    # Parse PARAMETER lines
    #
    if ($line_nospaces =~ /^\s*PARAMETER/i)
    {
	# Select stuff inside (...)
        # Bug: Multiple or nested (()) would confuse it
	# Should extract first and last paren instead
	my $parline = substr($line_nolabel,
			     index($line_nolabel, "(") + 1,
			     index($line_nolabel, ")") - 
			     index($line_nolabel, "(") - 1);
	# Split into individual parameter definitions
	my @params = split(",", $parline);
	my $parexpr;
	foreach $parexpr (@params)
	{
	    (my $parnam, my $parval) = split("=", $parexpr);
	    $parnam = uc(remove_whitespace($parnam));
	    my $type;
	    # This test shold go to const_f2c
	    if (($parval =~ /'/) || ($parval =~ /'/))
	    # The second condition is only to fix broken emacs highliting
	    {
		$type = "char";
	    }
	    else
	    {
		$type = "nonchar";
	    }
	    $parval = const_f2c($type, $parval);
	    print "#define $parnam ($parval)\n"
	}
	tdebug("Param: $line");
	next INPUT_LINE;
    }

    #
    # Parse the type lines
    #
    my $ftype;
    foreach $ftype (keys %types)
    {
	if ($line_nospaces =~ /^\s*($ftype)([^*].*)/i)
	{
	    my $char_dim;
	    my $type_part = $1;
	    my $variable_part = $2;
	    if ($types{$ftype} eq "char_err")
	    {
		terror("Illegal variable length $type_part");
		next INPUT_LINE;
	    }
	    if ($types{$ftype} eq "char_star")
	    {
		$type_part =~ /(CHARACTER[*])(.*)/i;
		$char_dim = $2;
		# Remove parens
		if (substr($char_dim, 0, 1) eq "(")
		{
		    $char_dim = substr($char_dim, 1, length($char_dim)-2);
		}
		$char_dim = uc($char_dim);
	    }
	    # now we are left with the list of variables
	    my $typeline2 = $variable_part;

	    # Split into individual variable definitions
	    my @varlist = split_varline($typeline2);

	    my $varname_par;
	    foreach $varname_par (@varlist)
	    {
		# Extract array dimensions
		(my $varname, my $varpar) = split_var_dim($varname_par);

		add_to_dim_hash($varname, $varpar);

		my $index = lc($varname);
		my $hash_value = $types{$ftype};

		# Add to character array size hash
		if ($hash_value eq "char_star")
		{
		    $hash_value = "char";
		    if (exists $charhash{$index} && 
			($charhash{$index} ne $char_dim))
		    {
			terror("Redefined char array size: " .
			       "$varname was *$charhash{$index} " .
			       "now *$char_dim");
		    }
		    $charhash{$index} = $char_dim;
		    tdebug("Added size: $index: $charhash{$index}");
		}

		# Add variable type to hash
		if (exists $varhash{$index} && 
		    ($varhash{$index} ne $hash_value))
		{
		    terror("Redefined type: $varname was $varhash{$index} " .
			   "now $hash_value");
		}
		$varhash{$index} = $hash_value;
		tdebug("Added type: $index: $varhash{$index}");
	    }
	    tdebug("Type: $line");
	    next INPUT_LINE;
	}
    }
    
    #
    # parse DIMENSION lines
    #
    if ($line_nospaces =~ /^\s*(DIMENSION)(.*)/i)
    {
	# Remove the word "DIMENSION"
	my $eaten_line = $2;
	my @varlist = split_varline($eaten_line);
        my $varname_par;
	foreach $varname_par (@varlist)
	{
	    (my $varname, my $varpar) = split_var_dim($varname_par);
	    add_to_dim_hash($varname, $varpar);
	}
	tdebug("Dimension: $line");
	next INPUT_LINE;
    }

    #
    # Parse COMMON blocks
    #
    if ($line_nospaces =~ /^\s*COMMON/i)
    {
	(my $nothing, my $common_name2, my $commonline2) = 
                                                     split("/", $line_nolabel);
	undef $nothing; # To avoid perl warning

	my $common_name = remove_whitespace($common_name2);
	
	# Split into individual variable definitions
	my @varlist = split_varline($commonline2);

	my @clean_varlist = ();
	my $varname_par;
	foreach $varname_par (@varlist)
	{
	    # Extract array dimensions
	    (my $varname, my $varpar) = split_var_dim($varname_par);
	    
	    @clean_varlist = (@clean_varlist, lc($varname));
	    
	    add_to_dim_hash($varname, $varpar);
	}

	$common_name = lc($common_name);
	if (exists $commonhash{$common_name})
	# Should check if they are identical (by variable names)
	{
	    terror("Redefined common: " . 
		   "$common_name was (@{$commonhash{$common_name}}) " .
		   "now (@clean_varlist)");
	}
	else
	{
	    @common_sequence = (@common_sequence, $common_name);
	}
	$commonhash{$common_name} = [ @clean_varlist ];
	tdebug("Added common: $common_name: @{$commonhash{$common_name}}");
	tdebug("Common: $line");
	tinfo("common $common_name was here");
	next INPUT_LINE;
    }

    #
    # Parse DATA statements
    # BUG: Does not work right with character or multi-dimensional data
    # BUG: I have one of two choices:
    #  1) spaces in strings would be lost.
    #  2) malicious spaces would break this code, 
    #  My choice is 2).
    #  No, wait, maybe I have fixed it completely now.
    #
    if ($line_nospaces =~ m{^\s*(DATA)(.*?)/(.*)/}i)
    {
	my $data_name = $2;
	(my $nothing, my $dataline2, my $nothing2) = split("/", $line_nolabel);
	undef $nothing ; # To avoid perl warning
	undef $nothing2; # To avoid perl warning

	# Split into individual constants
	my @datalist = split(",", $dataline2);

	$data_name = lc($data_name);
	if (exists $datahash{$data_name})
	{
	    terror("Redefined data: " . 
		   "$data_name was (@{$datahash{$data_name}}) " .
		   "now (@datalist)");
	}
	else
	{
	    @data_sequence = (@data_sequence, $data_name);
	}
	$datahash{$data_name} = [ @datalist ];
	tdebug("Added data: $data_name: @{$datahash{$data_name}}");
	tdebug("Data: $line");
	tinfo("data statement for $data_name was here");
	next INPUT_LINE;
    }

    #
    # Parse other simple statements
    #
    my $statement;
    foreach $statement (keys %statements)
    {
	if ($line_nospaces =~ /^\s*($statement)/i)
	{
	    my $match = $1;
	    if ($statements{$statement} eq "err")
	    {
		terror("'$match' does not belong in header file");
		next INPUT_LINE;
	    }
	    elsif ($statements{$statement} eq "warn")
	    {
		twarning("'$match' does not belong in header file");
		next INPUT_LINE;
	    }
	    elsif ($statements{$statement} eq "ignore")
	    {
		tinfo("Ignoring: $line");
		next INPUT_LINE;
	    }
	    elsif ($statements{$statement} eq "unimp")
	    {
		terror("'$match' is not implemented");
		next INPUT_LINE;
	    }
	    else
	    {
		terror("BUG: impossible action '$statements{$statement}'");
		next INPUT_LINE;
	    }
	}
    }

    #
    # Unrecognized statements
    #
    terror("Unrecognized statement");
}


$line = "";  # For twarning()


print "\n";
print "\n";
print 
"/*------ common blocks -------------------------------------------------*/\n";


#
# Make structures for common blocks
#
my $common_name;
foreach $common_name (@common_sequence)
{
    my $common_type_name = $common_name . "_common";
    my $common_underscore = $common_name . "_";
    my $common_addr = "&" . $common_underscore;
    print "\n";
    print "extern struct $common_type_name {\n";
    my $varname;
    foreach $varname (@{$commonhash{$common_name}})
    {
	my $type = get_type($varname);
	my $aligned_type = nice_type($type);
	my $safe_varname = check_reserved($varname);
	print "  $aligned_type $safe_varname";
	if ($dimhash{$varname})
	{
	    print "$dimhash{$varname}";
	}
	if ($charhash{$varname})
	{
	    print "[$charhash{$varname}]";
	}
	print ";\n";
    }
    print "} $common_underscore;\n";
    print "#ifndef NO_EXTERN_COMMON_POINTERS\n";
    # Bug: not checking for reserved words.
    print "extern struct $common_type_name *$common_name;\n";
    print "#endif\n";
    print "#ifdef STATIC_COMMON_POINTERS\n";
    # Should be #ifndef NO_STATIC_COMMON_POINTERS, but done this way
    # for historical reasons (for BU people, especially Chris).
    # Bug: not checking for reserved words.
    print "static struct $common_type_name *$common_name = $common_addr;\n";
    print "#endif\n";
}


print "\n";
print "\n";
print 
"/*------ data statements -----------------------------------------------*/\n";


print "\n";
print "\n";
print "#ifndef NO_STATIC_DATA\n";


#
# Make static variables for data statements
#
my $data_name;
foreach $data_name (@data_sequence)
{
    my $type = get_type($data_name);
    my $safe_data_name = check_reserved($data_name);
    print "\n";
    if ($charhash{$data_name})
    # Bug: should check dimensions instead,
    # after combining character hash with dimension hash
    {
	# Convert character arrays to pointers to character
        # i.e. traditional string type in C.
	$safe_data_name = "*" . $safe_data_name;
	# Make maximum string length visible in C
	my $max_size_name = uc($data_name . "_MAX_LENGTH");
	my $max_size = uc($charhash{$data_name});
	print "#define $max_size_name ($max_size)\n";
    }
    my $aligned_type = nice_type($type);
    print "static $aligned_type $safe_data_name";
    if ($dimhash{$data_name})
    {
	print "$dimhash{$data_name}";
	print " = {";
	my $first_const = 1;
	my $const;
	foreach $const (@{$datahash{$data_name}})
	{
	    if ($first_const)
	    {
		$first_const = 0;
	    }
	    else
	    {
		print ", ";
	    }
	    # Bug: would translate byte as string
	    my $const_c = const_f2c($type, $const);
	    print "$const_c";
	}
	print "}";
    }
    else
    {
	print " = ";
	my $const = ${$datahash{$data_name}}[0];
        my $const_c = const_f2c($type, $const);
        print "$const_c";
    } 
    print ";\n";
}


print "\n";
print "\n";
print "#endif  /* #ifndef NO_STATIC_DATA */\n";


#
# Print out symbol tables
#
if ($debug_comments)
{
    print "\n";
    print "\n";
    print 
"/*------ symbol tables -------------------------------------------------*/\n";

    print "\n/*** Type hash: ***/\n";
    my $varname;
    foreach $varname (keys %varhash)
    {
	tcmt(" $varname: $varhash{$varname} ");
    }
    print "\n/*** Dimension hash: ***/\n";
    foreach $varname (keys %dimhash)
    {
	tcmt(" $varname: $dimhash{$varname} ");
    }
    print "\n/*** Character hash: ***/\n";
    foreach $varname (keys %charhash)
    {
	tcmt(" $varname: $charhash{$varname} ");
    }
    print "\n/*** Common hash: ***/\n";
    foreach $varname (keys %commonhash)
    {
	tcmt(" $varname: @{$commonhash{$varname}} ");
    }
    print "\n/*** Data hash: ***/\n";
    foreach $varname (keys %datahash)
    {
	tcmt(" $varname: @{$datahash{$varname}} ");
    }
}


#
# "Post-amble"
#

print "\n";
print "\n";
print 
"/*------ end of fortran header -----------------------------------------*/\n";


print "\n";
print "\n";
print "#ifdef __cplusplus\n";
print "}\n";
print "#endif\n";



if (defined $ifdef_name)
{
    print "\n";
    print "\n";
    print "#endif  /* #ifndef $ifdef_name */\n";
}
