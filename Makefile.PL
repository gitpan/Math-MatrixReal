use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'		=> 'Math::MatrixReal',
    'DISTNAME'		=> 'Math-MatrixReal',
    'VERSION_FROM'	=> 'MatrixReal.pm',
    'dist'		=> { COMPRESS=>"gzip -9", SUFFIX=>"gz" },
    'LIBS'		=> [''],   # e.g., '-lm' 
    'DEFINE'		=> '',     # e.g., '-DHAVE_SOMETHING' 
    'INC'		=> ''      # e.g., '-I/opt/pkg/perl5.003/dist' 
);
