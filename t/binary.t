BEGIN { $| = 1; print "1..4\n"; }
END {print "not ok 1\n" unless $loaded;}
use File::Spec;
use lib File::Spec->catfile("..","lib");
use Math::MatrixReal;
$loaded = 1;

do 'funcs.pl';

$matrix = Math::MatrixReal->new_from_string(<<"MATRIX");
[ 1 0 0 0 1 ]
[ 0 2 0 0 2 ]
[ 0 0 3 0 0 ]
[ 0 0 0 4 0 ]
[ 0 0 0 0 5 ]
MATRIX
ok(1,  ! $matrix->is_binary );
$matrix->one();
ok(2, $matrix->is_binary );
$matrix->zero();
ok(3, $matrix->is_binary );
$matrix = Math::MatrixReal->new_from_string(<<"MATRIX");
[ 1 0 0 0 1 ]
MATRIX
ok(4, $matrix->is_binary );



