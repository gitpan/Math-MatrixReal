# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..4\n"; }
END {print "not ok 1\n" unless $loaded;}
use Math::MatrixReal;
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):

### First, some preparation

my $DEBUG2 = 0;
my $DEBUG = 0;
# Set this one if you want the REAL benchmarking to be done!
my $REALBENCH = 0;
my $bigsize = 40; # Size of big matrix REAL tests (be careful: n^3!)

use Benchmark;

# Useful test functions

sub ok ($$) {
    my($number, $result) = @_ ;
    
    print "ok $number\n"     if $result ;
    print "not ok $number\n" if !$result ;
}

sub ok_matrix ($$$)
{
  my ($no, $M1, $M2) = @_;
  my $tmp = $M1->shadow();
  $tmp->subtract($M1,$M2);
  my $v = $tmp->norm_one();
  ok($no, ($v < 1e-10));
  print " ($no: |Delta| = $v)\n" if $DEBUG;
}

sub ok_eigenvectors ($$$$)
  {
    my ($no, $M, $L, $V) = @_;
    # Now check that all of them correspond to eigenvalue * eigenvector
    my ($rows, $columns) = $M->dim();
    unless ($rows == $columns) {
	ok("$no", 0);
	return;
    }
    # Computes the result of all eigenvectors...
    my $test = $M * $V;
    my $test2 = $V->clone();
    for (my $i = 1; $i <= $columns; $i++)
    {
	my $lambda = $L->element($i,1);
	for (my $j = 1; $j <= $rows; $j++)
	{ # Compute new vector via lambda * x
	    $test2->assign($j, $i, $lambda * $test2->element($j, $i));
	}
      }
    ok_matrix("$no",$test,$test2);
    return;
  }

sub random_matrix ($)
{
    my ($size) = @_;
    my $M = Math::MatrixReal->new($size, $size);
    for (my $i=1; $i<=$size; $i++)
    {
	for (my $j=1; $j<=$size; $j++)
	{
	    $M->assign($i,$j,rand());
	}
    }
    return $M;
}

### We should use the black magic now...

if ($REALBENCH)
{
# test on random bigger matrix
    print "Matrix ".$bigsize."x$bigsize for eigenvalues & eigenvectors computation:\n" if $DEBUG;
# Creates a random matrix
    my $big = random_matrix($bigsize);

#
# Benchmark eigenvalues & eigenvectors computation
#
    $big = $big + ~$big;

    print "Householder reduction...\n" if $DEBUG;
    my ($Tbig, $Qbig);
    my $t = timeit(1, sub { ($Tbig, $Qbig) = $big->householder(); });
    print "Timing of ".$bigsize."x".$bigsize." Householder transformation:\n  ".timestr($t)."\n";
    print "Is Qbig orthogonal?\n" if $DEBUG;
    print "Diagonalization of tridiagonal...\n" if $DEBUG;
    my ($Lbig, $Vbig);
    my $t2 = timeit(1, sub { ($Lbig, $Vbig) = $Tbig->tri_diagonalize($Qbig); });
    print "Timing of ".$bigsize."x".$bigsize." QL-implicit diagonalization:\n  ".timestr($t2)."\n";

# We check the results anyway (just in case...:-)
    ok_eigenvectors(2, $big, $Lbig, $Vbig);

#
# Now test the eigenvalues only computations...
#
    print "Recomputing: Eigenvalues only.\n ".$bigsize."x".$bigsize."\n" if $DEBUG;
    my $altTbig;
    my $t3 = timeit(1, sub { $altTbig = $big->householder_tridiagonal(); });
    print "Timing of ".$bigsize."x".$bigsize." Householder transformation (tridiag. only):\n  ".timestr($t3)."\n";
    my $altLbig;
    my $t4 = timeit(1, sub { $altLbig = $altTbig->tri_eigenvalues(); });
    print "Timing of ".$bigsize."x".$bigsize." QL-implicit eigenvalues computation:\n  ".timestr($t4)."\n";

# We check the results anyway (just in case...:-)
    ok_matrix(3, $altTbig, $Tbig);
    ok_matrix(4, $altLbig, $Lbig);
}
else
{
    # We don't really do the benchmarks on those
    # in fact, this tests does nothing unless told to...
    ok(2, 1 == 1);
    ok(3, 1 == 1);
    ok(4, 1 == 1);
}


