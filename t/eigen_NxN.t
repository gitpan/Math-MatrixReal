# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..12\n"; }
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
my $bigsize = 30; # Size of big matrix tests (be careful: n^3!)

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

sub ok_matrix_orthogonal ($$)
{
  my ($no, $M) = @_;
  my $tmp = $M->shadow();
  $tmp->one();
  my $transp = $M->shadow();
  $transp->transpose($M);
  $tmp->subtract($M->multiply($transp), $tmp);
  my $v = $tmp->norm_one();
  ok($no, ($v < 1e-10));
  print " ($no: |M * ~M - I| = $v)\n" if $DEBUG;
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

# test on random bigger matrix
print "Matrix ".$bigsize."x$bigsize for eigenvalues & eigenvectors computation:\n" if $DEBUG;
# Creates a random matrix
my $big = random_matrix($bigsize);

#
# Tests eigenvalues & eigenvectors computation
#
$big = $big + ~$big; # Symmetric matrix

print "Householder reduction...\n" if $DEBUG;
my ($Tbig, $Qbig) = $big->householder();
print "Is Qbig orthogonal?\n" if $DEBUG;
ok_matrix_orthogonal(2, $Qbig);
ok_matrix(3, $big, $Qbig * $Tbig * ~$Qbig);
print "Diagonalization of tridiagonal...\n" if $DEBUG;
my ($Lbig, $Vbig) = $Tbig->tri_diagonalize($Qbig);
ok_eigenvectors(4, $big, $Lbig, $Vbig);
ok_matrix_orthogonal(5, $Vbig);

print "Direct diagonalization...\n" if $DEBUG;
my ($Lbig_2, $Vbig_2) = $big->sym_diagonalize();
print "eigenvalues L:\n$Lbig_2 eigenvectors:\n$Vbig_2" if $DEBUG2;
ok_eigenvectors(6,$big, $Lbig_2, $Vbig_2);
ok_matrix_orthogonal(7, $Vbig_2);
# Double check the equality
ok_matrix(8, $Lbig_2, $Lbig);
ok_matrix(9, $Vbig_2, $Vbig);

#
# Now test the eigenvalues only computations...
#
print "Recomputing: Eigenvalues only.\n ".$bigsize."x".$bigsize."\n" if $DEBUG;
my $altTbig = $big->householder_tridiagonal();
ok_matrix(10, $altTbig, $Tbig);
my $altLbig = $altTbig->tri_eigenvalues();
ok_matrix(11, $altLbig, $Lbig);
my $altLbig_2 = $big->sym_eigenvalues();
ok_matrix(12, $altLbig_2, $Lbig_2);


