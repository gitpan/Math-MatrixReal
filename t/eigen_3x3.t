# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..14\n"; }
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

### We should use the black magic now...

#
# Trying some matrixes creation extracted from the pod...
#
my $string = "[ 1 2 3 ]\n[ 2 2 -1 ]\n[ 1 1 1 ]\n";
my $matrix33 = Math::MatrixReal->new_from_string($string);
print "$matrix33" if $DEBUG;

#
# Tests eigenvalues & eigenvectors computation
#

#
# First the tridiagonal reduction (Householder) on the 3x3
#
my $symm = $matrix33 + ~$matrix33;
ok(2, not($matrix33->is_symmetric()));
ok(3, $symm->is_symmetric());
print "Matrix 3x3 for eigenvalues & eigenvectors computation:\n$symm"
  if $DEBUG;

print "Householder reduction...\n" if $DEBUG;
my ($T, $Q) = $symm->householder();
print "T=\n$T Q=\n$Q" if $DEBUG2;
print "Is Q orthogonal?\n" if $DEBUG;
print ($Q * ~$Q) if $DEBUG2;
ok_matrix_orthogonal(4, $Q);
ok_matrix(5, $symm, $Q * $T * ~$Q);
print "Diagonalization of tridiagonal...\n" if $DEBUG;
my ($L1, $V1) = $T->tri_diagonalize($Q);
print "eigenvalues L:\n$L1 eigenvectors:\n$V1" if $DEBUG2;
ok_eigenvectors(6,$symm, $L1, $V1);
	
# Get first eigenvector
my $aev1 = $V1->column(1);
my $al1 = $L1->element(1,1);
my $ap1_1 = $symm * $aev1; # A * x
my $ap1_2 = $al1 * $aev1;    # lambda *x
print "Original computation of A*ev1:\n$ap1_1 Scaled eigenvector:\n$ap1_2"
    if $DEBUG2;
ok_matrix(7, $ap1_1, $ap1_2);

print "Direct diagonalization...\n" if $DEBUG;
my ($L12, $V12) = $symm->sym_diagonalize();
print "eigenvalues L:\n$L12 eigenvectors:\n$V12" if $DEBUG2;
ok_eigenvectors(8,$symm, $L12, $V12);
ok_matrix_orthogonal(9, $V12);
# Double check the equality
ok_matrix(10, $L12, $L1);
ok_matrix(11, $V12, $V1);

#
# Now test the eigenvalues only computations...
#
print "Recomputing: Eigenvalues only.\n 3x3\n" if $DEBUG;
my $altT = $symm->householder_tridiagonal();
ok_matrix(12, $altT, $T);
my $altL1 = $altT->tri_eigenvalues();
ok_matrix(13, $altL1, $L1);
my $altL12 = $symm->sym_eigenvalues();
ok_matrix(14, $altL12, $L12);

__END__

# Attempt:
# Obtain the eigenvectors when eigenvalues are known
# using inverse iteration.
# We solve (M - lambda * I) * b(k+1) = b(k)
#  with b(0) a random unit vector and b(k) is
#  normalized at each step.
# This should converge towards the eigenvector,
# but there are problems:
#  - for some value, there is not convergence
#  (the above system is rather singular, so...)
#  - the solution can oscillate between v and -v (?)
# Rodolphe Ortalo, 99/06/14
sub obtain_eigenvector ($$)
{
  my ($M, $eigenvalue) = @_;
  # Form the linear system A - lamda1 * I
  my $inv_it = $M->shadow();
  $inv_it->one();
  $inv_it->multiply_scalar($inv_it, (-1.0 * $eigenvalue));
  $inv_it->add($M, $inv_it);
  print "Linear system matrix:\n $inv_it" if $DEBUG2;
  # Creates a random vector
  my ($rows, $cols) = $inv_it->dim();
  my $b = Math::MatrixReal->new($rows, 1);
  for (my $i = 1; $i <= $rows; $i++)
    {
      $b->assign($i, 1, rand());
    }
  # Normalize it
  my $l = $b->length();
  $b->multiply_scalar($b, (1.0 / $l));
  # Now do LR decomposition for linear system
  my $inv_it_LR = $inv_it->decompose_LR();
  # Check iterations
  my $iter = 0;
  my $delta;
  do {
    my ($dim, $b_base, $base) = $inv_it_LR->solve_LR($b);
    # Normalize
    my $l = $b_base->length();
    $b_base->multiply_scalar($b_base, (-1.0 / $l));
#    print "b_base=\n$b_base";
    $b->subtract($b_base,$b);
    $delta = $b->norm_one();
    print "delta=$delta\n";
    $b = $b_base;
  } while (($delta >= 1e-10) && ($iter++ <= 10));
  return $b;
}
#
# Now, try to find one eigenvector again...
# (Using Steffen's functions...:-)
#
my $ev = obtain_eigenvector($symm, $al1);
print "Real ev:\n $aev1 Found ev:\n $ev" if $DEBUG;
ok_matrix(15, $ev, $aev1);
ok_matrix(16, obtain_eigenvector($symm, $L1->element(2,1)),
	  $V1->column(2));
ok_matrix(17, obtain_eigenvector($symm, $L1->element(3,1)),
	  $V1->column(3));
