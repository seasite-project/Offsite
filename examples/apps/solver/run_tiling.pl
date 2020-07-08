#!/usr/bin/perl
################################################################################

use strict;
use warnings;
use File::Copy;
use POSIX qw(strftime);

################################################################################

sub create_filename
{
  return "table." . strftime("%Y-%m-%d.%H:%M", localtime) . ".pmd";
}

################################################################################

my $openmp = 1;

my $cc = "icc";
my $cflags =
  "-Wall -Wno-unknown-pragmas -Wno-unused-function -O3 -D_GNU_SOURCE -DBENCHMARK";
$cflags .= " -fopenmp" if ($openmp);

my $ldlibs = "-lm";

if ($cc eq 'icc')
{
    $cflags .= " -xHost -inline-max-size=10000 -inline-max-total-size=10000 -fno-alias";
}
else
{
    $cflags .= " -march=native";
}

################################################################################

my $table  = create_filename();
my $header = "header.pmd";

copy($header, $table);

################################################################################

my $variant_folder = "variants_inverter_radauIIIA7";

my @impls       = qw(
A_RHS_jl_LC_ijl_Approx_ij_Update_j
A_RHS_jl_LC_ijl_Approx_ji_Update_j
A_RHS_jl_LC_ilj_Approx_ij_Update_j
A_RHS_jl_LC_ilj_Approx_ji_Update_j
A_RHS_jl_LC_jil_Approx_ij_Update_j
A_RHS_jl_LC_jil_Approx_ji_Update_j
A_RHS_jl_LC_jli_Approx_ij_Update_j
A_RHS_jl_LC_jli_Approx_ji_Update_j
A_RHS_jl_LC_lij_Approx_ij_Update_j
A_RHS_jl_LC_lij_Approx_ji_Update_j
A_RHS_jl_LC_lji_Approx_ij_Update_j
A_RHS_jl_LC_lji_Approx_ji_Update_j
A_RHS_lj_LC_ijl_Approx_ij_Update_j
A_RHS_lj_LC_ijl_Approx_ji_Update_j
A_RHS_lj_LC_ilj_Approx_ij_Update_j
A_RHS_lj_LC_ilj_Approx_ji_Update_j
A_RHS_lj_LC_jil_Approx_ij_Update_j
A_RHS_lj_LC_jil_Approx_ji_Update_j
A_RHS_lj_LC_jli_Approx_ij_Update_j
A_RHS_lj_LC_jli_Approx_ji_Update_j
A_RHS_lj_LC_lij_Approx_ij_Update_j
A_RHS_lj_LC_lij_Approx_ji_Update_j
A_RHS_lj_LC_lji_Approx_ij_Update_j
A_RHS_lj_LC_lji_Approx_ji_Update_j
B_RHS_jl_LC_jli_ApproxUpdate_ji
B_RHS_jl_LC_lij_ApproxUpdate_ij
B_RHS_jl_LC_lij_ApproxUpdate_ji
B_RHS_jl_LC_lji_ApproxUpdate_ij
B_RHS_jl_LC_lji_ApproxUpdate_ji
B_RHS_jl_LC_ijl_ApproxUpdate_ij
B_RHS_jl_LC_ijl_ApproxUpdate_ji
B_RHS_jl_LC_ilj_ApproxUpdate_ij
B_RHS_jl_LC_ilj_ApproxUpdate_ji
B_RHS_jl_LC_jil_ApproxUpdate_ij
B_RHS_jl_LC_jil_ApproxUpdate_ji
B_RHS_jl_LC_jli_ApproxUpdate_ij
B_RHS_lj_LC_jli_ApproxUpdate_ji
B_RHS_lj_LC_lij_ApproxUpdate_ij
B_RHS_lj_LC_lij_ApproxUpdate_ji
B_RHS_lj_LC_lji_ApproxUpdate_ij
B_RHS_lj_LC_lji_ApproxUpdate_ji
B_RHS_lj_LC_ijl_ApproxUpdate_ij
B_RHS_lj_LC_ijl_ApproxUpdate_ji
B_RHS_lj_LC_ilj_ApproxUpdate_ij
B_RHS_lj_LC_ilj_ApproxUpdate_ji
B_RHS_lj_LC_jil_ApproxUpdate_ij
B_RHS_lj_LC_jil_ApproxUpdate_ji
B_RHS_lj_LC_jli_ApproxUpdate_ij
C_RHS_LC_ijl_RHS_jl_Approx_ij_Update_j
C_RHS_LC_jli_RHS_jl_Approx_ij_Update_j
C_RHS_LC_ijl_RHS_jl_Approx_ji_Update_j
C_RHS_LC_jli_RHS_jl_Approx_ji_Update_j
C_RHS_LC_ijl_RHS_lj_Approx_ij_Update_j
C_RHS_LC_jli_RHS_lj_Approx_ij_Update_j
C_RHS_LC_ijl_RHS_lj_Approx_ji_Update_j
C_RHS_LC_jli_RHS_lj_Approx_ji_Update_j
C_RHS_LC_ilj_RHS_jl_Approx_ij_Update_j
C_RHS_LC_ilj_RHS_jl_Approx_ji_Update_j
C_RHS_LC_ilj_RHS_lj_Approx_ij_Update_j
C_RHS_LC_ilj_RHS_lj_Approx_ji_Update_j
C_RHS_LC_jil_RHS_jl_Approx_ij_Update_j
C_RHS_LC_jil_RHS_jl_Approx_ji_Update_j
C_RHS_LC_jil_RHS_lj_Approx_ij_Update_j
C_RHS_LC_jil_RHS_lj_Approx_ji_Update_j
D_RHS_LC_ijl_RHS_jl_ApproxUpdate_ij
D_RHS_LC_jli_RHS_jl_ApproxUpdate_ij
D_RHS_LC_ijl_RHS_jl_ApproxUpdate_ji
D_RHS_LC_jli_RHS_jl_ApproxUpdate_ji
D_RHS_LC_ijl_RHS_lj_ApproxUpdate_ij
D_RHS_LC_jli_RHS_lj_ApproxUpdate_ij
D_RHS_LC_ijl_RHS_lj_ApproxUpdate_ji
D_RHS_LC_jli_RHS_lj_ApproxUpdate_ji
D_RHS_LC_ilj_RHS_jl_ApproxUpdate_ij
D_RHS_LC_ilj_RHS_jl_ApproxUpdate_ji
D_RHS_LC_ilj_RHS_lj_ApproxUpdate_ij
D_RHS_LC_ilj_RHS_lj_ApproxUpdate_ji
D_RHS_LC_jil_RHS_jl_ApproxUpdate_ij
D_RHS_LC_jil_RHS_jl_ApproxUpdate_ji
D_RHS_LC_jil_RHS_lj_ApproxUpdate_ij
D_RHS_LC_jil_RHS_lj_ApproxUpdate_ji
E_RHS_LC_ijl_RHS_Approx_ij_Update_j
E_RHS_LC_jli_RHS_Approx_ij_Update_j
E_RHS_LC_ijl_RHS_Approx_ji_Update_j
E_RHS_LC_jli_RHS_Approx_ji_Update_j
E_RHS_LC_ilj_RHS_Approx_ij_Update_j
E_RHS_LC_ilj_RHS_Approx_ji_Update_j
E_RHS_LC_jil_RHS_Approx_ij_Update_j
E_RHS_LC_jil_RHS_Approx_ji_Update_j
F_RHS_LC_ijl_RHS_Approx_Update_ij
F_RHS_LC_jli_RHS_Approx_Update_ij
F_RHS_LC_ijl_RHS_Approx_Update_ji
F_RHS_LC_jli_RHS_Approx_Update_ji
F_RHS_LC_ilj_RHS_Approx_Update_ij
F_RHS_LC_ilj_RHS_Approx_Update_ji
F_RHS_LC_jil_RHS_Approx_Update_ij
F_RHS_LC_jil_RHS_Approx_Update_ji
G_RHS_jl_LC_jli_RHS_Approx_ij_Update_j
G_RHS_jl_LC_jli_RHS_Approx_ji_Update_j
G_RHS_jl_LC_lij_RHS_Approx_ij_Update_j
G_RHS_jl_LC_lij_RHS_Approx_ji_Update_j
G_RHS_jl_LC_lji_RHS_Approx_ij_Update_j
G_RHS_jl_LC_lji_RHS_Approx_ji_Update_j
G_RHS_jl_LC_ijl_RHS_Approx_ij_Update_j
G_RHS_jl_LC_ijl_RHS_Approx_ji_Update_j
G_RHS_jl_LC_ilj_RHS_Approx_ij_Update_j
G_RHS_jl_LC_ilj_RHS_Approx_ji_Update_j
G_RHS_jl_LC_jil_RHS_Approx_ij_Update_j
G_RHS_jl_LC_jil_RHS_Approx_ji_Update_j
G_RHS_lj_LC_jli_RHS_Approx_ij_Update_j
G_RHS_lj_LC_jli_RHS_Approx_ji_Update_j
G_RHS_lj_LC_lij_RHS_Approx_ij_Update_j
G_RHS_lj_LC_lij_RHS_Approx_ji_Update_j
G_RHS_lj_LC_lji_RHS_Approx_ij_Update_j
G_RHS_lj_LC_lji_RHS_Approx_ji_Update_j
G_RHS_lj_LC_ijl_RHS_Approx_ij_Update_j
G_RHS_lj_LC_ijl_RHS_Approx_ji_Update_j
G_RHS_lj_LC_ilj_RHS_Approx_ij_Update_j
G_RHS_lj_LC_ilj_RHS_Approx_ji_Update_j
G_RHS_lj_LC_jil_RHS_Approx_ij_Update_j
G_RHS_lj_LC_jil_RHS_Approx_ji_Update_j
H_RHS_jl_LC_ijl_RHS_Approx_Update_ij
H_RHS_jl_LC_ijl_RHS_Approx_Update_ji
H_RHS_jl_LC_ilj_RHS_Approx_Update_ij
H_RHS_jl_LC_ilj_RHS_Approx_Update_ji
H_RHS_jl_LC_jil_RHS_Approx_Update_ij
H_RHS_jl_LC_jil_RHS_Approx_Update_ji
H_RHS_jl_LC_jli_RHS_Approx_Update_ij
H_RHS_jl_LC_jli_RHS_Approx_Update_ji
H_RHS_jl_LC_lij_RHS_Approx_Update_ij
H_RHS_jl_LC_lij_RHS_Approx_Update_ji
H_RHS_jl_LC_lji_RHS_Approx_Update_ij
H_RHS_jl_LC_lji_RHS_Approx_Update_ji
H_RHS_lj_LC_ijl_RHS_Approx_Update_ij
H_RHS_lj_LC_ijl_RHS_Approx_Update_ji
H_RHS_lj_LC_ilj_RHS_Approx_Update_ij
H_RHS_lj_LC_ilj_RHS_Approx_Update_ji
H_RHS_lj_LC_jil_RHS_Approx_Update_ij
H_RHS_lj_LC_jil_RHS_Approx_Update_ji
H_RHS_lj_LC_jli_RHS_Approx_Update_ij
H_RHS_lj_LC_jli_RHS_Approx_Update_ji
H_RHS_lj_LC_lij_RHS_Approx_Update_ij
H_RHS_lj_LC_lij_RHS_Approx_Update_ji
H_RHS_lj_LC_lji_RHS_Approx_Update_ij
H_RHS_lj_LC_lji_RHS_Approx_Update_ji
A_RHS_jl_LC_ijl_Approx_ij_Update_j_tiled
A_RHS_jl_LC_ijl_Approx_ji_Update_j_tiled
A_RHS_jl_LC_ilj_Approx_ij_Update_j_tiled
A_RHS_jl_LC_ilj_Approx_ji_Update_j_tiled
A_RHS_jl_LC_jil_Approx_ij_Update_j_tiled
A_RHS_jl_LC_jil_Approx_ji_Update_j_tiled
A_RHS_jl_LC_jli_Approx_ij_Update_j_tiled
A_RHS_jl_LC_jli_Approx_ji_Update_j_tiled
A_RHS_jl_LC_lij_Approx_ij_Update_j_tiled
A_RHS_jl_LC_lij_Approx_ji_Update_j_tiled
A_RHS_jl_LC_lji_Approx_ij_Update_j_tiled
A_RHS_jl_LC_lji_Approx_ji_Update_j_tiled
A_RHS_lj_LC_ijl_Approx_ij_Update_j_tiled
A_RHS_lj_LC_ijl_Approx_ji_Update_j_tiled
A_RHS_lj_LC_ilj_Approx_ij_Update_j_tiled
A_RHS_lj_LC_ilj_Approx_ji_Update_j_tiled
A_RHS_lj_LC_jil_Approx_ij_Update_j_tiled
A_RHS_lj_LC_jil_Approx_ji_Update_j_tiled
A_RHS_lj_LC_jli_Approx_ij_Update_j_tiled
A_RHS_lj_LC_jli_Approx_ji_Update_j_tiled
A_RHS_lj_LC_lij_Approx_ij_Update_j_tiled
A_RHS_lj_LC_lij_Approx_ji_Update_j_tiled
A_RHS_lj_LC_lji_Approx_ij_Update_j_tiled
A_RHS_lj_LC_lji_Approx_ji_Update_j_tiled
B_RHS_jl_LC_jli_ApproxUpdate_ji_tiled
B_RHS_jl_LC_lij_ApproxUpdate_ij_tiled
B_RHS_jl_LC_lij_ApproxUpdate_ji_tiled
B_RHS_jl_LC_lji_ApproxUpdate_ij_tiled
B_RHS_jl_LC_lji_ApproxUpdate_ji_tiled
B_RHS_jl_LC_ijl_ApproxUpdate_ij_tiled
B_RHS_jl_LC_ijl_ApproxUpdate_ji_tiled
B_RHS_jl_LC_ilj_ApproxUpdate_ij_tiled
B_RHS_jl_LC_ilj_ApproxUpdate_ji_tiled
B_RHS_jl_LC_jil_ApproxUpdate_ij_tiled
B_RHS_jl_LC_jil_ApproxUpdate_ji_tiled
B_RHS_jl_LC_jli_ApproxUpdate_ij_tiled
B_RHS_lj_LC_jli_ApproxUpdate_ji_tiled
B_RHS_lj_LC_lij_ApproxUpdate_ij_tiled
B_RHS_lj_LC_lij_ApproxUpdate_ji_tiled
B_RHS_lj_LC_lji_ApproxUpdate_ij_tiled
B_RHS_lj_LC_lji_ApproxUpdate_ji_tiled
B_RHS_lj_LC_ijl_ApproxUpdate_ij_tiled
B_RHS_lj_LC_ijl_ApproxUpdate_ji_tiled
B_RHS_lj_LC_ilj_ApproxUpdate_ij_tiled
B_RHS_lj_LC_ilj_ApproxUpdate_ji_tiled
B_RHS_lj_LC_jil_ApproxUpdate_ij_tiled
B_RHS_lj_LC_jil_ApproxUpdate_ji_tiled
B_RHS_lj_LC_jli_ApproxUpdate_ij_tiled
C_RHS_LC_ijl_RHS_jl_Approx_ij_Update_j
C_RHS_LC_jli_RHS_jl_Approx_ij_Update_j_tiled
C_RHS_LC_ijl_RHS_jl_Approx_ji_Update_j_tiled
C_RHS_LC_jli_RHS_jl_Approx_ji_Update_j_tiled
C_RHS_LC_ijl_RHS_lj_Approx_ij_Update_j_tiled
C_RHS_LC_jli_RHS_lj_Approx_ij_Update_j_tiled
C_RHS_LC_ijl_RHS_lj_Approx_ji_Update_j_tiled
C_RHS_LC_jli_RHS_lj_Approx_ji_Update_j_tiled
C_RHS_LC_ilj_RHS_jl_Approx_ij_Update_j_tiled
C_RHS_LC_ilj_RHS_jl_Approx_ji_Update_j_tiled
C_RHS_LC_ilj_RHS_lj_Approx_ij_Update_j_tiled
C_RHS_LC_ilj_RHS_lj_Approx_ji_Update_j_tiled
C_RHS_LC_jil_RHS_jl_Approx_ij_Update_j_tiled
C_RHS_LC_jil_RHS_jl_Approx_ji_Update_j_tiled
C_RHS_LC_jil_RHS_lj_Approx_ij_Update_j_tiled
C_RHS_LC_jil_RHS_lj_Approx_ji_Update_j_tiled
D_RHS_LC_ijl_RHS_jl_ApproxUpdate_ij_tiled
D_RHS_LC_jli_RHS_jl_ApproxUpdate_ij_tiled
D_RHS_LC_ijl_RHS_jl_ApproxUpdate_ji_tiled
D_RHS_LC_jli_RHS_jl_ApproxUpdate_ji_tiled
D_RHS_LC_ijl_RHS_lj_ApproxUpdate_ij_tiled
D_RHS_LC_jli_RHS_lj_ApproxUpdate_ij_tiled
D_RHS_LC_ijl_RHS_lj_ApproxUpdate_ji_tiled
D_RHS_LC_jli_RHS_lj_ApproxUpdate_ji_tiled
D_RHS_LC_ilj_RHS_jl_ApproxUpdate_ij_tiled
D_RHS_LC_ilj_RHS_jl_ApproxUpdate_ji_tiled
D_RHS_LC_ilj_RHS_lj_ApproxUpdate_ij_tiled
D_RHS_LC_ilj_RHS_lj_ApproxUpdate_ji_tiled
D_RHS_LC_jil_RHS_jl_ApproxUpdate_ij_tiled
D_RHS_LC_jil_RHS_jl_ApproxUpdate_ji_tiled
D_RHS_LC_jil_RHS_lj_ApproxUpdate_ij_tiled
D_RHS_LC_jil_RHS_lj_ApproxUpdate_ji_tiled
E_RHS_LC_ijl_RHS_Approx_ij_Update_j_tiled
E_RHS_LC_jli_RHS_Approx_ij_Update_j_tiled
E_RHS_LC_ijl_RHS_Approx_ji_Update_j_tiled
E_RHS_LC_jli_RHS_Approx_ji_Update_j_tiled
E_RHS_LC_ilj_RHS_Approx_ij_Update_j_tiled
E_RHS_LC_ilj_RHS_Approx_ji_Update_j_tiled
E_RHS_LC_jil_RHS_Approx_ij_Update_j_tiled
E_RHS_LC_jil_RHS_Approx_ji_Update_j_tiled
F_RHS_LC_ijl_RHS_Approx_Update_ij_tiled
F_RHS_LC_jli_RHS_Approx_Update_ij_tiled
F_RHS_LC_ijl_RHS_Approx_Update_ji_tiled
F_RHS_LC_jli_RHS_Approx_Update_ji_tiled
F_RHS_LC_ilj_RHS_Approx_Update_ij_tiled
F_RHS_LC_ilj_RHS_Approx_Update_ji_tiled
F_RHS_LC_jil_RHS_Approx_Update_ij_tiled
F_RHS_LC_jil_RHS_Approx_Update_ji_tiled
G_RHS_jl_LC_jli_RHS_Approx_ij_Update_j_tiled
G_RHS_jl_LC_jli_RHS_Approx_ji_Update_j_tiled
G_RHS_jl_LC_lij_RHS_Approx_ij_Update_j_tiled
G_RHS_jl_LC_lij_RHS_Approx_ji_Update_j_tiled
G_RHS_jl_LC_lji_RHS_Approx_ij_Update_j_tiled
G_RHS_jl_LC_lji_RHS_Approx_ji_Update_j_tiled
G_RHS_jl_LC_ijl_RHS_Approx_ij_Update_j_tiled
G_RHS_jl_LC_ijl_RHS_Approx_ji_Update_j_tiled
G_RHS_jl_LC_ilj_RHS_Approx_ij_Update_j_tiled
G_RHS_jl_LC_ilj_RHS_Approx_ji_Update_j_tiled
G_RHS_jl_LC_jil_RHS_Approx_ij_Update_j_tiled
G_RHS_jl_LC_jil_RHS_Approx_ji_Update_j_tiled
G_RHS_lj_LC_jli_RHS_Approx_ij_Update_j_tiled
G_RHS_lj_LC_jli_RHS_Approx_ji_Update_j_tiled
G_RHS_lj_LC_lij_RHS_Approx_ij_Update_j_tiled
G_RHS_lj_LC_lij_RHS_Approx_ji_Update_j_tiled
G_RHS_lj_LC_lji_RHS_Approx_ij_Update_j_tiled
G_RHS_lj_LC_lji_RHS_Approx_ji_Update_j_tiled
G_RHS_lj_LC_ijl_RHS_Approx_ij_Update_j_tiled
G_RHS_lj_LC_ijl_RHS_Approx_ji_Update_j_tiled
G_RHS_lj_LC_ilj_RHS_Approx_ij_Update_j_tiled
G_RHS_lj_LC_ilj_RHS_Approx_ji_Update_j_tiled
G_RHS_lj_LC_jil_RHS_Approx_ij_Update_j_tiled
G_RHS_lj_LC_jil_RHS_Approx_ji_Update_j_tiled
H_RHS_jl_LC_ijl_RHS_Approx_Update_ij_tiled
H_RHS_jl_LC_ijl_RHS_Approx_Update_ji_tiled
H_RHS_jl_LC_ilj_RHS_Approx_Update_ij_tiled
H_RHS_jl_LC_ilj_RHS_Approx_Update_ji_tiled
H_RHS_jl_LC_jil_RHS_Approx_Update_ij_tiled
H_RHS_jl_LC_jil_RHS_Approx_Update_ji_tiled
H_RHS_jl_LC_jli_RHS_Approx_Update_ij_tiled
H_RHS_jl_LC_jli_RHS_Approx_Update_ji_tiled
H_RHS_jl_LC_lij_RHS_Approx_Update_ij_tiled
H_RHS_jl_LC_lij_RHS_Approx_Update_ji_tiled
H_RHS_jl_LC_lji_RHS_Approx_Update_ij_tiled
H_RHS_jl_LC_lji_RHS_Approx_Update_ji_tiled
H_RHS_lj_LC_ijl_RHS_Approx_Update_ij_tiled
H_RHS_lj_LC_ijl_RHS_Approx_Update_ji_tiled
H_RHS_lj_LC_ilj_RHS_Approx_Update_ij_tiled
H_RHS_lj_LC_ilj_RHS_Approx_Update_ji_tiled
H_RHS_lj_LC_jil_RHS_Approx_Update_ij_tiled
H_RHS_lj_LC_jil_RHS_Approx_Update_ji_tiled
H_RHS_lj_LC_jli_RHS_Approx_Update_ij_tiled
H_RHS_lj_LC_jli_RHS_Approx_Update_ji_tiled
H_RHS_lj_LC_lij_RHS_Approx_Update_ij_tiled
H_RHS_lj_LC_lij_RHS_Approx_Update_ji_tiled
H_RHS_lj_LC_lji_RHS_Approx_Update_ij_tiled
H_RHS_lj_LC_lji_RHS_Approx_Update_ji_tiled
);

my $block_factors="-DB_RHS_jl=320 -DB_LC_ijl=320 -DB_Approx_ij=320 -DB_Update_j=320 -DB_Approx_ji=320 -DB_LC_ilj=320 -DB_LC_jil=320 -DB_LC_jli=320 -DB_LC_lij=320 -DB_LC_lji=320 -DB_RHS_lj=320 -DB_ApproxUpdate_ji=320 -DB_ApproxUpdate_ij=320 -DB_RHS_LC_jli=320 -DB_RHS_LC_jil=320 -DB_RHS_LC_ijl=320 -DB_RHS_LC_ilj=320 -DB_RHS_Approx_Update_ji=320 -DB_RHS_Approx_Update_ij=320 -DB_RHS_ApproxUpdate_ij=320 -DB_RHS_ApproxUpdate_ji=320 -DB_RHS_Approx_ij=320 -DB_RHS_Approx_ji=320";

my @threads     = ($openmp ? (8) : (1));

################################################################################

sub run($$$$)
{
    my ($i, $t, $N, $n) = @_;

    return if ($n % $t != 0);

    my $tpcs          = 2.5e-7;
    my $h             = 0.62;
    my $tps           = $tpcs * $n;
    my $steps_per_10s = 10.0 / $tps;
    if ($steps_per_10s < 5) { $steps_per_10s = 5; }
    my $te = $steps_per_10s * $h;

    print "* $i, n=$n, $t threads\n";

    print "  - Compiling ...\n";
    my $cmd =
"$cc $cflags $block_factors -DVARIANT_HEADER=\\\"$variant_folder/$i.h\\\" -I. -Dg=$N -Dn=$n -Dt0=0 -Dte=$te -Dh0=$h -DTHREADS=$t -DVARIANT_NAME=\\\"$i\\\" -o ./bench driver.c $ldlibs";

    system($cmd);
    print "  - Running benchmark ... \n";
    $cmd = "./bench >>$table";
    $cmd = "export KMP_AFFINITY=granularity=fine,compact,1,0; " . $cmd if ($openmp);
    system($cmd);
}

################################################################################

for my $t (@threads)
{
    my $N_increment = 20;
    if ($t == 1) { $N_increment = 100; }
    #for (my $N = 20; $N <= 20; $N += $N_increment)
    for (my $N = 6000; $N <= 6000; $N += $N_increment)
    {
        my $n = $N * $N;
        if ($N >= 1000)
        { 
            $N_increment = 100; 
            if ($t == 1) { $N_increment = 200; }  
        }
        if ($N >= 2000) { $N_increment = 200; }
        for my $i (@impls)
        {
            run($i, $t, $N, $n);
        }
    }
}
