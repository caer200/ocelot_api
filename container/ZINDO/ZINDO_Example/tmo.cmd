#!/bin/sh

ExeDir=/home/rajesh/VIrtual_Machine/Research/Simulations/ZINDO/ZINDO_Example

$ExeDir/zindo < A.inp > A.out
mv mo.out mo_A.out
$ExeDir/zindo < D.inp > D.out
mv mo.out mo_D.out
$ExeDir/zindo-ct < dimer.inp > dimer.out

