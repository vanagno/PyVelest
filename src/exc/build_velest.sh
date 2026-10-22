#!/bin/sh
# Compile VELEST with gfortran, in double precision. This is the executable
# PyVelest uses by default (GLOBAL_SETTINGS.velest_exec). The bundled
# src/exc/velest is a single-precision Intel build that needs the Intel
# runtime (libimf.so).
#
# Why double precision: VELEST solves its normal equations with a Cholesky
# decomposition that rejects G as "not positive definite" once a pivot keeps
# less than about 16*n*eps of its diagonal (n = number of unknowns, ~2000 for
# 500 events). In single precision that is ~0.2%, which the weak damping of
# the recipe's stage 1 (velocity 0.1, station delays 0.01) often breaks:
# "error in (damped) matrix G". Double precision leaves the algorithm as it
# is. Its static arrays exceed 2 GB (only the part in use is ever touched),
# hence -mcmodel=medium (x86-64).
#
#   sh src/exc/build_velest.sh            # -> src/exc/velest_gfortran
#   sh src/exc/build_velest.sh my_velest  # -> src/exc/my_velest
set -e
cd "$(dirname "$0")"
out="${1:-velest_gfortran}"
rm -f "$out"
gfortran -O2 -std=legacy -ffixed-line-length-none -fno-automatic \
    -fdefault-real-8 -fdefault-double-8 -mcmodel=medium -Wno-align-commons \
    velest.f -o "$out"
if [ -x "$out" ]; then
    echo "built $(pwd)/$out (double precision)"
else
    echo "build failed: $out was not created" >&2
    exit 1
fi
