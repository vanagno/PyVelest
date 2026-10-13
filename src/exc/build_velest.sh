#!/bin/sh
# Compile VELEST with gfortran. The bundled src/exc/velest was built with the
# Intel compiler and needs its runtime (libimf.so); use this when that is
# missing, then set GLOBAL_SETTINGS.velest_exec in configs/user.yaml.
#
#   sh src/exc/build_velest.sh            # -> src/exc/velest_gfortran
#   sh src/exc/build_velest.sh my_velest  # -> src/exc/my_velest
set -e
cd "$(dirname "$0")"
out="${1:-velest_gfortran}"
rm -f "$out"
gfortran -O2 -std=legacy -ffixed-line-length-none -fno-automatic \
    velest.f -o "$out"
if [ -x "$out" ]; then
    echo "built $(pwd)/$out"
else
    echo "build failed: $out was not created" >&2
    exit 1
fi
