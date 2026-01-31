#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Ask for mode/compiler
# -----------------------------
mode=${1:-}
[[ -z ${mode} ]] && read -p "Enter mode: [R]elease, [D]ebug: " mode

compiler=${2:-}
[[ -z ${compiler} ]] && read -p "Enter compiler: [d]efault, [g]cc, [i]ntel, [c]lang: " compiler

mode=${mode,,}
compiler=${compiler,,}

case "$mode" in
  r|release) mode=Release ;;
  d|debug)   mode=Debug ;;
  *) echo "Invalid mode '$mode'. Use R/Release or D/Debug."; exit 1 ;;
esac

have_cmd() { command -v "$1" >/dev/null 2>&1; }

# Pick first available command from a list
pick_cmd() {
  for c in "$@"; do
    if have_cmd "$c"; then
      echo "$c"
      return 0
    fi
  done
  return 1
}

# -----------------------------
# Decide compilers from PATH (no hardcoding versions)
# -----------------------------
compiler_C=""
compiler_CXX=""

case "$compiler" in
  d|default)
    # Prefer gcc/g++ if present, else clang/clang++, else system cc/c++
    if have_cmd gcc && have_cmd g++; then
      compiler_C=gcc; compiler_CXX=g++
    elif have_cmd clang && have_cmd clang++; then
      compiler_C=clang; compiler_CXX=clang++
    elif have_cmd cc && have_cmd c++; then
      compiler_C=cc; compiler_CXX=c++
    else
      echo "ERROR: No suitable C/C++ compiler found in PATH."
      exit 1
    fi
    ;;
  g|gcc)
    compiler_C=$(pick_cmd gcc cc) || { echo "ERROR: gcc/cc not found in PATH."; exit 1; }
    compiler_CXX=$(pick_cmd g++ c++) || { echo "ERROR: g++/c++ not found in PATH."; exit 1; }
    ;;
  c|clang)
    compiler_C=$(pick_cmd clang cc) || { echo "ERROR: clang/cc not found in PATH."; exit 1; }
    compiler_CXX=$(pick_cmd clang++ c++) || { echo "ERROR: clang++/c++ not found in PATH."; exit 1; }
    ;;
  i|intel)
    # oneAPI first (icx/icpx), then classic (icc/icpc)
    compiler_C=$(pick_cmd icx icc) || { echo "ERROR: Intel C compiler (icx/icc) not found in PATH."; exit 1; }
    compiler_CXX=$(pick_cmd icpx icpc) || { echo "ERROR: Intel C++ compiler (icpx/icpc) not found in PATH."; exit 1; }
    ;;
  *)
    echo "Invalid compiler '$compiler'. Use d/default, g/gcc, i/intel, c/clang."
    exit 1
    ;;
esac

# -----------------------------
# Basic tools checks
# -----------------------------
have_cmd cmake || { echo "ERROR: cmake not found in PATH."; exit 1; }

echo "Build type:   $mode"
echo "C compiler:   $compiler_C  -> $(command -v "$compiler_C")"
echo "C++ compiler: $compiler_CXX -> $(command -v "$compiler_CXX")"
echo "CMake:        $(cmake --version | head -n 1)"

# -----------------------------
# OpenMP retry flags (generic)
# -----------------------------
omp_c_flags=""
omp_cxx_flags=""

case "$compiler_C" in
  gcc|cc)
    # If cc is gcc underneath, -fopenmp is fine; if not, retry may still fail (harmless)
    omp_c_flags="-fopenmp"
    omp_cxx_flags="-fopenmp"
    ;;
  clang)
    # On many Linux setups, clang OpenMP needs libomp; -fopenmp is still the right first flag.
    omp_c_flags="-fopenmp"
    omp_cxx_flags="-fopenmp"
    ;;
  icx|icc)
    omp_c_flags="-qopenmp"
    omp_cxx_flags="-qopenmp"
    ;;
  *)
    # unknown compiler: don't force flags
    omp_c_flags=""
    omp_cxx_flags=""
    ;;
esac

# -----------------------------
# Configure & build
# -----------------------------
rm -rf cmake-build
mkdir -p cmake-build
cd cmake-build

# First configure attempt
set +e
cmake -S .. -B . \
  -DCMAKE_BUILD_TYPE="$mode" \
  -DCMAKE_C_COMPILER="$compiler_C" \
  -DCMAKE_CXX_COMPILER="$compiler_CXX"
rc=$?
set -e

# If OpenMP fails, retry once with explicit OpenMP flags (only if we have flags)
if [[ $rc -ne 0 ]] && grep -q "Could NOT find OpenMP" CMakeFiles/CMakeError.log 2>/dev/null; then
  if [[ -n "${omp_c_flags}" || -n "${omp_cxx_flags}" ]]; then
    echo
    echo "CMake failed to find OpenMP. Retrying with explicit flags:"
    echo "  OpenMP_C_FLAGS=${omp_c_flags}"
    echo "  OpenMP_CXX_FLAGS=${omp_cxx_flags}"
    echo

    rm -f CMakeCache.txt
    rm -rf CMakeFiles

    cmake -S .. -B . \
      -DCMAKE_BUILD_TYPE="$mode" \
      -DCMAKE_C_COMPILER="$compiler_C" \
      -DCMAKE_CXX_COMPILER="$compiler_CXX" \
      -DOpenMP_C_FLAGS="${omp_c_flags}" \
      -DOpenMP_CXX_FLAGS="${omp_cxx_flags}"
  else
    echo "CMake failed and OpenMP flags are unknown for compiler '$compiler_C'."
    exit $rc
  fi
elif [[ $rc -ne 0 ]]; then
  exit $rc
fi

# Parallel build level
if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
  JOBS="$SLURM_CPUS_PER_TASK"
else
  JOBS=4
fi

make -j "$JOBS"
