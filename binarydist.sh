#!/bin/sh -e
#
# Build a redistributable binary tarball of bpopsim.
#
# Produces bpopsim-<version>-<platform>-<arch>.tar.gz containing bin/bpopsim and
# share/man/man1/bpopsim.1, with GSL linked statically so the binary has no
# dependency on the build machine's library paths.
#
# Usage (from the repository root):
#   conda run -p $PWD/env ./binarydist.sh [--rebuild-prereqs]
#
# macOS: builds a universal (arm64 + x86_64) binary with the system Apple clang,
#   following the same approach as breseq's binarydist.sh. GSL is built from
#   source as a static universal archive into binarydist_build/ (once; reused
#   until the pinned version in environment.yml changes or --rebuild-prereqs is
#   given). The conda env's GSL and compiler are not used: they are single
#   architecture and would leave an rpath dependency on the env's libc++.
# Linux: builds for the host architecture with the env's compiler and links the
#   env's static libgsl.a.

REBUILD_PREREQS=0
if [ "$1" = "--rebuild-prereqs" ]; then
  REBUILD_PREREQS=1
fi

VERSION=$(sed -n 's/^[[:space:]]*VERSION[[:space:]][[:space:]]*\([0-9][0-9.]*\).*/\1/p' CMakeLists.txt | head -1)
if [ -z "$VERSION" ]; then
  echo "Could not read project VERSION from CMakeLists.txt" >&2
  exit 1
fi

GSL_VERSION=$(grep -E '^[[:space:]]*-[[:space:]]*gsl=' environment.yml | sed -E 's/.*=([0-9.]+).*/\1/')
if [ -z "$GSL_VERSION" ]; then
  echo "Could not read the gsl= pin from environment.yml" >&2
  exit 1
fi

PLATFORM=$(uname -s)
ARCH=$(uname -m)
CMAKE_EXTRA=""
PREREQDIR="${PWD}/binarydist_build"

if [ "$PLATFORM" = "Darwin" ]; then
  MACOSX_MIN=10.13
  ARCHS="arm64 x86_64"
  PLATFORM="MacOSX-${MACOSX_MIN}+"
  ARCH="universal"

  # System Apple clang can emit true universal binaries via multiple -arch flags.
  # Conda's compilers are single-architecture cross-compilers.
  CC=/usr/bin/clang
  CXX=/usr/bin/clang++

  # Discard any flags inherited from an active conda environment; they would
  # pull in the env's single-architecture headers and libraries.
  CFLAGS=""; CXXFLAGS=""; CPPFLAGS=""; LDFLAGS=""
  export CFLAGS CXXFLAGS CPPFLAGS LDFLAGS

  GSLDIR="${PREREQDIR}/gsl"
  if [ "$REBUILD_PREREQS" -eq 1 ] || [ ! -f "${GSLDIR}/.version" ] \
     || [ "$(cat "${GSLDIR}/.version")" != "$GSL_VERSION" ]; then
    echo "Building static universal GSL ${GSL_VERSION} into ${GSLDIR}"
    rm -rf "${GSLDIR}" "${GSLDIR}"-*
    mkdir -p "${PREREQDIR}/src"
    TARBALL="${PREREQDIR}/src/gsl-${GSL_VERSION}.tar.gz"
    if [ ! -f "$TARBALL" ]; then
      curl -L -o "$TARBALL" "https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz"
    fi

    for a in $ARCHS; do
      SRC="${PREREQDIR}/src/gsl-${GSL_VERSION}-${a}"
      rm -rf "$SRC"
      mkdir -p "$SRC"
      tar -xzf "$TARBALL" -C "$SRC" --strip-components 1
      HOSTFLAG=""
      if [ "$a" != "$(uname -m)" ]; then
        HOSTFLAG="--host=${a}-apple-darwin"
      fi
      (
        cd "$SRC"
        ./configure $HOSTFLAG --prefix="${GSLDIR}-${a}" --disable-shared --enable-static \
          CC="$CC" CFLAGS="-arch ${a} -mmacosx-version-min=${MACOSX_MIN} -O2" \
          > configure.log 2>&1
        make -j > make.log 2>&1
        make install > install.log 2>&1
      )
    done

    # Merge the per-architecture static archives into one universal tree.
    FIRST_ARCH=$(echo $ARCHS | cut -d' ' -f1)
    mkdir -p "${GSLDIR}/lib"
    cp -R "${GSLDIR}-${FIRST_ARCH}/include" "${GSLDIR}/include"
    cp -R "${GSLDIR}-${FIRST_ARCH}/bin" "${GSLDIR}/bin"
    for lib in libgsl.a libgslcblas.a; do
      INPUTS=""
      for a in $ARCHS; do INPUTS="$INPUTS ${GSLDIR}-${a}/lib/${lib}"; done
      lipo -create $INPUTS -output "${GSLDIR}/lib/${lib}"
    done
    echo "$GSL_VERSION" > "${GSLDIR}/.version"
  fi

  CMAKE_EXTRA="-DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} \
    -DCMAKE_OSX_ARCHITECTURES=arm64;x86_64 \
    -DCMAKE_OSX_DEPLOYMENT_TARGET=${MACOSX_MIN} \
    -DGSL_ROOT_DIR=${GSLDIR}"
fi

NAME="bpopsim-${VERSION}-${PLATFORM}-${ARCH}"
BUILD_DIR="build-dist"
STAGE_DIR="${PWD}/${NAME}"

rm -rf "${BUILD_DIR}" "${STAGE_DIR}" "${NAME}.tar.gz"

cmake -S . -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBPOPSIM_STATIC_GSL=ON \
  -DCMAKE_INSTALL_PREFIX="${STAGE_DIR}" \
  $CMAKE_EXTRA
cmake --build "${BUILD_DIR}" --parallel
ctest --test-dir "${BUILD_DIR}" --output-on-failure
cmake --install "${BUILD_DIR}"

cp COPYING "${STAGE_DIR}/"

tar -czf "${NAME}.tar.gz" "${NAME}"
rm -rf "${STAGE_DIR}"

echo "Built ${NAME}.tar.gz"
