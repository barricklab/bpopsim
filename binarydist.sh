VERSIONSTRING=`perl -ne 's/AC_INIT\(\[(.+?)\].+?\[(.+?)\].+/\1-\2/ && print' configure.ac`

./bootstrap.sh
BINARYPLATFORM=`uname`
BINARYARCH=`arch`
ARCHFLAGS=""
BINARYNAME=${BINARYPLATFORM}-${BINARYARCH}
if [ "$BINARYPLATFORM" == "Darwin" ]; then
	BINARYARCH="universal"
	ARCHFLAGS="-arch i386 -arch x86_64 -mmacosx-version-min=10.7"
	BINARYNAME="MacOSX-10.7+"
fi

BINARYLOCALDIR=${VERSIONSTRING}-${BINARYNAME}
BINARYDIR=${PWD}/${BINARYLOCALDIR}
rm -r ${BINARYDIR} ${BINARYDIR}.tgz

echo "${BINARYDIR}"
echo "./configure --prefix=\"${BINARYDIR}\" CFLAGS=\"${ARCHFLAGS} ${CFLAGS}\" CXXFLAGS=\"${ARCHFLAGS} ${CXXFLAGS}\" LDFLAGS=\"${ARCHFLAGS} ${LDFLAGS}\""
./configure --enable-static --prefix="${BINARYDIR}" CFLAGS="${ARCHFLAGS} ${CFLAGS}" CXXFLAGS="${ARCHFLAGS} ${CXXFLAGS}" LDFLAGS="${ARCHFLAGS} ${LDFLAGS}"

make clean
make -j 6
make test
make install


tar -czf ${BINARYLOCALDIR}.tar.gz ${BINARYLOCALDIR}
rm -r ${BINARYLOCALDIR}
