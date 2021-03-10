set "CC=gcc.exe"
set "CXX=g++.exe"
set "FORTRAN=%LIBRARY_BIN%\gfortran.exe"

cmake %SRC_DIR% ^
    -DCMAKE_BUILD_SHARED=ON ^
    -DCMAKE_CXX_FLAGS="--std=c++14 -D_USE_MATH_DEFINES -DCONDA_BUILD_M2W64_TOOLCHAIN" ^
    -G "MinGW Makefiles" ^
    -DCMAKE_INSTALL_PREFIX="%PREFIX%"

if errorlevel 1 exit 1

cmake --build . --target install
if errorlevel 1 exit 1
