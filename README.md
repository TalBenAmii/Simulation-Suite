# How to build
First go to desired tool via:
cd <tool-name>
then:
mkdir build
cd build
cmake ..
make

on windows:

mkdir build
cd build/
cmake ..
cmake --build . --config Release

exe should be inside /build/Release folder
