## Example analysis

#to compile
mkdir lib
mkdir bin
mkdir obj

make

#any .cc and .hh files will be compiled into obj files and into a library in lib/
#any .cpp file in test/ will be compiled as executable, in bin/

./bin/studyTransverseProperties -i ../ -v 20 -e 5
./bin/compareCaloStackPerformances -i ../ -v 20 -e 5
