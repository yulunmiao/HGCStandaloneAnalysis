###########################
## Example analysis #######
###########################

# make symbolic links to files needed from PFCalEE/include(src):

cd include
ln -s ../../include/G4SiHit.hh .
ln -s ../../include/HGCSSSimHit.hh .
ln -s ../../include/TransverseGeometry.hh .
ln -s ../../include/dict.h .

cd ../src
ln -s ../../src/G4SiHit.cc .
ln -s ../../src/HGCSSSimHit.cc .
ln -s ../../src/TransverseGeometry.cc .
ln -s ../../src/dict.cc .

make

#any .cc and .hh files will be compiled into obj files and into a library
#any .cpp file in test/ will be compiled as executable.

