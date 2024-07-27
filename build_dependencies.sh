cd code/deps/kahypar
mkdir build
cmake -B build -DCMAKE_BUILD_TYPE=RELEASE
cmake --build build
cd build
sudo make install.library