cd code/deps/kahypar
ls
mkdir build || cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
sudo make install.library