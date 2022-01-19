# TCrosser

C++ code to compute threshold crossings on a multichannel data stream and align to various minima/maxima.

## Installation

```
mkdir -p build/release
cd build/release
cmake -G "Unix Makefiles" -DBUILD_PLEXON=1 -DCMAKE_BUILD_TYPE=Release ../..
make
make install  # or sudo make install if on Unix
```
