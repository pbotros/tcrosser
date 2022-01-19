# TCrosser

C++ code to compute threshold crossings on a multichannel data stream and align to various minima/maxima.

## Installation

```
mkdir -p build/release
cd build/release
cmake -G "Unix Makefiles" -DBUILD_PLEXON=1 -DCMAKE_BUILD_TYPE=Release ../..
cd tcrosser
make
make install  # or sudo make install if on Unix
```

The only dependency needed is Python/Cython for this package. If Cython isn't installed, install `cython` via pip3 or conda.
