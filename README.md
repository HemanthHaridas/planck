![](https://github.com/HemanthHaridas/planck_cpp/actions/workflows/cmake-multi-platform.yml/badge.svg) 
### Planck

<p align="justify"> Planck is a rewrite of the previous project <a href="https://github.com/HemanthHaridas/plank.py"> plank.py</a>, which was a mixture of CPython and Cython. The older version suffered from poor performance due to its design and my limited programming skills at the time of its creation.</P>

<p align="justify"> Planck, on the other hand, is written entirely in C++ and incorporates OpenMP parallelism to improve performance in critical sections. However, it is important to note that Planck is not intended to replace commercially available software like Gaussian, GAMESS, or NWChem. Instead, it serves as a set of programs suitable for educational purposes. Planck is distributed under the GPL3 License. </p>

#### Installation Instructions
<p align="justify"> Before installing planck, you must install boost and <a href="https://github.com/mcodev31/libmsym.git"> mysm</a>, both of which are dependencies for Planck. </p>

##### Install instructions for msysm
```
git clone https://github.com/HemanthHaridas/planck_cpp.git
cd planck_cpp/src
mkdir external

git clone https://github.com/mcodev31/libmsym.git
cd libmsym
mkdir build install
cd build
cmake -DBUILD_SHARED_LIBS:BOOL=OFF -DMSYM_BUILD_EXAMPLES:BOOL=OFF -DCMAKE_INSTALL_PREFIX:PATH=<full path to the install directory that was created> ../.
make
make install
```

##### Install instructions for planck
```
cd planck_cpp
g++ src/base/planck_*.cpp src/lookup/planck_lookup.cpp src/planck_main.cpp -L<your path to lib inside libmsym install directory> -lmsym -o planck.o --std=c++17
```
#### Usage Instructions

<p align="justify"> To run Planck, specify an input file and call it using the following command: </p>

```bash 
planck input_file > output_file
``` 
