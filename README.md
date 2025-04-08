<!-- need to use divs to center images -->
<div align="center">
  <img src="./docs/images/planck.png">
</div>

### Planck

<p align="justify"> Planck is a rewrite of the previous project <a href="https://github.com/HemanthHaridas/plank.py"> plank.py</a>, which was a mixture of CPython and Cython. The older version suffered from poor performance due to its design and my limited programming skills at the time of its creation.</P>

<p align="justify"> Planck, on the other hand, is written entirely in C++ and incorporates OpenMP parallelism to improve performance in critical sections. However, it is important to note that Planck is not intended to replace commercially available software like Gaussian, GAMESS, or NWChem. Instead, it serves as a set of programs suitable for educational purposes. Planck is distributed under the GPL3 License. </p>

#### Installation Instructions
<p align="justify"> Before installing planck, you must install boost which is a dependency for Planck. </p>

```
git clone https://github.com/HemanthHaridas/planck.git
cd planck
mkdir build
cmake ..
make
make install
```

#### Usage Instructions

<p align="justify"> To run Planck, specify an input file and call it using the following command: </p>

```bash 
planck input_file > output_file
``` 
