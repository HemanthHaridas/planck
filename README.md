![example event parameter](https://github.com/HemanthHaridas/planck_cpp/actions/workflows/cmake-multi-platform.yml/badge.svg)  ![](https://tokei.rs/b1/github/HemanthHaridas/planck_cpp)

### Planck

<p align="justify"> Planck is a rewrite of the previous project <a href="https://github.com/HemanthHaridas/plank.py"> plank.py </a>, which was a mixture of CPython and Cython. The older version suffered from poor performance due to its design and my limited programming skills at the time of its creation.</P>

<p align="justify"> Planck, on the other hand, is written entirely in C++ and incorporates OpenMP parallelism to improve performance in critical sections. However, it is important to note that Planck is not intended to replace commercially available software like Gaussian, GAMESS, or NWChem. Instead, it serves as a set of programs suitable for educational purposes. Planck is distributed under the GPL3 License. </p>

#### Usage Instructions

<p align="justify"> To run Planck, specify an input file and call it using the following command: </p>

```bash 
planck input_file > output_file
``` 

<p align="justify"> It is mandatory to provide an input file when running Planck </p>