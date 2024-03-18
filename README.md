![example event parameter](https://github.com/HemanthHaridas/planck_cpp/actions/workflows/cmake-multi-platform.yml/badge.svg)  ![](https://tokei.rs/b1/github/HemanthHaridas/planck_cpp)

### Planck

<p style="text-align:justify"> Planck is a completely written version of the older <a href=""> plank.py </a> project of mine, which was written in a mixture of CPython and Cython. It was terribly slow, and my programming skills were also not great when I wrote plank.py two years back.</P>

<p style="text-align:justify"> Planck is written in pure C++ and have OpenMP parallelism included for performance critical sections. However, Planck is not designed to be a drop-in replacement for the commercially available codes like Gaussian, GAMESS or NWChem. Instead it is designed to be set of programs that can be used for pedagogical excercises. So, Planck is distributed under a GPL3 License.</p>

#### Usage Details

There are two ways to run planck. 

1. You can either specify an input file and call planck by using the following command

```bash 
planck input_file > output_file
``` 

2. you can also call planck without specifying an input file.  

``` bash 
planck > output_file 
``` 

In this case, The code will expect the input file to called **CORD**. 

It is recommended that you explicitly specify the input file for the calculation. 