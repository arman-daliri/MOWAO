# MOWAO 
MOWAO: A Multi-Objective metahuristic based on Water Optimization Algorithm

## Files
The algorithm implementation is in the `mowao` folder.  
The source code for the Python package is in the `src` folder.  
In the `example` folder, there are three examples that show how you can use the algorithm.

## How to run
To run the Python examples, you should first install the package and then run the commands below.
```sh
$ git clone https://github.com/arman-daliri/MOWAO
$ cd MOWAO/example/
$ python example1.py # python example2.py
```
To run the C example, you should have a C compiler like GCC or clang. And then in the example folder, run
```sh
$ make
$ ./example # example.exe
```
If you have clang instead of gcc, you should change the `CC = gcc` in `Makefile` to `CC = clang`.

## Installation
You can install the Python package from PyPI.
``` sh
$ pip install mowao
```
Or from this repository
``` sh
$ git clone https://github.com/arman-daliri/MOWAO
$ cd MOWAO
$ pip install . 
```

## Authors
- Mohammad Mahdi Mir [Contact](mailto:m.m.mir83hc@gmail.com)
- Arman Daliri [Contact](mailto:daliriwork2@gmail.com)
- Mahdieh Zabihimayvan [Contact](mailto:Zabihimayvan@ccsu.edu)

