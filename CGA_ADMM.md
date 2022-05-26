# CGA_ADMM
This code is the implementation of algorithms in our paper: CGA_ADMM: An efficient constraint generation based alternating direction method of multipliers for the metric nearness problem.

## Requirement
1. Microsoft Visual Studio 2019 or newer version. Our code uses Microsoft Concurrency Runtime thus this code can only be compiled by MSVC.
2. Eigen library 3.4 or newer version: https://eigen.tuxfamily.org/index.php?title=Main_Page
3. Intel MKL

## Build
1. Add the directory of Eigen library to the include path.
2. Add CGA_ADMM.lib to the linker.
3. Switch the option "Use oneMKL" to "Parallel". Note that if MKL library is not installed, add `/D "NO_MKL"` to the compiling command instead, which may get a slower version software.
3. Compile with MSVC.

## Run
Make sure CGA_ADMM.dll and CGA_ADMM.exe are in the same directory.
```
SYNOPSIS
        admm <norm> <filename>... [--maxits <max_iter>] [--ikt <inner_kkt_tol>] [--kt <kkt_tol>] [--gt <gaptol>] [--sig <sigma>] [--tau <tau>] [--mxvio] [-i <input directory>]

        dykstra <norm> <filename>... [--maxits <max_iter>] [-ct <contol>] [-gt <gaptol>] [-sf <statusfreq>] [--gam <gamma>] [-i <input directory>]


OPTIONS
        <norm>      0 for l_inf, 1 for l_1, 2 for l_2
        <max_iter>  maximum value of iterations, default:10000
        <inner_kkt_tol>
                    kkt tolerance of inner iterations, default:0.001
        <kkt_tol>   total kkt tolerance, default:0.001
        <gaptol>    tolerance of relative gap, default:0.01
        <sigma>     init value of sigma, default:1
        <tau>       value of tau, default:1.618
        --mxvio     add max(Ax-b)<=1e-2 to stopping criteria
        <max_iter>  maximum value of iterations, default:10000
        <contol>    tolerance of violation for constraints, default:0.01
        <statusfreq>frequence of checking status and exit, default:10
        <gamma>     parameter gamma for l_1 and l_inf, note that for l_2 gamma is fixed to 0.5, for
                    l_1 and l_inf default value=1
```

## Examples
To calculate l_1 norm problem by ADMM:

`CGA_ADMM.exe admm 1 minnesota.txt`

To calculate l_2 norm problem by ADMM, and add max(Ax-b)<=1e-2 to stopping criteria:

`CGA_ADMM.exe admm 2 minnesota.txt --mxvio`

To calculate l_\inf norm problem by Dykstra's Projection method, with multiple input files, and set gamma to 200:

`CGA_ADMM.exe dykstra 0 minnesota.txt ca-GrQc.txt --gam 200`

To calculate l_\inf norm problem by ADMM, with multiple input files in ./dataset directory, and set max_iter to 500:

`CGA_ADMM.exe admm 0 minnesota.txt ca-GrQc.txt -i dataset --maxits 500`
