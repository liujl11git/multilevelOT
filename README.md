Matlab codes for the optimal transport algorithm in [1].

[Authors]
Jialin Liu, liujl11@math.ucla.edu
Wuchen Li, wcli@math.ucla.edu
Matthew Jacobs, majaco@math.ucla.edu

[Note]
1. "otFunctions/W1PD_ML.m" is the Algorithm 1M in [1].
2. "otFunctions/W1PDHG_ML.m" is the Algorithm 2M in [1].
3. "demo/demo_cat_Algorithm_1M.m" and "demo/demo_img_Algorithm_1M.m" are two demos of "W1PD_ML.m".
4. "demo/demo_cat_Algorithm_2M.m" and "demo/demo_img_Algorithm_2M.m" are two demos of "W1PDHG_ML.m".

[Precompiled files]
1. Precompiled file "emdPDHG.mexa64" for linux 64 bit system is in "src_c" folder.
2. Precompiled file "emdPDHG.mexw64" for windows 64 bit system is in "src_c" folder.

[How to compile]
If the above precompiled files do not work for you:
1. "W1PD_ML.m" is implemented purely by MATLAB codes. So it can be run without compilation.
2. "W1PDHG_ML.m" is based on C source codes. We have to compile it before running.
3. To compile, please refer the file "how_to_compile".

[How to run]
1. Enter the current path.
2. Please run "initialize.m."
3. Run the script files in the "demo/" folders.

[References]
[1] J. Liu, W. Yin, W. Li, Y.T. Chow, "Multilevel Optimal Transport: a Fast Approximation of Wasserstein-1 distances", submitted, 2018.
[2] M. Jacobs, F. Leger, W. Li, S. Osher, "Solving Large-Scale Optimization Problems with a Convergence Rate Independent of Grid Size", submitted, 2018.
[3] W. Li, E. Ryu, S. Osher, W. Yin, W. Gangbo, "A Parallel Method for Earth Mover's Distance." Journal of Scientific Computing, 2018.
