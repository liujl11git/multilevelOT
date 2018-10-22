Matlab codes for the multilevel optimal transport (OT) algorithm in [1].

[Problem description]

Given two probability distributions (2D), or two images, rho0 and rho1, we want to find a transport between them with minimal energy:

min_{m} ||m||_p,
s.t. divergence(m,dx) = rho0-rho1,
     with zero-flux boundary condition.

p is 1,2, or inifity.
dx is the grid step size.

[Algorithms and Functions]
1. "otFunctions/W1PD_ML.m" implements the Algorithm 1M in [1]. 
2. "otFunctions/W1PDHG_ML.m" implements the Algorithm 2M in [1].

[Syntax]
1. [m,phi] = W1PD_ML(dx, rho0, rho1, p, opts);
2. [m,gphi,phi] = W1PDHG_ML(dx, rho0, rho1, p,opts);

[Input Description]
1. dx, rho0, rho1, p: defined in [Problem Description].
2. opts: parameter options. 
  (1) opts.MaxIter: Upper bound of the number of iterations. Default value is 100000.
  (2) opts.tol: Stopping tolerance for the fixed point residual. Default value is 1e-6 for W1PD_ML and 1e-5 for W1PDHG_ML.
  (3) opts.verbose: Flag for displaying the metrics. Default value is 1.
  (4) opts.L: Number of levels for calculation. Default value is 6.
  (5) opts = [] means we use the default options.

[Output Description]
1. m: The optimal transport between rho0 and rho1.
2. phi: The Kantorovich potential, see [1] for reference.
3. gphi: The gradient of phi.

[Demos]
1. "demo/demo_cat_Algorithm_1M.m" and "demo/demo_img_Algorithm_1M.m" are two demos of "W1PD_ML.m".
2. "demo/demo_cat_Algorithm_2M.m" and "demo/demo_img_Algorithm_2M.m" are two demos of "W1PDHG_ML.m".

[How to run the demos]
1. Enter the root path of this project.
2. Please run "initialize.m."
3. Run the script files in the "demo/" folders.

[Compilation]
1. "W1PD_ML.m" is implemented purely by MATLAB codes. So it can be run without compilation.
2. "W1PDHG_ML.m" is based on C source codes. We have to compile it before running.
3. To compile, please refer the file "how_to_compile".

[Precompiled files]
We provide precompiled files for the following three platforms:
1. Precompiled file "emdPDHG.mexa64" for linux 64 bit system is in "src_c" folder.
2. Precompiled file "emdPDHG.mexw64" for windows 64 bit system is in "src_c" folder.
3. Precompiled file "emdPDHG.mexmaci64" for windows 64 bit system is in "src_c" folder.

[Detailed Notes]

The first attempt to run "W1PDHG_ML" may take longer time than we reported in the paper. The reason is that FFT needs setup time. You may uncomment line 58 in "otFunctions/W1PDHG_ML.m" to display the FFT setup time.

[References]

[1] J. Liu, W. Yin, W. Li, Y.T. Chow, "Multilevel Optimal Transport: a Fast Approximation of Wasserstein-1 distances", submitted, 2018.
[2] M. Jacobs, F. Leger, W. Li, S. Osher, "Solving Large-Scale Optimization Problems with a Convergence Rate Independent of Grid Size", submitted, 2018.
[3] W. Li, E. Ryu, S. Osher, W. Yin, W. Gangbo, "A Parallel Method for Earth Mover's Distance." Journal of Scientific Computing, 2018.
