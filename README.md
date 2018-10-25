This is the Matlab code repo for Optimal Transport (OT) by a multilevel method that is introduced in [1].

OT plays crucial roles in many areas, including fluid dynamics, image processing, machine learning, and control. It is a well-posed distance measure of two probability distributions. This distance is also called Earth Mover’s distance (EMD) or the Wasserstein distance.

Our Matlab code computes the Wasserstein-1 distance between two distributions defined on a grid. It uses one of the following three ground measures: 1-norm, 2-norm, inf-norm.

Our code has the following advantages:
* Fast. It takes a few seconds on a single CPU to compute optimal transport on a 1024*1024 grid.
* User friendly.

Problem description
===================
Given two 2D probability distributions, or two images, rho0 and rho1, we want to find a transport from one to the other with the minimal energy:

![equation](https://latex.codecogs.com/gif.latex?%5Cmin_m%20%7E%7E%7E%5Cint_x%20%5C%7Cm%28x%29%5C%7C_p%20dx%5C%5C%20%5Ctext%7Bsubject%20to%7D%7E%7E%20%5Ctext%7Bdivergence%7D_h%28m%29%20%3D%20%5Crho%5E0%20-%20%5Crho%5E1%2C%5C%5C%20%7E%7E%7E%7E%7E%7E%7E%7E%7E%7E%7E%7E%7E%7E%7E%7E%5Ctext%7Bunder%20the%20zero-flux%20boundary%20condition%7D.)

Here, p can be 1, 2, or inifity, so ![equation](https://latex.codecogs.com/gif.latex?%5C%7Cm%28x%29%5C%7C_p) is the 1, 2, or infinity norms of m(x), respectively,
and h is the grid step size.

Algorithms and Functions
========================
1. "otFunctions/W1PD_ML.m" implements Algorithm 1M in [1];
2. "otFunctions/W1PDHG_ML.m" implements Algorithm 2M in [1].

## Syntax
```
[m,phi] = W1PD_ML(h, rho0, rho1, p, opts);

[m,gphi,phi] = W1PDHG_ML(h, rho0, rho1, p,opts);
```

## Input Description
1. h, rho0, rho1, p: defined in [Problem Description].
2. opts: parameter options. 
  * opts.MaxIter: Upper bound of the number of iterations. Default value is 100000.
  * opts.tol: Stopping tolerance for the fixed point residual. Default value is 1e-6 for W1PD_ML and 1e-5 for W1PDHG_ML.
  * opts.verbose: Flag for displaying the metrics. Default value is 1.
  * opts.L: Number of levels for calculation. Default value is 6.
  * set opts = [] to use the default options.

## Output Description
1. m: The optimal transport between rho0 and rho1.
2. phi: The Kantorovich potential; see [1] for reference.
3. gphi: The gradient of phi.

Demos
=====
1. "demo/demo_cat_Algorithm_1M.m" and "demo/demo_img_Algorithm_1M.m" are two demos that call "W1PD_ML.m".
2. "demo/demo_cat_Algorithm_2M.m" and "demo/demo_img_Algorithm_2M.m" are two demos that call "W1PDHG_ML.m".

## How to run the demos
1. Enter the root path of this project.
2. Run "initialize.m."
3. Run the script files in the "demo/" folders.

Compiling
===========
1. "W1PD_ML.m" is implemented purely in MATLAB. No compiling is needed.
2. "W1PDHG_ML.m" calls an external mex code. You can either use the precompiled files (below) or compile its C source code (cf. the file "how_to_compile").

## Precompiled files
We provide the following precompiled files in the folder "src_c":
1. For 64-bit Linux: use the precompiled file "emdPDHG.mexa64";
2. For 64-bit Windows: use the precompiled file "emdPDHG.mexw64";
3. For 64-bit Mac: use the precompiled file "emdPDHG.mexmaci64".

First call takes longer
==============
When you run "W1PDHG_ML" for the first time, some FFT codes takes an extra time to set itself up. You can uncomment Line 58 in "otFunctions/W1PDHG_ML.m" to display the FFT setup time.

References
==========
[1] Liu, J., Yin, W., Li, W., & Chow, Y. T. (2018). Multilevel Optimal Transport: a Fast Approximation of Wasserstein-1 distances. arXiv preprint arXiv:1810.00118.

[2] Jacobs, M., Léger, F., Li, W., & Osher, S. (2018). Solving Large-Scale Optimization Problems with a Convergence Rate Independent of Grid Size. arXiv preprint arXiv:1805.09453.

[3] Li, W., Ryu, E. K., Osher, S., Yin, W., & Gangbo, W. (2018). A parallel method for earth mover’s distance. Journal of Scientific Computing, 75(1), 182-197.
