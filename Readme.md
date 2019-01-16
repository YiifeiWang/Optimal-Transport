# Algorithms for large-scale optimal transport

This is a repository of the course project "Algorithms for large-scale optimal transport" for [Convex Optimizatin 2018 Fall](http://bicmr.pku.edu.cn/~wenzw/opt-2018-fall.html). It is a group project by Yifei Wang and Feng Zhu (lexicographically). 



## Numerical experiments

To reproduce the results in our report, please run the following Matlab programs:

- Figure 1 and Figure 2 in Section 3.1.2

  `plot_gmm.m`

- Figure 3 in Section 4.1.3

  `plot_ellipse.m`

- Figure 4 in Section 4.1.4

  `plot_caff.m`

- Perfomance of mosek and gurobi on randomly generated data in Section 4.2.1

  `Test_RGD_mb.m`

- Perfomance of mosek and gurobi on DOTmark in Section 4.2.2

  `Test_DOTmark_mb.m`

- Perfomance of mosek and gurobi on Ellipse Example and Caffarelli’s Example in Section 4.2.3

  `Test_ellipse_mb.m` `Test_caff_mb.m`

- Perfomance of first order methods on randomly generated data in Section 4.3.2

  `Test_RGD_fo.m`

- Perfomance of first order methods on DOTmark in Section 4.3.3

  `Test_DOTmark_fo.m`

- Perfomance of first order methods on Ellipse Example and Caffarelli’s Example in Section 4.3.4

  `Test_ellipse_fo.m` `Test_caff_fo.m`

- Figure 5 and Figure 6 in Section 4.4.1

  `plot_gmm2.m`

- Perfomance of algorithms for entropic regularization of OT on Gaussian mixture model in Section 4.4.1

  `Test_gmm_er.m`

- Perfomance of algorithms for entropic regularization of OT on randomly generated data in Section 4.4.2

  `Test_RGD_er.m`

- Perfomance of algorithms for entropic regularization of OT on DOTmark in Section 4.4.3

  `Test_RGD_er.m`

- Perfomance of algorithms for entropic regularization of OT on  Ellipse Example and Caffarelli’s Example in Section 4.4.4

  `Test_ellipse_er.m` `Test_caff_er.m`



## Explanation of our codes

### Model

 'model_*.m' generates an OT model based on different dataset and 'model_unified.m' is a unified program to generate OT models. 

### LP/LPER

'LP_*.m' implements different algorithms to solve the LP problem. 

'LPER_*.m' implements different algorithms to solve the entropic regularization of OT.  

### OT

'OT_*.m' provides necessary functions to create an OT model. 

### Plot

'Plot_*.m' plots figures in our reports.

### Test

'Test_*.m' implements numerical experiments in our reports.