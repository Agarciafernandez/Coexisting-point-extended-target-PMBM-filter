This repository contains Matlab implementations of the following multi-target filters for coexisting point and extended targets:

- Poisson multi-Bernoulli mixture (PMBM) filter.
- Poisson multi-Bernoulli (PMB) filter.
- Multi-Bernoulli mixture (MBM) filter.
- Multi-Bernoulli (MB) filter.

These filters are described in 

A. F. Garcia-Fernandez, J. L. Williams, L. Svensson and Y. Xia, "A Poisson multi-Bernoulli mixture filter for coexisting point and extended targets," in IEEE Transactions on Signal Processing, 2021 doi: 10.1109/TSP.2021.3072006.


https://arxiv.org/abs/2011.04464


The filters are evaluated using the generalised optimal subpattern-assignment (GOSPA) integrated with the Gaussian Wasserstein distance

A. S. Rahmathullah, A. F. Garcia-Fernandez, and L. Svensson, Generalized optimal sub-pattern assignment metric, in 20th International Conference on Information Fusion, 2017.

Video on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM


- main_PMBM_point_extended_MC.m runs the point-extended target PMBM/PMB filters

-main_MBM_point_extended_MC.m runs the point-extended target MBM/MB filter



