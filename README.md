# OPIT: A Simple but Effective Sparse Subspace Tracking Method

In this work, we propose a new provable effective method called OPIT (which stands for Online Power Iteration via Thresholding) for tracking the sparse principal subspace of data streams over time. Particularly, OPIT introduces a new adaptive variant of power iteration with space and computational complexity linear to the data dimension. In addition, a new column-based thresholding operator is developed to regularize the subspace sparsity. Utilizing both advantages of power iteration and thresholding operation, OPIT is capable of tracking the underlying subspace
in both classical and high dimensional regimes.  


## Demo
Run 
+ "demo_effect_forgetting_factor.m": To illustrate the effect of the forgetting factor on the performance of OPIT
+ "demo_noise_effect.m": To illustrate the effect of noise on the performance of OPIT
+ "demo_nonstationary.m": To illustrate the performance of OPIT in nonstationary environments
+ "demo_low_dimension_comparison.m": To illustrate the performance of subspace tracking algorithms in the classical setting 
+ "demo_high_dimension_comparison.m": To illustrate the performance of subspace tracking algorithms in high dimension


## State-of-the-art algorithms for comparison

+ **LORAF**:  P. Strobach, “[*Low-rank adaptive filters*](https://ieeexplore.ieee.org/abstract/document/553469/)”. **IEEE Trans. Signal Process.**, 1996.
+ **OPAST**:  K. Abed-Meraim et al. “[*Fast orthonormal PAST algorithm*](https://ieeexplore.ieee.org/abstract/document/823526/)”. **IEEE Signal Processing Letters**, 2000.
+ **FAPI**: R. Badeau et al. “[*Fast approximated power iteration subspace tracking*](https://ieeexplore.ieee.org/abstract/document/1468483/)”. **IEEE Trans. Signal Process.**, 2015.
+ **GFAPI**: M. Arjomandi-Lari et al. “[*Generalized YAST algorithm for signal subspace tracking*](https://www.sciencedirect.com/science/article/pii/S0165168415001607)”. **Signal Process.**, 2015.
+ **SSPCA**: W. Yang and H. Xu, “[*Streaming sparse principal component analysis*](https://proceedings.mlr.press/v37/yangd15.html)”. **Proc. ICML**, 2015.
+ **L1-PAST**: X. Yang et al. “[*Fast STAP method based on PAST with sparse constraint for airborne phased array radar*](https://ieeexplore.ieee.org/abstract/document/7470515/)”. **IEEE Trans. Signal Process.**, 2016.
+ **SSFAPI**: N. Lassami et a. “[Low cost sparse subspace tracking algorithms](https://www.sciencedirect.com/science/article/pii/S0165168420300657)”. **Signal Process.**, 2020.
+ **PETRELS-ADMM**: L.T. Thanh et al. "[*Robust Subspace Tracking with Missing Data and Outliers: Novel Algorithm with Convergence Guarantee*](https://ieeexplore.ieee.org/document/9381678)". **IEEE Trans. Signal Process.**, 2021.


## Some Experimental Results

+ Effect of the sparsity level

  ![SST_Sparse](https://github.com/thanhtbt/SST/assets/26319211/108f3b0f-6648-4fb6-bfed-bebd39c75bc3)

+ OPIT vs SOTA Algorithms in High Dimension
  
![SST_Compare](https://github.com/thanhtbt/SST/assets/26319211/da5edddd-b9c7-4695-bf6c-820ba4264252)

+ Performance of SST algorithms with different data dimensions and sample sizes

![SST_Compare_v2](https://github.com/thanhtbt/SST/assets/26319211/5fcc2299-b0ca-414a-9e71-5ea5a25f1aa6)

## Reference


This code is free and open source for research purposes. If you use this code, please acknowledge the following paper.

[1] L.T. Thanh, K. Abed-Meraim, N.L. Trung, A. Hafiance. "[*Sparse Subspace Tracking in High Dimensions*](https://ieeexplore.ieee.org/document/9746546)". **Proc. 47th IEEE ICASSP**, 2022. [[PDF](https://thanhtbt.github.io/files/2022_ICASSP%20-%20Sparse%20Subspace%20Tracking%20in%20High%20Dimensions.pdf)].

[2] L.T. Thanh, K. Abed-Meraim, N. L. Trung, & A. Hafiane. "[*OPIT: A Simple and Effective Method for Sparse Subspace Tracking in High-dimension and Low-sample-size Context*](https://doi.org/10.1109/TSP.2023.3349070)". To appear in **IEEE Trans. Signal Process.**, 2024. 
