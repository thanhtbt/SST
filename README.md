# OPIT: A simple but effective sparse subspace tracking method


In this work, we propose a new provable effective method called OPIT (which stands for Online Power Iteration via Thresholding)
for tracking the sparse principal subspace of data streams over time. Particularly, OPIT introduces a new adaptive variant of
power iteration with space and computational complexity linear to the data dimension. In addition, a new column-based thresholding operator is developed to regularize the subspace sparsity. Utilizing both advantages of power iteration and thresholding operation, OPIT is capable of tracking the underlying subspace
in both classical regime and high dimensional regime.  


## Demo
+ Run demo_xyz.m for synthetic data


## Experimental Results

+ OPIT vs The State-of-the-art Subspace Tracking Algorithms 

![compare](https://user-images.githubusercontent.com/26319211/203633980-0870621d-9cf0-4166-8838-442c4dfad109.jpg)

## Reference

This code is free and open source for research purposes. If you use this code, please acknowledge the following paper.

[1] L.T. Thanh, K. Abed-Meraim, N.L. Trung, A. Hafiance. "[*Sparse Subspace Tracking in High Dimensions*](https://ieeexplore.ieee.org/document/9746546)". **Proc. 47th IEEE ICASSP**, 2022. [[PDF](https://thanhtbt.github.io/files/2022_ICASSP%20-%20Sparse%20Subspace%20Tracking%20in%20High%20Dimensions.pdf)].
