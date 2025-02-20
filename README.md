# Automatic Clustering Using an Improved Differential Evolution Algorithm

## Overview
This repository contains a MATLAB implementation of **Automatic Clustering Using an Improved Differential Evolution Algorithm (ACDE)**. The algorithm is designed to automatically determine the optimal number of clusters and cluster data efficiently. It has been tested on the **Iris dataset**.

## Paper Reference
This implementation is based on the paper:

> "Automatic clustering using an improved differential evolution algorithm," IEEE International Conference on Systems, Man and Cybernetics, 2007.  
> **Link:** [IEEE Xplore](https://ieeexplore.ieee.org/document/4390004)

## Features
- Uses **Differential Evolution (DE)** for clustering.
- Determines the optimal number of clusters **automatically**.
- Tested on the **Iris dataset**.
- Produces clustering results for multiple runs to ensure robustness.

## Installation
1. Clone the repository:
   ```sh
   git clone https://github.com/your-repo/acde-clustering.git
   cd acde-clustering
   ```
2. Open MATLAB and navigate to the project directory.

## Usage
Run the following command in MATLAB:
```matlab
acde_clustering()
```
This will execute the clustering algorithm on the Iris dataset and display the results.

## Files
- **acde_clustering.m** - Main script implementing the ACDE clustering algorithm.
- **mutation.m, crossover.m, fitness.m, etc.** - Helper functions for the algorithm.

## Results
The algorithm outputs visualized clustering results for 10 independent runs, highlighting the identified clusters and their centroids.

## License
This project is released under the MIT License.

---
Feel free to contribute or raise issues if you encounter any problems! 🚀

