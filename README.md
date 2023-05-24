# IBR-Data-Analytics

This repository has the following data analytics algorithms for post processing. 
1. FFT: The code is in Julia. Input data: PSCAD files of time-series data in CSV format; output data: phasors at the measuring frequency point. The particular example is to extract dq admittance. 
2. Prony: MATLAB code
3. ERA: MATLAB code
4. Matrix Pencil: MATLAB code
5. Dynamic mode decomposition: MATLAB code
6. Frequency-domain data fitting: MATALB code. MATLAB's system identification toolbox is needed. Input: frequency-domain data (nonparametric). Output: transfer functions (parametric). 

Algorithms 2-5 have input data as time-series data with equal intervals. The output are the modes or eigenvalues, and/or the mode shapes or the eigenvectors.  
