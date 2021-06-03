## Mining Order-Preserving Submatrices Under Data Uncertainty: A Possible-World Approach and Efficient Approximation Methods

This is a algorithm for Order-Preserving Submatrices Mining(OPSM) based on the possible world semantics. OPSM means if there is a permutation of the columns of S, under which the entry values of each row in S are strictly increasing. The possible world semantics
means the data are often noisy and assumes that each matrix entry is a random variable following certain distributions. 

And in this algorithm, we converts the replicates for each entry in the matrix into an interval-based data model. And we are interested in those long patterns supported by sufficient rows, which exhibitstatistical significance rather than occurring by chance. We assume the data are follow uniform distributions(U) or Gaussian Distribution(G) and implement two version algorithms. In each algorithm, we design two different definitions of OPSM significance: (S1) expected support(ES) and (S2) probabilistic frequentness(PF). Using the above definitions to check the significance of each pattern P, two pattern-growth based OPSM mining algorithms are proposed using different techniques, (1) prefix-projection(We mark it as DFS) and (2) Apriori. And we also design two approximation techniques, Poisson Distribution(PFP) and Gaussian Distribution(PFG, also PFA) to speed up the examination of pattern frequentness. More information can see our paper.

There are 8 applications on top of our framework:

* Apri-ES: the exact algorithm using expected support definitions and Apriori techniques
* DFS-ES: the exact algorithm using expected support definitions and prefix-projection techniques
* Apri-PF: the exact algorithm using probabilistic frequentness definitions and Apriori techniques
* DFS-PF: the exact algorithm using probabilistic frequentness definitions and prefix-projection techniques
* Apri-PFP: the approximation algorithm based on Poisson Distribution using expected support definitions and Apriori techniques
* DFS-PFP: the approximation algorithm based on Poisson Distribution using probabilistic frequentness definitions and prefix-projection techniques
* Apri-PFA: the approximation algorithm based on Gaussian Distribution using expected support definitions and Apriori techniques
* DFS-PFA: the approximation algorithm based on Gaussian Distribution using probabilistic frequentness definitions and prefix-projection techniques


## code structure
**Uniform**: The uniform distributions version algorithm, which assumes the data entry follow uniform distributions. 
**Gaussian**: The the Gaussian distributions version algorithm, which assumes the data entry follow gaussian distributions.
**data**: The data folder. We provide a RFID dataset. 

The codes can run in both Window and Linux system.

## Uniform algorithm

### input foramt:

row col

l,r,l,r,***

l,r,l,r,***

****

the *row* and *col* refer to the mining matrix's size. And each matrix entry is defined by a set of [l,r] below the first line.

### How to compile:

g++ -std=c++11 -O2 -m64 -c -o Fourier.o Fourier.cpp

g++ -std=c++11 -O2 -m64 -o opsm Fourier.o ydTest.cpp -lpthread

### How to run:
For Method Apri_ES and DFS_ES:

opsm.exe(./opsm) input_file_name method_name min_row min_col

For Method: Apri_PF  DFS_PF  Apri_PFA  DFS_PFA   Apri_PFP  DFS_PFP

opsm.exe(./opsm) input_file_name method_name min_row min_col th_prob

Please refer to the paper/read_me.txt for method names and more info." << endl;



## Gaussian algorithm

### input foramt:

row col replicates

x,x,**,y,y,***

x,x,**,y,y,***

****

the *row* and *col* refer to the mining matrix's size. The *replicates* means how many times matrix entry repeats. And each matrix entry is defined by a set of [x,x,***] below the first line.

### How to compile:

g++ -std=c++11 -O2 -m64 -c -o Fourier.o Fourier.cpp

g++ -std=c++11 -O2 -m64 -o opsm Fourier.o ydTest.cpp -lpthread

### How to run:
For Method Apri_ES and DFS_ES:

opsm.exe(./opsm) input_file_name method_name min_row min_col

For Method: Apri_PF  DFS_PF  Apri_PFA  DFS_PFA   Apri_PFP  DFS_PFP

opsm.exe(./opsm) input_file_name method_name min_row min_col th_prob


### Contact
Wenwen qu

Email: wenwenqu@stu.ecnu.edu.cn

### Contributors
Ji Cheng

Da Yan 

Wenwen Qu

