# MOCHI

## Required Package
MOCHI requires:
* Python (tested 2.7.13)
* numpy (tested 1.13.1)
* pandas (tested 0.20.3)

To recreate certain figures in the paper. Requires:
* Matplotlib(2.0.2)
* Seaborn (0.8)


## Usage
The parameters are as followed:
* "-c","--cell",default = 'gm12878' :Controls the data it runs
* "-p", "--parent_threshold",default = 0.45 :Controls the threshold of the overlapping TF
* "-t", "--total",default = False :Whether to run the code all over again. By default it uses the cached data to accelerate.
* "-m", "--motif",default = 4 :4-node motif or 3-node motif
