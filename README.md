[![Python 3.7](https://img.shields.io/badge/python-3.7-green.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Contact Info](https://img.shields.io/badge/Email-fangfeili0525@gmail.com-orange.svg)]()


## PyFitSeq

### 1. What is PyFitSeq?

PyFitSeq is a Python-based fitness estimation tool for pooled amplicon sequencing studies. PyFitSeq is Python re-coded version of the MATLAB tool FitSeq https://github.com/sashaflevy/Fit-Seq. If you use this software, please reference: [F. Li, et al. Unbiased Fitness Estimation of Pooled Barcode or Amplicon Sequencing Studies. Cell Systems, 7: 521-525 (2018)](https://www.sciencedirect.com/science/article/pii/S2405471218303909?via%3Dihub). PyFitSeq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

It currently has two main functions:
* `evo_simulator.py` performs simulations of competitve pooled growth of a population of genotypes.
* `pyfitseq.py` calculates the fitness of each genotype from read-count time-series data.
    
A walk-through is included as the jupyter notebook [here](https://github.com/FangfeiLi05/PyFitSeq/blob/master/PyFitSeq_Walk_Through.ipynb).



### 2. How to install PyFitSeq?
* Python 3 is required. This version has been tested on a MacBook Pro (3.1 GHz Intel Core i5), with Python 3.7.4.
* Clone this repository by running `git clone https://github.com/FangfeiLi05/PyFitSeq.git` in terminal.
* `cd` to the root directory of the project (the folder containing `README.md`).
* Install dependencies by running `pip install -r requirements.txt` in terminal.
* Install pyfitseq by running `pip install .`

**OR**
* Run `python3 -m pip install git+https://github.com/darachm/PyFitSeq.git` to install
  without cloning the repository.

#### 2a Alternative, use in a Singularity container

With the closing of Singularity Hub there aren't yet publicly available
containers for this, but you can build your own with a command like:

    singularity build fitseq-latest.simg Singularity.fitseq-latest

Then you can run on any 
[computer running Singularity](https://sylabs.io/guides/3.8/user-guide/quick_start.html#quick-installation-steps),
such as your local HPC, using a command like:

    singularity exec fitseq-latest.simg pyfitseq.py -h

### 3. How to use PyFitSeq?

#### 3.1. Evolution Simulation
`evo_simulator.py` models competative pooled growth of a population of genotypes with different fitnesses. This simulation can be made to include sources of noise, including growth noise, noise from cell transfers, DNA extraction, PCR, and sequencing.

##### OPTIONS
* `--input` or `-i`: a .csv file, with
  + 1st column of .csv: fitness of each genotype, [x1, x2, ...]
  + 2nd column .csv: initial cell number of each genotype at generation 0, [n1, n2, ...]
* `--t_seq` or `-t`: time-points evaluated in number of generations (`format: 0 t1 t2 ...`)
* `--read_num_average_seq` or `-r`: average number of reads per genotype for each time-point (`format: 0 r1 r2 ...`)
* `--noise_option` or `-n`: which types of noise to include in the simulation, default is all sources of noise (`default: growth bottleneck_transfer DNA_extraction PCR sequencing`)
* `--dna_copies` or `-d`: average genome copy number per genotype used as template in PCR (`default: 500`)
* `--pcr_cycles` or `-p`: number of cycles of PCR (`default: 25`) 
* `--fitness_type` or `-f`: type of fitness: Wrightian fitness (w), or Malthusian fitness (m)' (`default: m`)
* `--output_filename` or `-o`: prefix of output .csv files (`default: output`)

##### OUTPUTS
* `output_filename_EvoSimulation_Read_Number.csv`: read number per genotype for each time-point
* `output_filename_EvoSimulation_Mean_Fitness.csv`: mean fitness for each time-point
* `output_filename_EvoSimulation_Input_Log.csv`: a record of all inputs

##### For Help
```
python evo_simulator.py --help
```

##### Examples
```
python evo_simulator.py -i input_EvoSimulation.csv -t 0 3 6 9 12 -r 50 50 50 50 50 -o output
python evo_simulator.py -i input_EvoSimulation.csv -t 0 2 4 6 8 -r 75 75 75 75 50 -n DNA_extraction PCR sequencing -d 300 -p 27 -f w -o output
```      


#### 3.2. Fitness Estimation
`pyfitseq.py` estimates the fitness of each genotype from read-count time-series data.

##### OPTIONS
* `--input` or `-i`: a .csv file, with each column being the read number per genotype at each sequenced time-point
* `--t_seq` or `-t`: sequenced time-points in number of generations (`format: 0 t1 t2 ...`)
* `--max_iter_num` or `-m`: maximum number of iterations in the optimization (Small numbers can reduce running time and decrease accuracy.) (`default: 10`)
* `--kappa` or `-k`: a noise parameter that characterizes the total noise introduced by growth, cell transfer, DNA extraction, PCR, and sequencing (To measure kappa empirically, see the reference: [S. F. Levy, et al. Quantitative Evolutionary Dynamics Using High-resolution Lineage Tracking. Nature, 519: 181–186 (2015)].) (`default: 2.5`)
* `--regression_num` or `-g`: number of points used in the initial linear-regression-based fitness estimate (`default: 2`)
* `--fitness_type` or `-f`: type of fitness: Wrightian fitness (w), or Malthusian fitness (m) (`default: m`)
* `--output_filename` or `-o`: prefix of output .csv files (`default: output`)

##### OUTPUTS
* `output_filename_FitSeq_Result.csv`: a .csv file, with
  + 1st column of .csv: estimated fitness of each genotype, [x1, x2, ...]
  + 2nd column of .csv: log likelihood value of each genotype, [f1, f2, ...]
  + 3rd column of .csv: estimated mean fitness per sequenced time-point, [x_mean(0), x_mean(t1), ...]
  + 4th+ columns of .csv: estimated read number per genotype per time-point, with each time-point being a column

##### For Help
```
python pyfitseq.py --help
```  

##### Examples
```
python pyfitseq.py -i output_EvoSimulation_Read_Number.csv -t 0 3 6 9 12 -o output
python pyfitseq.py -i output_EvoSimulation_Read_Number2.csv -t 0 2 6 8 -m 12 -k 2 -g 3 -f w -o output
```  



