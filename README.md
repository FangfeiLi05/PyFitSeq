## FitSeq

### What is FitSeq?

FitSeq is a Python-based fitness estimation tool for pooled amplicon sequencing studies. FitSeq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


### Installing
* Tested on MacBook Pro (3.1 GHz Intel Core i5), with Python 3.7.4.
* Clone this repository by running `git clone https://github.com/FangfeiLi05/FitSeq.git` in terminal.
* `cd` to the root directory of the project (the folder containing README.md)
* Install dependencies by running `pip install -r requirements.txt` in terminal.

### Evolution Simulations
Models competative pooled growth of a population of genotypes with different fitnesses. This simulation may include many sources of noise, including growth noise, noise from cell transfers, DNA extraction, PCR, and sequencing.

#### OPTIONS
+ `--input` or `-i`: a .csv file, with the 1st column being fitness of each genotype (x1 x2 ...), and the 2nd column being initial cell number of each genotype at generation 0 (n1 n2 ...)
+ `--t_seq` or `-t`: sequenced time-points (0 t1 t2 ...)
+ `--read_num_average_seq` or `-r`: average number of reads per genotype per sequencing time-point, [0, r1, r2, ...]
+ `--noise_option` or `-n`: five types of possible noise (cell growth, bottleneck transfer, DNA extraction, PCR, sequencing (default: growth bottleneck_transfer DNA_extraction PCR sequencing)
+ `--dna_copies` or `-d`: average copy number of genome DNA per genotype as template in PCR (default: 500)
+ `--pcr_cycles` or `-p`: number of cycles in PCR (default: 25) 
+ `--output_filename` or `-o`: prefix of output .csv files (default: output)

#### OUTPUTS
+ `output_filename_EvoSimulation_Read_Number.csv`: read number per genotype per sequencing time-point
+ `output_filename_EvoSimulation_Mean_Fitness.csv`: mean fitness per sequencing time-point
+ `output_filename_EvoSimulation_Input_Log.csv`: a record of all inputs

#### EXAMPLES
```
python evo_simulator.py --help

python evo_simulator.py -i Input.csv -t 0 3 6 9 12 -r 50 50 50 50 50 -o result
```      

### Fitness Estimation

#### OPTIONS

#### OUTPUTS

#### EXAMPLES




