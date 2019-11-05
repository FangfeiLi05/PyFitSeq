import numpy as np
import pandas as pd
import itertools
import csv
import argparse


def main():
    # ------------------------------------------------------------------------------------------------------------------
    # SIMULATED COMPETITIVE POOLED GROWTH OF A POPULATION OF GENOTYPES WITH DIFFERENT FITNESS (WRIGHTIAN FITNESS).
    # THESE SIMULATIONS INCLUDE EXPERIMENTAL NOISE SUCH AS GROWTH NOISE, SAMPLING DURING BOTTLENECKS, DNA EXTRACTION,
    # PCR, AND SAMPLING ON SEQUENCER.
    #
    # OPTIONS
    # --input: a .csv file, with the 1st column being fitness of each genotype, [x1,x2,...], and the 2nd column being
    #           initial cell number of each genotype at generation 0, [n1,n2,...]
    # --t_seq: sequenced time-points (format: 0 t1 t2 ...)
    # --read_num_average_seq: average number of reads per genotype per sequencing time-point (format: 0 r1 r2 ...)
    # --noise_option: five types of possible noise (cell growth, bottleneck transfer, DNA extraction, PCR, sequencing)
    #                  (default: growth bottleneck_transfer DNA_extraction PCR sequencing)
    # --dna_copies: average copy number of genome DNA per genotype as template in PCR (default: 500)
    # --pcr_cycles: number of cycles in PCR (default: 25)
    # --output_filename: prefix of output .csv files (default: output)
    #
    # OUTPUTS
    # output_filename_EvoSimulation_Read_Number.csv: read number per genotype per sequencing time-point
    # output_filename_EvoSimulation_Mean_Fitness.csv: mean fitness per sequencing time-point
    # output_filename_EvoSimulation_Input_Log.csv: a record of all inputs
    # ------------------------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Simulated competitive pooled growth of a population of genotypes '
                                                 'with different fitness',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str,
                        help='a .csv file: 1st column is fitness of each genotype, '
                             '2nd column is initial cell number of each genotype at generation 0')
    parser.add_argument('-t', '--t_seq', nargs='*', type=int, help='sequenced time-points')
    parser.add_argument('-r', '--read_num_average_seq', nargs='*', type=int,
                        help='average number of reads per genotype per sequencing time-point')
    parser.add_argument('-n', '--noise_option', nargs='*', type=str,
                        default=['growth', 'bottleneck_transfer', 'DNA_extraction', 'PCR', 'sequencing'],
                        help='five types of possible noise: cell growth, bottleneck transfer, DNA extraction, PCR, '
                             'sequencing')
    parser.add_argument('-d', '--dna_copies', type=int, default=500,
                        help='average copy number of genome DNA per genotype')
    parser.add_argument('-p', '--pcr_cycles', type=int, default=25, help='number of cycles in PCR')
    parser.add_argument('-o', '--output_filename', type=str, default='output', help='prefix of output .csv files')

    args = parser.parse_args()
    csv_input = np.array(pd.read_csv(args.input, header=None), dtype=float)
    x = csv_input[:, 0]
    cell_num_ini = csv_input[:, 1].astype(int)
    t_seq = np.array(args.t_seq, dtype=int)
    read_num_average_seq = np.array(args.read_num_average_seq, dtype=int)
    noise_option = args.noise_option
    dna_copies = args.dna_copies
    pcr_cycles = args.pcr_cycles
    output_filename = args.output_filename

    delta_t = t_seq[1] - t_seq[0]
    lineages_num = np.max(cell_num_ini.shape)
    seq_num = np.max(t_seq.shape)
    evo_num = t_seq[-1]

    # Growth phase, with two possible noise involved: cell growth noise, bottleneck cell transfer noise
    if ('growth' in noise_option) and ('bottleneck_transfer' in noise_option):
        f1 = np.random.poisson
        f2 = np.random.poisson
    elif ('growth' in noise_option) and ('bottleneck_transfer' not in noise_option):
        f1 = np.random.poisson
        f2 = np.round
    elif ('growth' not in noise_option) and ('bottleneck_transfer' in noise_option):
        f1 = np.round
        f2 = np.random.poisson
    else:
        f1 = np.round
        f2 = np.round

    # -- Minimum cell number expected for a genotype, that is not growing and being lost to dilution
    cell_num_min_seq = np.zeros((lineages_num, evo_num + 1)).astype(int)
    cell_num_min_seq[:, 0] = cell_num_ini
    for j in range(1, t_seq[-1] + 1):
        cell_num_min_seq[:, j] = f2(np.true_divide(cell_num_ini, 2 ** (np.floor((j - 1) / delta_t) * delta_t)))

    cell_num_seq = np.zeros((lineages_num, evo_num + 1)).astype(int)
    cell_num_seq[:, 0] = cell_num_ini
    cell_num_saturated_seq = np.zeros((lineages_num, seq_num)).astype(int)
    cell_num_saturated_seq[:, 0] = cell_num_ini * (2 ** delta_t)
    x_mean = np.zeros(evo_num + 1)
    x_mean[0] = np.dot(cell_num_ini, x) / np.sum(cell_num_seq[:, 0])

    for j in range(1, evo_num + 1):
        pos = cell_num_seq[:, j - 1] != 0
        x_rela = np.max([np.true_divide(1 + x[pos], 1 + x_mean[j - 1]), np.zeros(np.sum(pos))], axis=0)
        cell_num_prev = cell_num_seq[pos, j - 1]
        cell_num_current = f1(np.multiply(2 * cell_num_prev, x_rela))
        cell_num_current = np.max([cell_num_current, cell_num_min_seq[pos, j]], axis=0)

        if np.mod(j, delta_t) == 0:
            ind = int(j / delta_t)
            cell_num_saturated_seq[pos, ind] = cell_num_current  # !!!
            cell_num_current = f2(cell_num_current / (2 ** delta_t))
        cell_num_seq[pos, j] = cell_num_current
        x_mean[j] = np.dot(x, cell_num_seq[:, j]) / np.sum(cell_num_seq[:, j])

    # After-growth phase, with three possible noise involved: DNA extraction noise, PCR noise, sequencing noise
    if ('DNA_extraction' in noise_option) and ('PCR' in noise_option) and ('sequencing' in noise_option):
        f3 = np.random.poisson
        f4 = np.random.poisson
        f5 = np.random.poisson
    elif ('DNA_extraction' in noise_option) and ('PCR' in noise_option) and ('sequencing' not in noise_option):
        f3 = np.random.poisson
        f4 = np.random.poisson
        f5 = np.round
    elif ('DNA_extraction' in noise_option) and ('PCR' not in noise_option) and ('sequencing' in noise_option):
        f3 = np.random.poisson
        f4 = np.round
        f5 = np.random.poisson
    elif ('DNA_extraction' not in noise_option) and ('PCR' in noise_option) and ('sequencing' in noise_option):
        f3 = np.round
        f4 = np.random.poisson
        f5 = np.random.poisson
    elif ('DNA_extraction' in noise_option) and ('PCR' not in noise_option) and ('sequencing' not in noise_option):
        f3 = np.random.poisson
        f4 = np.round
        f5 = np.round
    elif ('DNA_extraction' not in noise_option) and ('PCR' in noise_option) and ('sequencing' not in noise_option):
        f3 = np.round
        f4 = np.random.poisson
        f5 = np.round
    elif ('DNA_extraction' not in noise_option) and ('PCR' not in noise_option) and ('sequencing' in noise_option):
        f3 = np.round
        f4 = np.round
        f5 = np.random.poisson
    else:
        f3 = np.round
        f4 = np.round
        f5 = np.round

    dna_num_seq = f3(np.true_divide(cell_num_saturated_seq, np.sum(cell_num_saturated_seq, axis=0))
                     * dna_copies * lineages_num)

    pcr_num_seq = dna_num_seq
    for i in range(pcr_cycles):
        pcr_num_seq = f4(2 * pcr_num_seq)

    read_num_seq = np.multiply(np.true_divide(pcr_num_seq, np.sum(pcr_num_seq, axis=0)),
                               read_num_average_seq * lineages_num)
    read_num_seq = f5(read_num_seq)

    evo_simulator_output = {'Read_Number': read_num_seq,
                            'Mean_Fitness': x_mean[t_seq],
                            'Input_Log': {'Fitness': x,
                                          'Cell_Number_Initial': cell_num_ini,
                                          'Time_Points': t_seq,
                                          'Average_Read_Depth': read_num_average_seq,
                                          'Noise': noise_option,
                                          'gDNA_Copies': [dna_copies],
                                          'PCR_cycles': [pcr_cycles]}}

    tempt = pd.DataFrame(evo_simulator_output['Read_Number'])
    tempt.to_csv(output_filename + '_EvoSimulation_Read_Number.csv', index=False, header=False)
    tempt = pd.DataFrame(evo_simulator_output['Mean_Fitness'])
    tempt.to_csv(output_filename + '_EvoSimulation_Mean_Fitness.csv', index=False, header=False)
    tempt = list(itertools.zip_longest(*list(evo_simulator_output['Input_Log'].values())))
    with open(output_filename + '_EvoSimulation_Input_Log.csv', 'w') as f:
        w = csv.writer(f)
        w.writerow(evo_simulator_output['Input_Log'].keys())
        w.writerows(tempt)


if __name__ == "__main__":
    main()
