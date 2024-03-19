import pandas as pd
import matplotlib.pyplot as plt

df_diamond = pd.read_csv('diamond_ko_results_6.4_seed_10')
df_sourmash = pd.read_csv('sourmash_gather_6.4_seed_10_k_11')
df_ground_truth = pd.read_csv('ko_ground_truth_6.4_seed_10')

diamond_kos = set(df_diamond['ko_id'])
ground_truth_kos = set(df_ground_truth['ko_id'])
common_kos_diamond = diamond_kos.intersection(ground_truth_kos)

diamond_abund_columns = ['relative_abundance_by_nucleotides_covered', 'relative_abundance_by_num_reads']
ground_truth_abund_columns = ['abund_by_num_reads', 'abund_by_num_nts', 'abund_by_mean_cov']

# iterate over all combinations of diamond and ground truth abundance columns
for diamond_col, ground_truth_col in [(diamond_col, ground_truth_col) for diamond_col in diamond_abund_columns for ground_truth_col in ground_truth_abund_columns]:
    # create a dict: ko id to diamond abundance
    diamond_ko_to_abund = dict(zip(df_diamond['ko_id'], df_diamond[diamond_col]))

    # create a dict: ko id to ground truth abundance
    ground_truth_ko_to_abund = dict(zip(df_ground_truth['ko_id'], df_ground_truth[ground_truth_col]))

    # plot rel abund of diamond vs ground truth for the common kos
    diamond_abundances = [diamond_ko_to_abund[ko] for ko in common_kos_diamond]
    ground_truth_abundances = [ground_truth_ko_to_abund[ko] for ko in common_kos_diamond]

    plt.scatter(ground_truth_abundances, diamond_abundances)
    plt.xlabel(f'Ground Truth Abundance ({ground_truth_col})')
    plt.ylabel(f'Diamond Abundance ({diamond_col})')
    plt.title(f'Ground Truth vs Diamond Abundance ({ground_truth_col} vs {diamond_col})')
    plt.savefig(f'ground_truth_vs_diamond_abundance_{ground_truth_col}_vs_{diamond_col}.png')
    plt.clf()

sourmash_kos = set(df_sourmash['name'])
common_kos_sourmash = sourmash_kos.intersection(ground_truth_kos)

sourmash_abund_columns = ['f_unique_weighted', 'average_abund', 'median_abund', 'unique_intersect_bp', 'f_orig_query']

# iterate over all combinations of sourmash and ground truth abundance columns
for sourmash_col, ground_truth_col in [(sourmash_col, ground_truth_col) for sourmash_col in sourmash_abund_columns for ground_truth_col in ground_truth_abund_columns]:
    # create a dict: ko id to sourmash abundance
    sourmash_ko_to_abund = dict(zip(df_sourmash['name'], df_sourmash[sourmash_col]))

    # create a dict: ko id to ground truth abundance
    ground_truth_ko_to_abund = dict(zip(df_ground_truth['ko_id'], df_ground_truth[ground_truth_col]))

    # plot rel abund of sourmash vs ground truth for the common kos
    sourmash_abundances = [sourmash_ko_to_abund[ko] for ko in common_kos_sourmash]
    ground_truth_abundances = [ground_truth_ko_to_abund[ko] for ko in common_kos_sourmash]

    plt.scatter(ground_truth_abundances, sourmash_abundances)
    plt.xlabel(f'Ground Truth Abundance ({ground_truth_col})')
    plt.ylabel(f'Sourmash Abundance ({sourmash_col})')
    plt.title(f'Ground Truth vs Sourmash Abundance ({ground_truth_col} vs {sourmash_col})')
    plt.savefig(f'ground_truth_vs_sourmash_abundance_{ground_truth_col}_vs_{sourmash_col}.png')
    plt.clf()