"""Design oligopool for scalable gene assembly."""

import os
import sys
from os.path import join, exists
from typing import Optional, Union
import random

from jsonargparse import CLI
import pandas as pd
import numpy as np

# Cartographer imports
from data_classes import Enzyme, EnzymeTypes, LigationDataOpt, define_ligation_data, PrimerIterator, define_enzyme
from library_classes import Library
from junctions import optimize_junctions
from pricing import (
    annotate_pool_stats_with_primer_cost,
    cost_summary,
    write_cost_summary,
    TwistOligoPoolPricing,
    DEFAULT_TWIST_TABLE,
)
from Bio import SeqIO


def genes(
        input_seqs: str,
        njunctions: int,
        upstream_bbsite: str,
        downstream_bbsite: str,
        primers: str,
        output_dir: str = 'output',

        other_used_sites: Union[list[str], None] = None,

        # enzyme and ligation data arguments
        enzyme: Union[EnzymeTypes, Enzyme] = 'BsaI',
        illegal_dna_sequences: Union[tuple, None] = None,
        ligation_data: LigationDataOpt = 'T4_18h_37C',

        # oligo packaging arguments
        add_primers: bool = True,
        pad_oligos: bool = True,

        # optimization arguments
        nopt_steps: int = 1000,
        nopt_runs: Optional[int] = 5,
        opt_seeds: Optional[list[int]] = None,
        njobs: int = 1,

        # miscellaneous
        oligo_len: int = 350,
        min_size: int = 40,
        optimization: str = 'simulated_annealing',
        dev: bool = False,

        # pricing
        pricing_enabled: bool = True,
        twist_quote: bool = False,
) -> None:
    """
    Design library for pooled golden gate assembly.

    Args:
        input_seqs: File path to library sequences in fasta format.
        upstream_bbsite: GG site in str format if applicable.
        downstream_bbsite: GG site in str format if applicable.
        njunctions: Number of junctions used per pool - should include the upstream
            and downstream sites.
        other_used_sites: Other GG sites that are included in the assembly that are not
            upstream or downstream vector ligation sites.
        enzyme: Enzyme used to digest DNA fragments.
        ligation_data: Ligation data used to optimize the library.
        oligo_len: Oligo len used in oligopool.
        add_primers: If True, primer sequences are added to the ends of each oligo
        pad_oligos: If True, spacer DNA is added between primer sequences and restrction
            binding sites to make all oligos the same length to reduce amplification bias.
        min_size: Minimum fragment size allowed in fragmentation patterns.
        output_dir: Filepath to where optimization output is saved to.
        nopt_steps: Number of steps to run each optimization run.
        nopt_runs: Number of times each pool is optimized. Each new optimization run
            uses a different random seed.
        njobs: Number of CPUs to use for simultaneous optimization runs. Each CPU executes
            one run.

    """
    #pylint: disable=too-many-arguments, too-many-locals

    # if output directory doesn't exist, write it
    if not exists(output_dir):
        os.makedirs(output_dir)

    assembly_enzyme = define_enzyme(enzyme)
    #! Right now this only recognizes other TypeIIS enzymes
    illegal_dna_sequences = illegal_dna_sequences or ("")

    # format input sequences into Gene objects
    input_seqs_recs = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(input_seqs, 'fasta')]
    other_used_sites = other_used_sites or np.array([])
    primers = PrimerIterator(primers, assembly_enzyme)
    ligation_data = define_ligation_data(ligation_data, assembly_enzyme)

    # print input parameters
    print(f"Read {len(input_seqs_recs)} from {input_seqs}.\nSequence lengths are {min(map(lambda x: len(x[1]), input_seqs_recs))} - {max(map(lambda x: len(x[1]), input_seqs_recs))} bp")
    print(f"Using {enzyme.name} and {ligation_data.name} data to assemble library.")
    print(f"Maximum number of GG sites allowed per pool: {njunctions}")
    print(f"Upstream backbone site: {upstream_bbsite}")
    print(f"Downstream backbone site: {downstream_bbsite}")


    library = Library(
        genes=input_seqs_recs,
        primers=primers,
        oligo_len=oligo_len,
        enzyme=assembly_enzyme,
        upstream_bbsite=upstream_bbsite,
        downstream_bbsite=downstream_bbsite,
        other_used_sites=other_used_sites,
        illegal_dna_sequences=illegal_dna_sequences,
        njunctions=njunctions,
        min_size=min_size
    )
    
    # assign optimization seeds - use nopt_runs to get random_opt seeds
    random_seeds = None
    if nopt_runs is not None:
        rng = np.random.default_rng(seed=42)
        random_seeds = rng.integers(1000, size=nopt_runs)
    elif (nopt_runs is None) and (opt_seeds is not None):
        random_seeds = opt_seeds
    else:
        raise ValueError('Cannot provide values for both `nopt_runs` and `opt_seeds`.')

    # assign genes to pools
    use_plan = None
    better_plan = library.suggest_better_plan()
    if better_plan and sys.stdin.isatty():
        print("\n--- OMEGA Cost Optimization ---")
        current_cost = library.estimate_cost(library.plan_fragmentation())
        new_cost = library.estimate_cost(better_plan)
        savings = current_cost['total_cost'] - new_cost['total_cost']

        print(f"Standard fragmentation cost: ${current_cost['total_cost']:,.2f}")
        print(f"A cheaper plan exists: ${new_cost['total_cost']:,.2f} (Save ${savings:,.2f}!)")
        print("This plan increases fragmentation for some groups to drop into a cheaper Twist length tier.")

        try:
            ans = input("Would you like to use the cost-optimized plan? [Y/n]: ").strip().lower()
            if ans in ('', 'y', 'yes'):
                use_plan = better_plan
                print("Using cost-optimized plan.")
        except EOFError:
            pass

    library.optimize_pools(
        njunctions=njunctions,
        nopt_steps=nopt_steps,
        opt_seeds=random_seeds,
        njobs=njobs,
        ligation_data=ligation_data.data,
        optimization=optimization,
        planned_groups=use_plan
    )

    optimized_library = library.package_library(add_primers=add_primers, pad_oligo=pad_oligos)
    optimized_library.to_csv(os.path.join(output_dir, 'optimization_results.csv'))
    print(f"Finished optimization. Saved to {os.path.join(output_dir, 'optimization_results.csv')}")

    oligopool = library.package_oligos(add_primers=add_primers, pad_oligo=pad_oligos)
    oligopool.to_csv(os.path.join(output_dir, 'oligo_order.csv'))
    print(f"Designed {len(oligopool)} oligos. Saved to {os.path.join(output_dir, 'oligo_order.csv')}")

    pool_stats = pd.DataFrame.from_dict(
        [{
            'pool':p.name,
            'fidelity':fidelity,
            'min_gene':p.min_gene_fidelity,
            'min_site':p.min_site_fidelity,
            'n_genes':len(p.genes),
            'n_sites':(p.nfrags-1)*len(p.genes) + 2 + len(p.other_used_sites),
            'seed':seed,
            'enzyme':p.enzyme.name,
            'pfwd_name':p.fprimer.name,
            'pfwd_sequence':p.fprimer.sequence,
            'prev_name':p.rprimer.name,
            'prev_sequence':p.rprimer.sequence
        } for p, fidelity, seed in library.optimized_pools]
    )
    if pricing_enabled:
        try:
            pool_stats = annotate_pool_stats_with_primer_cost(pool_stats)
        except Exception as exc:
            print(f"Primer pricing on pool_stats skipped: {exc}")

    pool_stats.to_csv(join(output_dir, 'pool_stats.csv'))
    print(f"Pool-level statistics saved to {join(output_dir, 'pool_stats.csv')}")

    with open(join(output_dir, 'experiment_details.txt'), 'w') as f:
        f.write(ligation_data.experiment_information())

    if pricing_enabled:
        try:
            summary = cost_summary(
                oligopool,
                twist_quote=twist_quote,
                pool_stats_df=pool_stats,
            )
            out = write_cost_summary(output_dir, summary)
            offline = summary['offline_pool_price_usd']
            print(f"Twist offline list price: ${offline:,.2f} "
                  f"(tier {summary['offline_tier']}, {summary['offline_length_bin']})")
            
            # Efficiency alert
            eff = TwistOligoPoolPricing.from_csv(DEFAULT_TWIST_TABLE).get_tier_efficiency(len(oligopool))
            if eff and eff.get('over_prev_tier', 0) > 0 and eff['over_prev_tier'] < (eff['tier_max'] - eff['tier_min']) * 0.05:
                print(f"\n[ADVISORY] You are only {eff['over_prev_tier']} oligos over the previous Twist pricing tier boundary ({eff['prev_tier_max']}).")
                print("Consider dropping a few non-essential genes or increasing junction counts to consolidate pools and save money.")
            if twist_quote and 'live_pool_subtotal_usd' in summary:
                print(f"Twist live quote subtotal: ${summary['live_pool_subtotal_usd']:,.2f}  "
                      f"total: ${summary['live_total_price_usd']:,.2f}")
            if 'primers_total_usd' in summary:
                print(f"IDT primer pairs: ${summary['primers_total_usd']:,.2f} "
                      f"({summary['n_pools']} pools at "
                      f"${summary['primers_per_pool_avg_usd']:,.2f}/pool)")
            if 'wetlab_total_steps' in summary:
                print(f"Wet-lab steps: {summary['wetlab_total_steps']} total "
                      f"({summary['wetlab_pcrs']} PCRs + "
                      f"{summary['wetlab_pcr_cleanups']} PCR cleanups + "
                      f"{summary['wetlab_quants']} PicoGreen quants + "
                      f"{summary['wetlab_assembly_reactions']} GG assemblies + "
                      f"{summary['wetlab_final_cleanups']} final cleanup + "
                      f"{summary['wetlab_transformations']} transformation)")
            print(f"Cost summary saved to {out}")
        except Exception as exc:
            print(f"Pricing skipped: {exc}")

    # if dev, save extra data
    if dev:
        # make trajectory save dir
        trajectory_dir = join(output_dir, 'opt_trajectories')
        if not exists(trajectory_dir):
            os.mkdir(trajectory_dir)
        for pool in library.all_opt_runs:
            # print(pool[0].opt_trajectory)
            np.save(join(trajectory_dir, f'pool_{pool[0].name}_seed-{pool[2]}.npy'), np.array(pool[0].opt_trajectory))

def junctions(
        set_size: int,
        fixed_sites: Optional[list[str]],
        excluded_sites: Optional[list[str]],
        enzyme: Union[EnzymeTypes, Enzyme] = 'BsaI',
        ligation_data: LigationDataOpt = 'T4_01h_25C',
        output_dir: str = 'output',

        nopt_steps: int = 5000,
        nopt_runs: Optional[int] = None,
        opt_seeds: Optional[list[int]] = None,
        njobs: int = 1,

        dev: bool = False
) -> None:
    """
    TODO: Make this a better docstring
    Design library for pooled golden gate assembly.

    Args:
        input_seqs: File path to library sequences in fasta format.
        upstream_bbsite: GG site in str format if applicable.
        downstream_bbsite: GG site in str format if applicable.
        njunctions: Number of junctions used per pool - should include the upstream
            and downstream sites.
        other_used_sites: Other GG sites that are included in the assembly that are not
            upstream or downstream vector ligation sites.
        enzyme: Enzyme used to digest DNA fragments.
        ligation_data: Ligation data used to optimize the library.
        oligo_len: Oligo len used in oligopool.
        add_primers: If True, primer sequences are added to the ends of each oligo
        pad_oligos: If True, spacer DNA is added between primer sequences and restrction
            binding sites to make all oligos the same length to reduce amplification bias.
        min_size: Minimum fragment size allowed in fragmentation patterns.
        output_dir: Filepath to where optimization output is saved to.
        nopt_steps: Number of steps to run each optimization run.
        nopt_runs: Number of times each pool is optimized. Each new optimization run
            uses a different random seed.
        njobs: Number of CPUs to use for simultaneous optimization runs. Each CPU executes
            one run.

    """
    #pylint: disable=too-many-arguments, too-many-locals

    # if output directory doesn't exist, write it
    if not exists(output_dir):
        os.makedirs(output_dir)

    assembly_enzyme = define_enzyme(enzyme)
    ligation_data = define_ligation_data(ligation_data, assembly_enzyme)

    # assign optimization seeds - use nopt_runs to get random_opt seeds
    random_seeds = None
    if nopt_runs is not None:
        rng = np.random.default_rng(seed=42)
        random_seeds = rng.integers(1000, size=nopt_runs)
    elif (nopt_runs is None) and (opt_seeds is not None):
        random_seeds = opt_seeds
    else:
        raise ValueError('Cannot provide values for both `nopt_runs` and `opt_seeds`.')

    # set random seed
    random.seed(42)

    # optimize pools
    optimization_results = optimize_junctions(
        name=0,
        enzyme=assembly_enzyme,
        set_size=set_size,
        fixed_sites=fixed_sites,
        excluded_sites=excluded_sites,
        ligation_data=ligation_data.data,
        opt_seeds=random_seeds,
        nopt_steps=nopt_steps,
        njobs=njobs
    )

    # save best optimized set
    optimized_set, _, random_seed = optimization_results[-1]

    opt_stats = optimized_set.package_sites() | {'random_seed':random_seed}
    pd.DataFrame.from_dict([opt_stats]).to_csv(join(output_dir, 'optimized_junction_set.csv'))

    if dev:
        # make trajectory save dir
        trajectory_dir = join(output_dir, 'opt_trajectories')
        if not exists(trajectory_dir):
            os.mkdir(trajectory_dir)
        for jset, _, seed in optimization_results:
            np.save(join(trajectory_dir, f'set_{jset.name}_seed-{seed}.npy'), np.array(optimized_set.opt_trajectory))

def costs(
        output_dir: str,
        twist_quote: bool = False,
) -> None:
    """Re-price an existing OMEGA output directory.

    Reads ``oligo_order.csv`` from ``output_dir`` and writes (or overwrites)
    ``cost_summary.csv``. Useful for repricing past runs without re-running
    optimization, or for adding ``--twist-quote`` after the fact.
    """
    oligo_path = join(output_dir, 'oligo_order.csv')
    if not exists(oligo_path):
        raise FileNotFoundError(f"No oligo_order.csv in {output_dir}")
    oligopool = pd.read_csv(oligo_path)
    pool_stats_df = None
    pool_stats_path = join(output_dir, 'pool_stats.csv')
    if exists(pool_stats_path):
        pool_stats_df = pd.read_csv(pool_stats_path)
    summary = cost_summary(
        oligopool, twist_quote=twist_quote, pool_stats_df=pool_stats_df,
    )
    out = write_cost_summary(output_dir, summary)
    offline = summary['offline_pool_price_usd']
    print(f"Twist offline list price: ${offline:,.2f} "
          f"(tier {summary['offline_tier']}, {summary['offline_length_bin']})")
    if twist_quote and 'live_pool_subtotal_usd' in summary:
        print(f"Twist live quote subtotal: ${summary['live_pool_subtotal_usd']:,.2f}  "
              f"total: ${summary['live_total_price_usd']:,.2f}")
    if 'primers_total_usd' in summary:
        print(f"IDT primer pairs: ${summary['primers_total_usd']:,.2f} "
              f"({summary['n_pools']} pools at "
              f"${summary['primers_per_pool_avg_usd']:,.2f}/pool)")
    if 'wetlab_total_steps' in summary:
        print(f"Wet-lab steps: {summary['wetlab_total_steps']} total "
              f"({summary['wetlab_pcrs']} PCRs + "
              f"{summary['wetlab_pcr_cleanups']} PCR cleanups + "
              f"{summary['wetlab_quants']} PicoGreen quants + "
              f"{summary['wetlab_assembly_reactions']} GG assemblies + "
              f"{summary['wetlab_final_cleanups']} final cleanup + "
              f"{summary['wetlab_transformations']} transformation)")
    print(f"Cost summary saved to {out}")


if __name__ == "__main__":
    CLI([genes, junctions, costs], as_positional=False)
