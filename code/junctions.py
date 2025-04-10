"""Code to optimize set of junctions without library sequences. Used to design junction sets tested in Fig 2."""

from typing import Any, Optional
from itertools import product
import random

import numpy as np
from tqdm import tqdm

from data_classes import Enzyme
from predict_fidelity import predict_fidelity

from joblib import Parallel, delayed


def optimize_junctions(
    name: Any, enzyme: Enzyme, set_size: int, ligation_data: np.ndarray, opt_seeds: list[int],
    fixed_sites: Optional[list[str]] = None, excluded_sites: Optional[list[str]] = None,
    nopt_steps: int = 1000, njobs: int = 1
):

    opt_results = Parallel(n_jobs=njobs)(delayed(single_optimization)(
        name=name,
        enzyme=enzyme,
        set_size=set_size,
        ligation_data=ligation_data,
        fixed_sites=fixed_sites,
        excluded_sites=excluded_sites,
        nopt_steps=nopt_steps,
        random_seed=rand_seed
    ) for rand_seed in opt_seeds)

    opt_results.sort(key=lambda x: x[1])

    return opt_results

def single_optimization(
    name: Any, enzyme: Enzyme, set_size: int, ligation_data: np.ndarray,
    fixed_sites: Optional[list[str]] = None, excluded_sites: Optional[list[str]] = None,
    nopt_steps: int = 1000, random_seed: int = 42
):
    jset = JunctionSet(
        name=name,
        enzyme=enzyme,
        set_size=set_size,
        ligation_data=ligation_data,
        fixed_sites=fixed_sites,
        excluded_sites=excluded_sites
    )

    fidelity = jset.optimize(nopt_steps=nopt_steps, random_seed=random_seed, disable_progress=True)

    return jset, fidelity, random_seed

class JunctionSet:

    """Class object to design sets of junctions."""

    def __init__(
        self,
        name: Any,
        enzyme: Enzyme,
        set_size: int,
        ligation_data,
        fixed_sites: Optional[list[str]] = None,
        excluded_sites: Optional[list[str]] = None
    ):
        
        self.name = name
        self.enzyme = enzyme
        self.set_size = set_size
        self.fixed_sites = fixed_sites or []
        self.excluded_sites = excluded_sites or []
        self.ligation_data = ligation_data

        # site options to pick from
        self.site_options = self.__site_options()
        # assign starting set of sites
        self.current_sites: np.ndarray = self.__assign_start_sites()

        self.optimized_sites = None
        self.opt_trajectory = None
        self.optimized_fidelity = None

    def shuffle_site(self, sites: np.ndarray) -> np.ndarray:
        """Replace a random site with a new site."""

        candidate_sites = sites.copy()
        site_idx = random.choice(range(len(candidate_sites)))
        
        options = [s for s in self.site_options if s != candidate_sites[site_idx]]
        candidate = random.choice(options)
        candidate_sites[site_idx] = candidate

        return candidate_sites

    def optimize(self, nopt_steps: int, random_seed: int, disable_progress: False) -> float:
        random.seed(random_seed)
        start_temp: float = 0.005
        end_temp: float = 0.00001

        current_state: np.ndarray = self.current_sites.copy()
        best_state: np.ndarray = self.current_sites.copy()

        start_fidelity = predict_fidelity(
            np.hstack([
                current_state,
                np.array(self.fixed_sites)
            ]),
            self.ligation_data
        )

        opt_trajectory = [(start_fidelity, start_fidelity)]  # best, current

        for temp in tqdm(np.linspace(start_temp, end_temp, nopt_steps), ncols=100, total=nopt_steps, disable=disable_progress, leave=True):
            candidates: np.array = self.shuffle_site(current_state)
            candidate_fidelity = predict_fidelity(candidates.tolist()+self.fixed_sites, self.ligation_data)

            energy_diff = candidate_fidelity - opt_trajectory[-1][1]
            thing_to_beat = random.random()  # threshold that determines if candidates are accepted
            contender = np.exp(min([0, energy_diff / temp]))  # energy of new sites adjusted by current temperature

            if contender > thing_to_beat:

                # if we've gotten too far away from the best, then reset to the best
                if (opt_trajectory[-1][0] - candidate_fidelity) > 0.01:
                    current_state = best_state.copy()

                    # append the fidelity - best and current are now best, candidate is what triggered reset to best
                    opt_trajectory.append((opt_trajectory[-1][0], opt_trajectory[-1][0]))

                else:
                    current_state = candidates.copy()

                    # if new fidelity is better than best fidelity, update best fidelity
                    if candidate_fidelity > opt_trajectory[-1][0]:
                        best_state = current_state.copy()
                        opt_trajectory.append((candidate_fidelity, candidate_fidelity))
                    
                    else:
                        # don't update best
                        opt_trajectory.append((opt_trajectory[-1][0], candidate_fidelity))

            # if change is rejected, do nothing except update trajectory with old fidelity
            else:
                opt_trajectory.append((opt_trajectory[-1][0], opt_trajectory[-1][1]))

        self.optimized_sites = best_state
        self.opt_trajectory = opt_trajectory
        self.optimized_fidelity = predict_fidelity(self.optimized_sites.tolist()+self.fixed_sites, self.ligation_data)

        return self.optimized_fidelity
    
    def package_sites(self) -> dict:
        opt_sites = self.fixed_sites + self.optimized_sites.tolist()
        opt_sites = {f'ggsite_{i}':ggsite for i, ggsite in enumerate(opt_sites)}
        return {'pool_id':self.name, 'fidelity':self.optimized_fidelity} | opt_sites

    def __site_options(self) -> np.ndarray:
        """Retrieve site options for junction set."""

        palindromic_sites = ['GCGC', 'CGCG', 'ATAT', 'TATA', 'GGCC', 'GGCC', 'AATT', 'TTAA', 'TGCA', 'AGCT','TCGA', 'ACGT', 'GATC', 'GTAC', 'CATG', 'CTAG']
        illegal_sites = self.fixed_sites + palindromic_sites
        options = np.array([
            s for s in map("".join, product(*[list("ATGC")]*4)) if s not in illegal_sites
        ])
        return options


    def __assign_start_sites(self) -> np.ndarray:
        """Assign a random set of start sites."""
        nsites = self.set_size - len(self.fixed_sites)
        return np.array(random.sample(self.site_options.tolist(), nsites) + self.fixed_sites)