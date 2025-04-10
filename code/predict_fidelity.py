
from typing import Iterable
from itertools import chain
from Bio.Seq import Seq

import numpy as np
import pandas as pd


def correct_ligations(site: str, wc_site: str, data_matrix: pd.DataFrame) -> int:
    """Retrieves number of correct ligations observed for a GG site."""

    return data_matrix.loc[site, wc_site]


def total_ligations(site: str, sites: Iterable[str], data_matrix: pd.DataFrame) -> int:
    """Retrieves total number of ligations observed for a set of GG sites
    in reference to a query golden gate site.
    """
    wc_sites = [str(Seq(s).reverse_complement()) for s in sites]

    return data_matrix.loc[site, np.concatenate((sites, wc_sites))].sum()


def site_probability(site: str, sites: Iterable[str], data_matrix: pd.DataFrame) -> float:
    """Calculate probability of site orthogonality in reference to proposed set of sites."""

    wc_site = str(Seq(site).reverse_complement())

    # get correct ligations for the site and it's watson crick pair
    site_correct_ligations = correct_ligations(site=site, wc_site=wc_site, data_matrix=data_matrix)
    wc_correct_ligations = correct_ligations(site=wc_site, wc_site=site, data_matrix=data_matrix)

    # total correct ligations for both the site and WC pair
    correct_total = site_correct_ligations + wc_correct_ligations

    # all liagation events
    total = total_ligations(site, sites, data_matrix) + total_ligations(wc_site, sites, data_matrix)

    return correct_total / total


def predict_fidelity(sites: Iterable[str], data: pd.DataFrame) -> float:
    """Calculate the predicted fidelity for a set of GG sites."""

    return np.prod([site_probability(site, sites, data) for site in sites])

def predict_minimum_site(sites: Iterable[list[str]], data: pd.DataFrame) -> float:
    """Get least orthogonal site from set."""

    return min(site_probability(site, sites, data) for site in sites)

def geneset_fidelity(gene_sites: Iterable[list[str]], data: pd.DataFrame) -> float:
    """Get fidelities for every gene in a set of genes (typically a pool).
    
    Args:
        gene_sites (Iterable): an iterable containing a list-like container of gene sets.

    Returns:
        A list of fidelities whose index matches the index of the gene. 
    """

    pool_sites = list(chain(*gene_sites))

    return [
        np.prod([site_probability(site, pool_sites, data) for site in sets]) for sets in gene_sites
    ]


def predict_minimum(gene_sites: Iterable[np.ndarray], data: pd.DataFrame) -> float:
    """Return minimum fidelity of things in a pool.
    
    Args:
        pool_sites (Iterable[np.ndarray]): iterable of ndarray's containing gg sites.
                                           each array corresponds to a gene.

    Returns:
        Lowest fidelity for a gene in a pool.
    
    """
    return min(geneset_fidelity(gene_sites, data))


def predict_average(gene_sites: Iterable[np.ndarray], data: pd.DataFrame) -> float:
    """Return average fidelity for genes in a pool.

    Args:
        pool_sites (Iterable[np.ndarray]): iterable of ndarray's containing gg sites.
                                           each array corresponds to a gene.

    Returns:
        Average fidelity for a gene in a pool.
    """

    return sum(geneset_fidelity(gene_sites, data)) / len(gene_sites)
