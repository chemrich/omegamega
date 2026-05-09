"""Random helper functions."""

from functools import lru_cache
from itertools import chain
from typing import Iterable, Any
from random import shuffle

import regex as re
import numpy as np
import pandas as pd

from Bio.Seq import Seq
from Bio import SeqIO


_RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def _rc(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def unique_orthogonal(sites) -> bool:
    """Checks if GG sites are orthogonal with the Watson Crick pairs."""

    wc_sites = [_rc(s) for s in sites]

    return len(set(sites.tolist() + wc_sites)) == len(sites)*2


@lru_cache(maxsize=128)
def _compile_dna_pattern(seq_elements: tuple[str, ...], reverse_complement: bool):
    """Compile a single regex matching any of seq_elements (and their RCs)."""
    if reverse_complement:
        elements = list(chain.from_iterable((s, _rc(s)) for s in seq_elements))
    else:
        elements = list(seq_elements)
    return re.compile("|".join(elements), flags=re.I)


def dna_contains_seq(dna: str, *seq_elements: str, reverse_complement: bool = True) -> bool:
    """Checks if dna sequence contains defined sequences.
    True if sequence is present. If reverse_complement is True, then
    searches for reverse complement as well.
    """
    pattern = _compile_dna_pattern(tuple(seq_elements), reverse_complement)
    return pattern.search(dna) is not None

def count_sequence_element(dna: str, seq_element: str) -> int:
    """Returns number of times a DNA element is found in a DNA
    sequence. Considers all overlapping regions of DNA."""

    matches = re.findall(
        f'{seq_element}|{str(Seq(seq_element).reverse_complement())}',
        dna,
        flags=re.IGNORECASE,
        overlapped=True
    )

    return len(matches)

def dynamic_chunker(iterable: Iterable[Any], chunk_sizes: list[int]) -> Iterable[list[Any]]:
    """Break up list into variable-sized chunks specified in `chunk_sizes`"""

    it = iter(iterable)

    for chunk in chunk_sizes:
        yield [next(it) for i in range(chunk)]

_DNA_BASES = np.array(list("ATCG"))
_DNA_PROBS = np.array([0.3, 0.3, 0.2, 0.2])


def random_dna(size: int) -> str:
    """Generate random string of DNA with 40% GC content."""
    if size <= 0:
        return ""
    bases = np.random.choice(_DNA_BASES, size=size, p=_DNA_PROBS)
    return bases.view((str, size)).item()

def flatten(iterable: list) -> list:
    """Remove one level of a nested list."""

    return list(chain(*iterable))

def index_array(iterable: Iterable[int], length: int) -> np.ndarray:
    """Return array with 1's at the indexes included in iterable."""

    return np.array([1 if i in iterable else 0 for i in range(length)])

def read_fasta(file: str) -> pd.DataFrame:
    """Returns dataframe of sequences from fasta file. Preserves
    all information included in fasta file.
    
    Args:
        file (str): relative file path

    Returns:
        pd.DataFrame with `sequence` column holding sequences. Index is automatically
        0 to n-sequences, meaning that this format does not accept custom labels.
    """

    records = list(SeqIO.parse(file, format='fasta'))
    records = [
        pd.DataFrame({
            'sequence':[str(r.seq)], 'fasta_id':[r.id], 'fasta_description':[r.description]
        }) for r in records
    ]

    return pd.concat(records, ignore_index=True)

def overlap_split(iterable, overlap=1) -> list:
    """Returns a nested list of split input iterable with n elements
    overlapping between each sublist. No unique elements except
    first and last element of input iterable.
    
    Args:
        iterable: the iterable to split
        overlap: number of elements to repeat between each sublist
    
    Returns:
        nested list
    """

    return [[iterable[i:i+1+overlap]] for i in range(len(iterable)-overlap)]

def shuffler(iterable):
    iterable = iterable.copy()
    shuffle(iterable)

    return iterable
