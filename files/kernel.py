"""
Spectrum Class.

This module allows the user to
apply a Spectrum Kernel to Kmers.

Classes
-------
Spectrum
"""

from sequence import Sequence


class Spectrum:
    """A class to represent Spectrum Kernel for Sequence comparison."""

    def __init__(self, size: int) -> None:
        """Construct all attributes for Spectrum."""
        self.size = size

    @property
    def size(self) -> int:
        """Size of kmer."""
        return self._size

    @size.setter
    def size(self, size: int) -> None:
        self._size = size

    def _union(self, kmers1: list[str], kmers2: list[str]) -> list[str]:
        """Find union of kmers."""
        kmerSet: set[str] = set(kmers1)
        union: list[str] = list(kmerSet.union(kmers2))
        return union

    def _featVec(self, seq: Sequence, union: list[str]) -> None:
        """Compute feature vector."""
        count: int = 0
        v: list[int] = list()
        for i in range(len(union)):
            for j in range(len(seq.kmers)):
                if seq.kmers[j] == union[i]:
                    count += 1
            v.append(count)
            count = 0
        seq.features = v

    def dotProduct(self, v1: list[int], v2: list[int]) -> int:
        """Compute dot product between feature vectors."""
        dot: int = 0
        for i in range(len(v1)):
            dot += v1[i] * v2[i]
        return dot
