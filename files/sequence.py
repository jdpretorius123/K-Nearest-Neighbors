"""
Sequence Class.

This module allows the user to store one
sequence from a fasta file in a Sequence instance.

Classes
-------
Sequence
"""

from __future__ import annotations


class Sequence:
    """A class to represent a sequence."""

    def __init__(self, seq: str, id: str) -> None:
        """Construct all attributes for Sequence."""
        self.seqStr = seq
        self.seqLst: list[str] = list(seq)
        self.id = id
        self.predId: str = ""
        self.kmers: list[str] = list()
        self.features: list[int] = list()
        self.neighbors: list[tuple[Sequence, float]] = list()

    @property
    def seqStr(self) -> str:
        """Sequence as string."""
        return self._seqStr

    @seqStr.setter
    def seqStr(self, seq: str) -> None:
        if isinstance(seq, str):
            self._seqStr = seq
        else:
            raise ValueError('"seq" must be a str')

    @property
    def seqLst(self) -> list[str]:
        """Sequence as list."""
        return self._seqLst

    @seqLst.setter
    def seqLst(self, seq: list[str]) -> None:
        self._seqLst = seq

    @property
    def id(self) -> str:
        """Sequence id."""
        return self._id

    @id.setter
    def id(self, id: str) -> None:
        self._id = id

    @property
    def predId(self) -> str:
        """Sequence predicted id."""
        return self._predId

    @predId.setter
    def predId(self, predId: str) -> None:
        self._predId = predId

    @property
    def kmers(self) -> list[str]:
        """Sequence kmers."""
        return self._kmers

    @kmers.setter
    def kmers(self, kmers: list[str]) -> None:
        self._kmers = kmers

    @property
    def features(self) -> list[int]:
        """Feature vector of Sequence."""
        return self._features

    @features.setter
    def features(self, features: list[int]) -> None:
        self._features = features

    @property
    def neighbors(self) -> list[tuple[Sequence, float]]:
        """Neighbors of Sequence."""
        return self._neighbors

    @neighbors.setter
    def neighbors(self, neighbors: list[tuple[Sequence, float]]) -> None:
        self._neighbors = neighbors

    def _buildKmers(self, size: int) -> None:
        """Build self.kmers."""
        kmers: list[str] = list()
        length: int = self.getLength()
        stop: int = length - (size - 1)
        for i in range(stop):
            end: int = i + size
            kmer: str = self.seqStr[i:end]
            kmers.append(kmer)
        self.kmers = kmers

    def toList(self) -> None:
        """Convert Sequence to list."""
        self.seqLst = list(self.seqStr)

    def reverse(self) -> None:
        """Reverse Sequence as list."""
        self.seqLst.reverse()

    def getBase(self, pos: int) -> str:
        """Return base pair in Sequence."""
        length: int = self.getLength()
        if pos >= length:
            raise ValueError(f"{pos} greater than length {length}")
        base: str = self.seqStr[pos]
        return base

    def getLength(self) -> int:
        """Return Sequence length."""
        length: int = len(self.seqStr)
        return length

    def buildKmers(self, size: int) -> None:
        """Build kmers."""
        self._buildKmers(size)

    def print(self, rep: str) -> None:
        """Print Sequence as string, list, or kmers."""
        match rep:
            case "string":
                print(self.seqStr)
            case "list":
                print(" ".join(self.seqLst))
            case "kmers":
                print(" ".join(self.kmers))
