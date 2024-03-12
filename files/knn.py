"""
Knn Class.

This module allows the user to run
the K Nearest Neighbors algorithm.

Classes
-------
Knn
"""

from typing import TextIO
from kernel import Spectrum
from sequence import Sequence

from file import FastaFile

NEIGHBORS = list[tuple[Sequence, float]]


class Knn:
    """A class to represent KNN."""

    def __init__(
        self,
        train: list[Sequence],
        test: list[Sequence],
        kernel: Spectrum,
        k: int,
    ) -> None:
        """Construct all attributes for Knn."""
        self.train = train
        self.test = test
        self.kernel = kernel
        self.k = k

    @property
    def train(self) -> list[Sequence]:
        """All Sequences for training."""
        return self._train

    @train.setter
    def train(self, train: list[Sequence]) -> None:
        self._train = train

    @property
    def test(self) -> list[Sequence]:
        """All Sequences for testing."""
        return self._test

    @test.setter
    def test(self, test: list[Sequence]) -> None:
        self._test = test

    @property
    def kernel(self) -> Spectrum:
        """Kernel for KNN algorithm."""
        return self._kernel

    @kernel.setter
    def kernel(self, kernel: Spectrum) -> None:
        self._kernel = kernel

    @property
    def k(self) -> int:
        """Number of neighbors."""
        return self._k

    @k.setter
    def k(self, k: int) -> None:
        self._k = k

    def _createKmers(self, seqs: list[Sequence]) -> None:
        """Create Sequence kmers."""
        size: int = self.kernel.size
        for seq in seqs:
            seq.buildKmers(size)

    def _create(self):
        """Initialize training and testing sets."""
        self._createKmers(self.train)
        self._createKmers(self.test)

    def _computeFeatVecs(
        self, seq1: Sequence, seq2: Sequence, union: list[str]
    ) -> None:
        """Compute Sequence feature vectors."""
        self.kernel._featVec(seq1, union)
        self.kernel._featVec(seq2, union)

    def _computeDot(self, seq: Sequence, seqs: list[Sequence]) -> None:
        """Compute all dot products between seq and seqs."""
        neighbors: NEIGHBORS = list()
        for seq2 in seqs:
            union: list = self.kernel._union(seq.kmers, seq2.kmers)
            self._computeFeatVecs(seq, seq2, union)
            dot: float = self.kernel.dotProduct(seq.features, seq2.features)
            neighbors.append((seq2, dot))
        seq.neighbors = neighbors

    def _sortByDot(self, t) -> float:
        """Return neighbor's dot product."""
        return t[1]

    def _sort(self, seq: Sequence) -> NEIGHBORS:
        """Sort k nearest neighbors of seq."""
        neighbors: NEIGHBORS = seq.neighbors[:]
        neighbors.sort(reverse=True, key=self._sortByDot)
        neighbors = neighbors[: self.k]
        # for n in neighbors:
        #     print(n[0].id, n[1], sep = " ", end = " ")
        return neighbors

    def _predId(self, neighbors: NEIGHBORS) -> str:
        """Return majority id among neighbors."""
        exon: int = 0
        intron: int = 0
        for n in neighbors:
            if n[0].id == "exon":
                exon += 1
            else:
                intron += 1
        if exon >= intron:
            return "exon"
        return "intron"

    def _evaluate(self, seq: Sequence) -> bool:
        """Evaluate accuracy of Sequence predId."""
        # print(f"Pred: {seq.predId} Real: {seq.id}\n")
        return seq.predId == seq.id

    def _execute(self) -> float:
        """Excecute KNN algorithm."""
        total: int = len(self.test)
        correct: int = 0
        self._create()
        for seq in self.test:
            self._computeDot(seq, self.train)
            neighbors: NEIGHBORS = self._sort(seq)
            seq.predId = self._predId(neighbors)
            if self._evaluate(seq):
                correct += 1
        accuracy: float = correct / total
        return accuracy

    def write(self, outfile: str) -> None:
        """Write KNN algorithm results to outfile."""
        out: TextIO = open(outfile, "a")
        n: int = len(self.train) // 2
        kmer: int = self.kernel.size
        acc: float = self._execute()
        acc = round(acc, 3)
        line: str = f"{n}\t\t\t{kmer}\t\t\t{self.k}\t\t\t{acc}\n"
        out.write(line)
        out.close()


# train: str = "KNN/train10.fasta"
# trainF: FastaFile = FastaFile(train)
# trainS: list[Sequence] = trainF.getSequences()

# test: str = "KNN/test.fasta"
# testF: FastaFile = FastaFile(test)
# testS: list[Sequence] = testF.getSequences()

# size: int = 8
# kernel: Spectrum = Spectrum(size)

# k: int = 1
# knn: Knn = Knn(trainS, testS, kernel, k)
# knn.write("KNN/knn.txt")
# acc: float = knn._execute()
# print(f"Acc: {acc}")
