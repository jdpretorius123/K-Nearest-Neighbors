"""
Write KNN Output.

This module allows the user to setup
and run the KNN algorithm.

Functions
---------
write(argv: list[str]) -> None:
    Writes classification results.
"""

import os
from typing import TextIO
from file import FastaFile
from sequence import Sequence
from kernel import Spectrum
from knn import Knn


def write(argv: list[str]) -> None:
    """Write clustering results."""
    trainPath: str = argv[0]
    trainFile: FastaFile = FastaFile(trainPath)
    train: list[Sequence] = trainFile.getSequences()

    testPath: str = argv[1]
    testFile: FastaFile = FastaFile(testPath)
    test: list[Sequence] = testFile.getSequences()

    size: int = int(argv[2])
    kernel: Spectrum = Spectrum(size)

    k: int = int(argv[3])
    knn: Knn = Knn(train, test, kernel, k)

    outfile: str = argv[4]
    knn.write(outfile)
    

def runKnn() -> None:
    """Run K-Means algorithm."""
    trainIn: list[str] = ["KNN/train10.fasta", "KNN/train30.fasta", "KNN/train100.fasta"]
    testIn: str = "KNN/test.fasta"
    # trainIn: list[str] = ["KNN/train-small.fasta"]
    # testIn: str = "KNN/test-small.fasta"
    sizes: list[str] = ["2", "4", "6", "8"]
    neighbors: list[str] = ["1", "3", "5", "7"]
    outfile: str = "KNN/knn.txt"

    if os.path.isfile(outfile):
        os.remove(outfile)

    columns: str = "SampleSize\tKmerSize\tNeighbors\tAccuracy\n"
    out: TextIO = open(outfile, "a")
    out.write(columns)
    out.close()

    argv: list[str] = list()
    for file in trainIn:
        for size in sizes:
            for n in neighbors:
                argv = [file, testIn, size, n, outfile]
                write(argv)
