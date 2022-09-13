"""Parallel (multi-threaded) map function for python. 
 
Uses multiprocessing.Pool with error-resistant importing. There are two map
functions:
 
1) pmap(function, iterable) -> rapid fork-based multi-threaded map function.
 
2) low_memory_pmap(function, iterable) -> a more memory-efficient version
    intended for function calls that are individually long & memory-intensive.
 
"""
import os
from warnings import warn
from pickle import PicklingError
from time import sleep 
try:
    import multiprocessing
except ImportError:
    print("Cannot import 'multiprocessing' module. Parallelization not possible.")
    pmap = map
    low_memory_pmap = map
    larger_iter_pmap = map
    CPUs = 1
finally:
    CPUs = multiprocessing.cpu_count()
    CHUNKS = 50*CPUs
    def pmap(func, Iter, processes=CPUs):
        with multiprocessing.Pool(processes=processes) as P:
            return P.map(func, Iter)
 
    def low_memory_pmap(func, Iter, processes=int(round(CPUs/2)), chunksize=1):
        with multiprocessing.Pool(processes=processes) as P:
            return [result for result in P.imap(func, Iter)]
         
    def large_iter_pmap(func, Iter, processes=CPUs, status_bar=True, nice=True, wait_interval=1):
        if nice:
            os.nice(10)
        try:
            with multiprocessing.Pool(processes=processes) as P:
                size = max(1, int(round(len(Iter)/CHUNKS)))
                rs = P.map_async(func, Iter, chunksize=size)
                while not rs.ready():
                    sleep(wait_interval)
                return rs.get()
 
        except PicklingError:
            warn("Lambda functions cannot be Pickled for Parallelization. Using single Process.", RuntimeWarning)
            return list(map(func, Iter))
 
class _mid_fastq_iter(object):
    def __init__(self, filename, seq_id, start, stop):
        self.filename = filename
        self.seq_id = seq_id
        self.start = start
        self.stop = stop
 
    def __iter__(self): 
        self.f = self.filename.open('rb')
        self.f.seek(self.start)
        # Find the beginning of the next FASTQ header
        lines = []
        for i, line in zip(range(5), self.f):
            lines.append(line)
            if line.startswith(self.seq_id):
                break
        else:
            raise RuntimeError("Took more than 4 lines to find header in middle of FASTQ file (start pos: {:}, stop pos: {:}, file length: {:}):\n".format(
                                    self.start, self.stop, self.filename.stat().st_size)+'\n'.join(lines))
        self.f.seek(self.f.tell() - len(line))
        return self
         
    def __next__(self):
        if self.f.tell() < self.stop:
            header = self.f.readline()
            dna = self.f.readline()
            self.f.readline()
            qc = self.f.readline()
            return header, dna, qc
        else:
            self.f.close()
            raise StopIteration
         
    def __exit__(self):
        self.f.close()
 
compressions = dict(gzip='gunzip',gz='gunzip', bz2='bunzip2', lzma='unxz', xz='unxz') 