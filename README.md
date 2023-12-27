# Overview
CoroGraph: Bridging Work Efficiency and Cache Efficiency for Graph Algorithm Execution

CoroGraph is a graph framework implemented based on Galois parallel library.
CoroGraph aims to optimize the cache locality for irregular access patterns while maintaining the graph's preferred execution model (work efficiency).

# Build CoroGraph

## Dependencies

CoroGraph builds, runs, and has been tested on GNU/Linux.

- A modern C++ compiler compliant with the C++-20 standard (gcc >= 11)
- CMake (>= 3.13)
- Boost library (>= 1.58.0, we recommend building/installing the full library)
- libfmt (>= 4.0)

## Compile and Run CoroGraph

```
SRC_DIR=`pwd` # Or top-level CoroGraph source dir
BUILD_DIR=<path-to-your-build-dir>

mkdir -p $BUILD_DIR
cmake -S $SRC_DIR -B $BUILD_DIR -DCMAKE_BUILD_TYPE=Release
cd $BUILD_DIR
make ..
```

## Graph Input Format
We use the adjacent graph format from the Problem Based Benchmark Suite (http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html).
Example Graph Format:
```
AdjacencyGraph
<n>
<m>
<o0>
<o1>
...
<o(n-1)>
<e0>
<e1>
...
<e(m-1)>
```


## Run Graph Algorithms

```
cd $BUILD_DIR
./app/sssp/crg-sssp -f /path/to/your/graph.adj -t <threads num> -delta <delta> 
./app/k-core/crg-kcore -f /path/to/your/graph.adj  -t <threads num>
./app/pr/crg-pr -f /path/to/your/graph.adj  -t <threads num>
./app/cc/crg-cc -f /path/to/your/graph.adj  -t <threads num>
```

The latest version of CoroGraph will be maintained [here](https://github.com/xiangyuzhi/corograph).



