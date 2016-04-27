# KADABRA: an ADaptive Algorithm for Betweenness via Random Approximation

## Synopsis

This repository contains the source code of KADABRA, an adaptive probabilistic algorithm for 
approximating the betweenness centrality of nodes in a network. 
The algorithm is presented in the paper *ARXIV id coming soon*.

## Installation 

The software requires the [OpenMP API](http://openmp.org/wp/). After cloning this repository,
build the software by issuing the `make` command inside the project folder: 
    $ make
    g++ -fopenmp -std=c++11 -Ofast -Wall -pedantic -g -Iinclude main.cpp src/* -lm -o kadabra

To test the software, see the [code examples below](#usage-and-code-example).

## License

This work is licensed under the Apache License 2.0. 

## Usage and Code Example 

The network has to be provided as a file containing two space-separated
integers `u v` per line, for each edge `(u,v)` of the network. The labels of
the vertices are assumed to be consecutive. 

The following examples are based on the networks provided in the `example_input` folder of
the project, taken from the [SNAP dataset](http://snap.stanford.edu/data/index.html).

Computing all centralities in the undirected graph facebook_combined, with
maximum error 0.0005 with probability 0.9 (one output every 60 seconds, as
default).

    $ ./kadabra 0.002 0.1 example_input/facebook_combined.txt 
    Undirected graph
    Number of nodes: 4039
    Number of edges: 176468
    
    Finished after 722480 iterations.
    Edges visited: 4550877633
    Average edges visited: 6298
    Total time: 5.41213
    Time bfs: 5.16474
    Time critical: 0.0596783
    Time compute finished: 0.104325
    (Printing thread: 0)
    Maximum confidence interval: 0.00234721


Computing all centralities in the directed graph p2p-Gnutella08, with maximum
error 0.0005 with probability 0.9 (one output every second).

    $ ./kadabra -v 1 -d 0.0005 0.1 example_input/p2p-Gnutella08.txt 
    Directed graph
    Number of nodes: 6301
    Number of edges: 20777
    
    Situation after 2217695 iterations.
    Edges visited: 176258926
    Average edges visited: 79
    Total time: 1.00005
    Time bfs: 0.78371
    Time critical: 0.0585964
    Time compute finished: 0.113125
    (Printing thread: 7)
    Maximum confidence interval: 0.000810047
    
    Finished after 3498902 iterations.
    Edges visited: 271524930
    Average edges visited: 77
    Total time: 1.53687
    Time bfs: 1.14729
    Time critical: 0.0852268
    Time compute finished: 0.180261
    (Printing thread: 0)
    Maximum confidence interval: 0.000499986


Computing the 3 most central vertices in the undirected graph facebook_combined, with maximum error 0.001 with probability 0.9 (one output every 5 seconds).

    $ ./kadabra -v 5 -k 3 0.001 0.1 example_input/facebook_combined.txt 
    Undirected graph
    Number of nodes: 4039
    Number of edges: 176468
    
    Situation after 649056 iterations.
    Edges visited: 4224362767
    Average edges visited: 6508
    Total time: 5.00007
    Time bfs: 4.82385
    Time critical: 0.0560443
    Time compute finished: 0.0922923
    (Printing thread: 2)
           1)      107 0.473286 0.481347 0.489408
           2)     1684 0.329856 0.337917 0.345978
         ? 3)     3437 0.230624 0.238348 0.241497
      Others) <0.236979
    
    Finished after 1010328 iterations.
    Edges visited: 6489038902
    Average edges visited: 6422
    Total time: 7.68489
    Time bfs: 7.21628
    Time critical: 0.0838349
    Time compute finished: 0.142478
    (Printing thread: 0)
           1)      107 0.476197 0.481375 0.486554
           2)     1684 0.332754 0.337933 0.343111
           3)     3437 0.233782 0.238778 0.240798


Computing the 10 most central vertices in the directed graph p2p-Gnutella08, with maximum error 0.0005 with probability 0.9 (one output every second).

    $ ./kadabra -d -v 1 -k 10 0.0005 0.1 example_input/p2p-Gnutella08.txt 
    Directed graph
    Number of nodes: 6301
    Number of edges: 20777
    
    Situation after 2131545 iterations.
    Edges visited: 169802374
    Average edges visited: 79
    Total time: 1.00005
    Time bfs: 0.784716
    Time critical: 0.0571349
    Time compute finished: 0.115175
    (Printing thread: 6)
         ? 1)     1317 0.017343 0.0188347 0.0201337
         ? 2)        3 0.0160498 0.0172678 0.0188218
         ? 3)      146 0.0138388 0.0145167 0.0152279
         ? 4)      175 0.0133574 0.0140494 0.0147774
         ? 5)      390 0.0132523 0.0139275 0.014637
         ? 6)      559 0.0111538 0.011828 0.012543
         ? 7)     1534 0.0110589 0.0117309 0.0124439
         ? 8)      250 0.0105703 0.0112501 0.0119736
         ? 9)      700 0.00967956 0.010362 0.0110926
        ? 10)      264 0.00955875 0.0102503 0.010992
      Others) <0.00915037
    
    Finished after 4060826 iterations.
    Edges visited: 313388533
    Average edges visited: 77
    Total time: 1.81202
    Time bfs: 1.34299
    Time critical: 0.09832
    Time compute finished: 0.218097
    (Printing thread: 0)
           1)     1317 0.0180754 0.0188747 0.0195466
           2)        3 0.0166251 0.0172758 0.0180752
           3)      146 0.0142079 0.0145685 0.0149382
           4)      390 0.0137179 0.0140784 0.0144485
           5)      175 0.0136059 0.0139725 0.014349
           6)      559 0.0115459 0.011906 0.0122773
           7)     1534 0.0113632 0.0117208 0.0120896
           8)      250 0.0108522 0.0112137 0.0115873
           9)      700 0.0100112 0.0103755 0.010753
          10)      264 0.00987839 0.0102474 0.0106303
