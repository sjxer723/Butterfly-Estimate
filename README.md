# Butterfly-Estimate

## Setup Enviroment (on Ubuntu)

You can download the dependencies using `apt`
```bash
$ sudo apt install libomp-dev
```
Then build the source codes by
```bash
$ cd src/
$ make
```

## Prepare Dataset
Our datasets are from the two websites [KONECT](https://konect.cc/) and [SNAP](https://snap.stanford.edu/snap/install.html). To download a specific dataset and decompress the graph file, please execute the following command.
```bash
$ cd dataset/
# Take http://konect.cc/networks/lastfm_band/ as example
$ python3 download_dataset.py -i http://konect.cc/networks/lastfm_band/
```

## Run Butterfly Counting Algorithms

* (Step 1): Convert the text file to binary format
Since it takes a long time to read a graph from a text file, we save it to a binary file by the following command,
    ```bash
    $ cd src/
    $ ./btfc -p ../dataset/[Graph Name]/ -m txt2bin 
    ```
    If the input graph is unipartite graph and you would like to convert it to a bipartite graph by partioning the vertices by partity, execute the following command with extra option `--uni-to-bip`,
    ```
    $ ./btfc -p ../dataset/[Graph Name]/ -m txt2bin --uni-to-bip
    ```
* (Step 2): Run butterfly counting algorithm. The meaning of common options are as follow,
    
    * `-p`: the path to the graph
    * `-m`: counting method, including `exact`, `espar`, `tls`, `wps`
    * `-q`: whether to count the number of queries
    * `-d`: the density of the graph and range from `0` to `10`, where `10` represents all edges are kept
    * `--s1`: the value of parameter $s_1$. If it is set to `3`, it means $s_1$ is set as $3\cdot \sqrt{m}$
    * `--s2`: the value of parameter $s_2$. If it is set to `3`, it means $s_2$ is set as $3\cdot \sqrt{m}$
    * `--uni-to-bip`: converting the unipartite graph to bipartite graph
    * `--rbase`, `--rexp`, `--rround`: run the sampling method for $r_{base} \cdot r_{exp}^{r_{round}}$ rounds and get the error rate after rounds $r_{base}\cdot r_{exp}^i$ for $i=0,\ldots, r_{round}$.
    
    Take the graph `lastfm_band` as an example

    **Exact Couting**:
    ```
    $ ./btfc -p ../dataset/lastfm_band -m exact -q
    ```
    **ESpar**:
    ```
    $ ./btfc -p ../dataset/lastfm_band -m espar -q
    ```
    **Our method (TLS, Two Leveled Sampling)**
    
    Set $s_1=\sqrt{m}$ and $s_2 = 5\sqrt{m}$ and $r=200$:
    ```
    $ ./btfc -p ../dataset/lastfm_band -m tls -q --s1 1 --s2 5 --rbase 200
    ```

    **WPS (Weighted Pair Sampling)**
    Set $r=20000$:
    ```
    $ ./btfc -p ../dataset/lastfm_band -m wps -q --rbase 20000
    ```
    If the graph is converted from a unipartite graph, run the following command instead,
    Set $r=20000$:
    ```
    $ ./btfc -p ../dataset/lastfm_band -m wps -q --rbase 20000 --uni-to-bip
    ```

## Varying Graph Density
To see the performance on graphs with different density, run the following command to randomly keep $20\%$ edges,
```bash
$ ./btfc -p ../dataset/[Graph Name]/ -m txt2bin -d 2
```
Then when running counting algorithms, also add the extra option `-d 2`, as follows,
```bash
$ ./btfc -p ../dataset/lastfm_band -m espar -q -d 2
```
    