# DynaHubs
DynaHubs reveals critical residues and their frequency of occurrence that contribute to allosteric communication throughout the MD simulation. DynaHubs computes the contact topology for every frame produced along the MD and identifies the residues with a betweenness value in the top 5%. This calculation determines the frequency of high betweenness residues over the whole molecular dynamics simulation for each frame. 

### Background
DynaHubs utilizes Residue Interaction Network (RIN) to identify hub rezidues. Using a bidirectional graph to characterize the contact topology, RIN can pinpoint allosteric residues/regions of a protein complex [reference 1] (https://journals.tubitak.gov.tr/biology/vol42/iss5/3/) . According to RIN, the protein structure is a network with nodes connected by edges. The local interaction magnitude a<sub>ij </sub>between adjacent residue pairs (i, j) determines the lengths of the edges.

![aij](https://github.com/user-attachments/assets/8b7efb6c-3bb9-4ffa-b053-156a3ac2938d)

where N<sub>ij </sub>is the total number of ith and jth residue nonhydrogen atom-atom interactions within a 4.5 Å cutoff distance. The number of nonhydrogen atoms in the ith and jth residues, respectively, is represented by N<sub>i </sub>and N<sub>j </sub>, which are applied to weight N<sub>ij </sub>with the goal to diminish the impact of amino acid size on the magnitude of the local interaction. Close nearby residue (node) pairs are thought to have high connections in this model and are able to share information, such as fluctuations. When calculating the edge length between nodes i and j, 1/a<sub>ij </sub>is used to remove the substantial bias toward the covalent interactions.

RIN depends on protein topology. Accordingly, centrality measurements play a very useful role in exposing the protein network's topological characteristics. In order to create allosteric communication pathways between far-off locations, betweenness centrality (C<sub>B) </sub>measurement is employed  to identify hub residues with a high potential for information transmission and reception through tertiary interactions. C<sub>B </sub>centrality is determined in this way:

![cb](https://github.com/user-attachments/assets/49dd7c69-2630-4882-954a-907ae0436102)

The shortest number of routes between nodes i and j is denoted by σ<sub>ij </sub>, while the shortest number of routes that traverse node l is denoted by σ<sub>ij</sub>(l). High betweenness hub residue regions (top 5% of C<sub>B</sub>) are considered as allosteric locations that may be assessed as potential drug target regions [reference 2] (https://pubs.acs.org/doi/10.1021/acs.jpcb.4c00925).

### DynaHubs Architecture

![image](https://github.com/user-attachments/assets/efc20187-3b86-4fe3-a435-1cc19ee3420b)


## Usage 
### System Dependencies
- python3 (3.8 or higher)
- ipkernel
### Python Dependencies
- glob
- itertools
- matplotlib.pyplot as plt
- mdanalysis
- networkx
- numpy 

## Installation

**We recommend you use [git](https://git-scm.com/downloads) clone instead of download with .zip while installing the package. 

### Clone the repository
```
git clone https://github.com/inantu/DynaHubs.git
```
```
cd DynaHubs
```
## Run from terminal
```
python dynahubs.py --pdb  --dcd  --stride 10 --top_betweenness 5

```
