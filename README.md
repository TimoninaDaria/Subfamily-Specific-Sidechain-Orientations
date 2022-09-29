Subfamily-Specific-Sidechain-Orientations
=========================================

Usage:   python3 main.py aligned_pdbs=\</path/to/folder\>
aligned_fasta=\</path/to/file\> output=\</path/to/folder\> [options]

Example: python3 main.py aligned_pdbs=./input_pdbs aligned_fasta=./input.fasta
output=./output

 

Mandatory input parameters:

===========================

aligned_pdbs=\<string\>  \# Path to folder with aligned protein 3D-structures as
separate files in the PDB format (each file should represent one chain)

aligned_fasta=\<string\> \# Path to the corresponding sequence representation of
the alignment in the FASTA format

output=\<string\>        \# Path to folder to store results

 

Utilization of computing resources:

===================================

cpu_threads=\<int\> \# Number of parallel CPU threads to utilize (the default is
"all" physically available)

 

Cluster analysis methods:

=========================

method=hdbscan              \# Use HDBSCAN automatic method (default)

method=optics               \# Use OPTICS automatic method

method=dbscan eps=\<float\>   \# Use DBSCAN method for manual fine-tuning of the
results by specifying the ‘eps’ value (eps \> 0)

 

Cluster size:

=========================

min_samples=\<int\>           \# The 'min_samples' parameter of HDBSCAN, OPTICS,
and DBSCAN (the number of points in a neighborhood

for a point to be considered as the cluster core)

min_cluster_size=\<int\>      \# The HDBSCAN 'min_cluster_size' parameter to
regulate the minimal size of a cluster

 

Selection of common core positions:

===================================

max_content_of_gaps=\<int\>      \# Define the allowed gap content in alignment
column, in % (at most 5%, by default)

max_content_of_mismatch=\<int\>  \# Define the allowed 3D-mismatch content in
alignment column, in % (at most 5%, by default)

mismatch_threshold=\<float\>     \# Define the cut-off value to discriminate
spatially aligned from misaligned residues, in angstroms (selected
automatically)

 

Printing result:

===================================

number_of_result_resids=\<int\>     \# Number of residues with the greatest
score to print in RESULT-files (10, by default)

ref=\<string\>                      \# Set the reference protein by its name in
the FASTA alignment (default is the first protein)

 

PyMol parameters:

=================

compile_pymol_pse=false \# Disable compilation of PyMol sessions with
3D-annotation of results
