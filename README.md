# broadinstitute-long-reads
## Dana Rocha 

First created: 2/28/21 
Last updated: 3/17/21

compute_kmer.py addresses the 4 criteria outlined in the ACS Take-home question: 
  1. Allows a user to choose k-mer sizes of 31, 41, 51, 61, 71, 81, 91, 101
  2. Computes and plots the k-mer coverages for the chosen k-mer size
  3. Determines where the local minimum between the first and second peak is 
  4. Displays the number of k-mers that have coverage greater than the local minimum value 
  
I used Python 3.7 and the PyCharm IDE to write this program. 



compute_kmer.py should be executed on the command line with the following commands:

> python3 compute_kmer.py -i reads.fixed.txt

or 

> python3 compute_kmer.py --infile reads.fixed.txt




This program assumes a regular text file as the input - not a gzip-compressed text file. 

The user will be prompted to enter a K-mer size. 
The program will output a figure that plots K-mer Coverage vs. Count. 
The graph will indicate where the threshold value is with a dashed red line. 

The threshold value and number of k-mers with coverage greater than the threshold is shown on the graph and command line. 
