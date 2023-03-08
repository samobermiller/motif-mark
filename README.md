# motif-mark README

python script using object-oriented code to visualize motifs on a sequence using Pycairo

## Environment:

Runs on base python with pycairo

## Input (arg parse):

-m "file of possible motifs, max 5 with at most 10 bases each, one motif per line"
example:
ygcy
GCAUG
catag
YYYYYYYYYY

-f "fasta file, max 10 sequences with at most 1000 bases each, where exons are indicated by capitol letters"
example:
>INSR chr19:7150261-7150808 (reverse complement)
atgtccacatgtagtcacgtttgacatcccagggccacctcagcaggccgtctctgggga
gaattttctctgatttcttccccttcccttgctggacccctgcacctgctggggaagatg

## Output:

{filename}_one_line.fa
"one line version of input fasta file where exons are indicated by capital letters"

{filename}.png and {filename}.svg
"image visualizing exons, introns and motifs in a sequence"

## How The Code Works:
This code accounts for ambigous nucleotides in either the motif or the sequence by using Object Oriented Programming (OOP) to create a regex expression for each of the given motifs. The key of ambiguous nucleotides can be found here: https://en.wikipedia.org/wiki/Nucleic_acid_notation

The input fasta file is parsed through and rewritten into a one line fasta file named {filename}_one_line.fa. 

The code parses each sequence in the new fasta file, distinguishing introns from exons based on whether the base is capitalized and identifying motifs based on the list of regex expressions created earlier. Introns and exons are drawn in black on the Pycairo canvas and distinguished by line thickness (exons thicker than introns). The motif matches are drawn slightly thicker than the exons and colored by the original motif.

A key is created at the bottom of the image. The image is then exported as both an svg and png file. 

Overlapping motifs will be colored by last match. All elements (motifs, introns, exons) are to scale by base length. 

## Example Image Output:
![image](https://user-images.githubusercontent.com/105182636/223596075-351d2290-79f8-4385-85e6-bd528c2593aa.png)

