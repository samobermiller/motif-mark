#!/usr/bin/env python
import argparse
def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--file", help="fasta file")
    parser.add_argument("-m", "--motif", help="motifs file")
    return parser.parse_args()
class motif:
    def __init__(self, sequence):
        '''This is how a motif is made'''
        ## Data ##
        self.sequence = sequence
    ## Methods ##
    def fix_ambiguity(self,sequence):
        "take motif and generate regex expression that can be used to search sequence while accounting for any ambiguity from motif or sequence"
        sequence=sequence.upper()
        fix = {'Y': "[CTUctu]", 'T':"[tTUu]",'A':"[Aa]",'G':"[Gg]",'C':"[Cc]",'U':"[UuTt]",'W':"[AaTtUu]",'S':"[CcGg]",'M':"[AaCc]",'K':"[GgTtUu]",'R':"[AaGg]",'B':"[CcGgTtUu]",'D':"[AaGgTtUu]",'H':"[AaCcTtUu]",'V':"[AaCcGg]",'N':"[AaCcGgTtUu]"}
        sequence=list(sequence)
        fixed_sequence=str()
        for letter in sequence:
            fixed_sequence+=fix[letter]
        #print(fixed_sequence)
        return[fixed_sequence]
        #creates regular expression based on given motif that you can use to search input sequence
#TEST for class motif
# practice=motif('catag')
# amb=practice.fix_ambiguity(practice.sequence)
# print(amb[0])

class gene:
    def __init__(self, sequence):
        '''This is how an intron is made'''
        ## Data ##
        self.sequence = sequence
    ## Methods ##
    #def draw(self):
        #draw thin line for intron, draw thick line for exon

##Generate motif regex statements
args=get_args()
motifs=[]
with open(args.motif, "r") as motif_list:
    for line in motif_list:
        motifs.append(line.strip())
#print(motifs)
possible_motifs=[]
for item in motifs:
    item_fix=motif(item)
    #print(item_fix.fix_ambiguity(item))
    possible_motifs+=item_fix.fix_ambiguity(item)
#print(possible_motifs[1])
#possible_motifs now stores regex expressions for each motif given in args.motif
import re
match=re.findall(r'[ \w-]+?(?=\.)', args.file)
output_filename=match[0]
#print(output_filename)
output=open(f'{output_filename}.fa', "w")
count_lines=0
#Open canvas
import cairo

##Parse motif and fasta file by header
start=10
total_bases=0
lower_bases=0
upper_bases=0
count_nomotif=0
count_motif=0
surface = cairo.SVGSurface(f"{output_filename}_plot.svg", 1000, 100)
context = cairo.Context(surface)
with open(args.file,"r") as fasta:
    for line in fasta:
        if line.startswith(">"):
            #print(line)
            output.write(line)
        else:
            line=line.strip('\n')
            #print(line)
            motif_match=()
            for entry in possible_motifs:
                motif_match=re.findall(f'{entry}', line)
                #if motif_match not empty then:
            if len(motif_match)!=0:
                count_motif+=1
                #print(motif_match)
                #FIGURE OUT WAY TO GO THROUGH LINE AND DRAW INTRONS OR EXONS AND DIFFERENTIATE MOTIFS. MAYBE THIS NEEDS TO BE ANOTHER OOP???
                for motif_seq in motif_match:
                    #for each motif match in the line
                    #print(seq)
                    #print(line)
                    pre_motif=re.findall(f'.*(?={motif_seq})', line)
                    #sequence before the motif
                    pre_motif=pre_motif[0]
                    #print(sequence)
                    for base in pre_motif:
                        total_bases+=1
                        start+=1
                        context.move_to(start,40)
                        if base.islower():
                            #print(base)
                            lower_bases+=1
                            context.set_line_width(5)
                            context.move_to(start,40)
                            context.line_to(start+1,40)
                            context.stroke()
                        if base.isupper():
                            upper_bases+=1
                            context.set_line_width(30)
                            context.move_to(start,40)
                            context.line_to(start+1,40)
                            context.stroke()
                    for base in motif_seq:
                        total_bases+=1
                        start+=1
                        context.move_to(start,40)
                        if base.islower():
                            #print(base)
                            lower_bases+=1
                            context.set_line_width(5)
                            context.move_to(start,40)
                            context.set_source_rgba(4, 0, 4, 0.5)
                            context.line_to(start+1,40)
                            context.stroke()
                        if base.isupper():
                            upper_bases+=1
                            context.set_line_width(30)
                            context.move_to(start,40)
                            context.line_to(start+1,40)
                            context.stroke()                        
            if len(motif_match)==0:
                count_nomotif+=1
                context.set_source_rgba (0, 0, 0, 1)
                for base in line:
                    total_bases+=1
                    start+=1
                    #context.move_to(leftright,updown)
                    context.move_to(start,40)
                    if base.islower():
                        #print(base)
                        lower_bases+=1
                        context.set_line_width(5)
                        context.move_to(start,40)
                        context.line_to(start+1,40)
                        context.stroke()
                    if base.isupper():
                        upper_bases+=1
                        context.set_line_width(30)
                        context.move_to(start,40)
                        context.line_to(start+1,40)
                        context.stroke()
        count_lines+=1
surface.write_to_png("ooca_motif_plot.png")
surface.finish()
# # print(count_nomotif)
# # print(count_motif)
# print(upper_bases)
# print(lower_bases)
# print(count_lines)
# print(total_bases)