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
    ## Method ##
    def fix_ambiguity(self,sequence):
        "take motif and generate regex expression that can be used to search sequence while accounting for any ambiguity from motif or sequence"
        sequence=sequence.upper()
        fix = {'Y': "[CTUctu]", 'T':"[tTUu]",'A':"[Aa]",'G':"[Gg]",'C':"[Cc]",'U':"[UuTt]",'W':"[AaTtUu]",'S':"[CcGg]",'M':"[AaCc]",'K':"[GgTtUu]",'R':"[AaGg]",'B':"[CcGgTtUu]",'D':"[AaGgTtUu]",'H':"[AaCcTtUu]",'V':"[AaCcGg]",'N':"[AaCcGgTtUu]"}
        sequence=list(sequence)
        fixed_sequence=str()
        for letter in sequence:
            fixed_sequence+=fix[letter]
        return[fixed_sequence]
        #creates regular expression based on given motif that you can use to search input sequence
#TEST for class motif
# practice=motif('catag')
# amb=practice.fix_ambiguity(practice.sequence)
# print(amb[0])

##Generate motif regex statements
args=get_args()
motifs=[]
with open(args.motif, "r") as motif_list:
    for line in motif_list:
        motifs.append(line.strip())
possible_motifs=[]
for item in motifs:
    item_fix=motif(item)
    possible_motifs+=item_fix.fix_ambiguity(item)
#print(possible_motifs[1])
#possible_motifs now stores regex expressions for each motif given in args.motif
import re
match=re.findall(r'[ \w-]+?(?=\.)', args.file)
output_filename=match[0]
output=open(f'{output_filename}.fa', "w")
count_lines=0
#Open canvas
import cairo
##Parse motif and fasta file by header
headers=0
#number of headers
start=10
#start point, running sum
total_bases=0
#total number of nucleotides
lower_bases=0
#intron bases
upper_bases=0
#exon bases
motif_base=0
#bases in motifs
count_nomotif=0
count_motif=0
#number of motifs found
surface = cairo.SVGSurface(f"{output_filename}_plot.svg", 1000, 100)
context = cairo.Context(surface)
with open(args.file,"r") as fasta:
    #PER LINE
    for line in fasta:
        if line.startswith(">"):
            headers+=1
            output.write(line)
        else:
            line=line.strip('\n')
            motif_match=()
            for entry in possible_motifs:
                motif_match=re.findall(f'{entry}', line)
                #list of motif matches in the current line based on regex created with class motif
            #if the number of motifs in the line is not 0 then:
            if len(motif_match)!=0:
                #print(motif_match)
                #PER MOTIF
                for motif_seq in motif_match:
                    #for each motif match in the line
                    count_motif+=1
                    pre_motif=re.findall(f'.*(?={motif_seq})', line)
                    pre_motif=pre_motif[0]
                    #pre_motif is the sequence before the current motif
                    #PER BASE
                    for base in pre_motif:
                    #differentiate exons and introns before motif match
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
                        motif_base+=1
                        context.move_to(start,40)
                        context.set_line_width(30)
                        context.move_to(start,40)
                        context.set_source_rgba(4, 0, 4, 0.5)
                        context.line_to(start+1,40)
                        context.stroke()                   
            if len(motif_match)==0:
                count_nomotif+=1
                context.set_source_rgba (0, 0, 0, 1)
                #PER BASE
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
surface.write_to_png(f"{output_filename}_plot.png")
surface.finish()
print(f"# of Headers= {headers}")
print(f"Total Lines in File= {count_lines}")
print(f"# of Motifs= {count_motif}")
print(f"# of Motif Bases= {motif_base}")
print(f"# of Exon Bases= {upper_bases}")
print(f"# of Introns Bases= {lower_bases}")
print(f"# Total Bases= {total_bases}")


##Current issue: not all motifs are depicted in figure getting drawn. keep making counts to find where issue is, total bases checks out