#!/usr/bin/env python
import argparse
def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--file", help="fasta file")
    parser.add_argument("-m", "--motif", help="motifs file")
    return parser.parse_args()
##OOP
class motif:
    def __init__(self, sequence,color):
        '''This is how a motif is made'''
        ## Data ##
        self.sequence = sequence
        self.color=color
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

##Test class motif generation
# practice=motif('catag')
# amb=practice.fix_ambiguity(practice.sequence)
# print(amb[0])

args=get_args()
possible_motifs=[]
i=0
motif_classes=[]
color=[[0,0.2,0.8],[0.6,0,0.4],[1,0,0],[0.5,0.2,0.3],[0.3,0.7,0.2]]
with open(args.motif, "r") as motif_list:
    for line in motif_list:
        line=line.strip('\n')
        item_fix=motif(line,color[i])
        i+=1
        possible_motifs+=item_fix.fix_ambiguity(line)
        motif_classes.append(item_fix)
#possible_motifs now stores regex expressions for each motif given in args.motif


##Open files and start counters
import re
match=re.findall(r'[ \w-]+?(?=\.)', args.file)
output_filename=match[0]

#open canvas
import cairo
surface = cairo.SVGSurface(f"{output_filename}_plot.svg", 1000, 500)
context = cairo.Context(surface)

##Make Fasta file input into single line sequence
sequence=""
header=""
with open(f"./{output_filename}_one_line.fa", "w") as output:
    with open(args.file,"r") as fasta:
        for line in fasta:
            line=line.strip('\n')
            if ">" in line:
                if header != "":
                    output.write(header+"\n")
                    output.write(sequence+"\n")
                    sequence=""
                header=line
            else:
                sequence+=line
        output.write(header+"\n")
        output.write(sequence+"\n")
motif_dict={}
#motif_seq: start,end
##Parse new fasta file
line_value=10
with open(f"./{output_filename}_one_line.fa", "r") as new_fasta:
    for line in new_fasta:
        if line.startswith(">"):
            continue
        else:
            line=line.strip('\n')
            line_value+=40
            sequence_length=len(line)
            new_line=re.split(f'[a-z][A-Z]|[A-Z][a-z]',line)
            #print(new_line)
            start=10
            end=start
            for segment in new_line:
                #print(segment)
                segment_length=len(segment)
                end=start+segment_length
                context.set_source_rgb (0, 0, 0)
                if segment.isupper():
                    #print(f"exon length {segment_length}")
                    #print(f"start position {start}")
                    #print(f"end position {end}")
                    context.set_line_width(30)
                    context.move_to(start,line_value)
                    context.line_to(end,line_value)
                    context.stroke()
                if segment.islower():
                    #print(f"intron length {segment_length}")
                    #print(f"start position {start}")
                    #print(f"end position {end}")
                    context.set_line_width(5)
                    context.move_to(start,line_value)
                    context.line_to(end,line_value)
                    context.stroke()
                start+=segment_length
            motif_match=()
            for single in motif_classes:
                possible_motifs+=single.fix_ambiguity(single.sequence)
                #print(single.color)
                for entry in possible_motifs:
                    motif_match=re.findall(f'{entry}', line)
                    for motif_seq in motif_match:
                        #for each motif match in the segment
                        pre_motif=re.findall(f'.*(?={motif_seq})', line)
                        #motif_color=
                        #print(motif_seq)
                        #print(pre_motif[0])
                        motif_length=len(motif_seq)
                        motif_start=len(pre_motif[0])
                        motif_end=motif_start+motif_length
                        motif_dict[motif_seq]=motif_start,motif_end
                        #print(motif_start)
                        #print(motif_length)
                        #print(f"motif length {motif_length}")
                        #print(f"start position {motif_start}")
                        #print(f"end position {motif_end}")
                        context.set_line_width(30)
                        color1=single.color[0]
                        color2=single.color[1]
                        color3=single.color[2]
                        print(single.color)
                        
                        context.set_source_rgb(color1, color2, color3)
                        context.move_to(motif_start,line_value)
                        context.line_to(motif_end,line_value)
                        context.stroke()
                        #print(start)
                        #whats the length of the segment?
                        #is it upper or lower case?
surface.write_to_png(f"{output_filename}_plot.png")
surface.finish()
# for single in motif_classes:
#     print(single.color)