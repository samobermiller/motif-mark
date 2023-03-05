#!/usr/bin/env python
import argparse
def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--file", help="fasta file")
    parser.add_argument("-m", "--motif", help="motifs file")
    return parser.parse_args()
##OOP
class motif:
    def __init__(self, sequence,color,legend_text,legend_location):
        '''This is how a motif is made'''
        ## Data ##
        self.sequence = sequence
        self.color=color
        self.legend_text=legend_text
        self.legend_location=legend_location
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
locations=[[35,675],[35,700],[35,725],[35,750],[35,775]]
with open(args.motif, "r") as motif_list:
    for line in motif_list:
        line=line.strip('\n')
        item_fix=motif(line,color[i],line,locations[i])
        i+=1
        motif_classes.append(item_fix)
#possible_motifs now stores regex expressions for each motif given in args.motif


##Open files and start counters
import re
match=re.findall(r'[ \w-]+?(?=\.)', args.file)
output_filename=match[0]

#open canvas
import cairo
surface = cairo.SVGSurface(f"{output_filename}.svg", 1000, 800)
context = cairo.Context(surface)
context.set_source_rgb(1,1,1)
context.paint()

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
header_location=50
with open(f"./{output_filename}_one_line.fa", "r") as new_fasta:
    for line in new_fasta:
        if line.startswith(">"):
            header=line
            print(header)
            context.set_source_rgb (0, 0, 0)
            context.set_font_size(13)
            context.move_to(15,header_location)
            context.show_text(header)
            context.stroke()
            header_location+=60
        else:
            line=line.strip('\n')
            line_value+=60
            sequence_length=len(line)
            new_line=re.split(f'[a-z][A-Z]|[A-Z][a-z]',line)
            #print(new_line)
            start=10
            end=start
            for segment in new_line:
                #print(segment)
                segment_length=len(segment)
                end=start+segment_length
                if segment.isupper():
                    #print(f"exon length {segment_length}")
                    #print(f"start position {start}")
                    #print(f"end position {end}")
                    context.set_line_width(20)
                    context.move_to(start,line_value)
                    context.line_to(end,line_value)
                    context.set_source_rgb (0, 0, 0)
                    context.stroke()
                if segment.islower():
                    #print(f"intron length {segment_length}")
                    #print(f"start position {start}")
                    #print(f"end position {end}")
                    context.set_line_width(5)
                    context.move_to(start,line_value)
                    context.line_to(end,line_value)
                    context.set_source_rgb (0, 0, 0)
                    context.stroke()
                start+=segment_length
            motif_match=()
            location=()
            for single in motif_classes:
                #print(single.sequence)
                possible_motifs=[]
                possible_motifs+=single.fix_ambiguity(single.sequence)
                #print(single.sequence)
                for entry in possible_motifs:
                    #color assignment issue here
                    motif_match = ()
                    motif_match=re.findall(f'{entry}', line)
                    #print(motif_match)
                    for motif_seq in motif_match:
                        #print(motif_seq, single.color)
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
                        location=single.legend_location
                        location1=location[0]
                        location2=location[1]
                        #print(single.color)
                        context.set_line_width(25)
                        color1=single.color[0]
                        color2=single.color[1]
                        color3=single.color[2]
                        #print(single.color)
                       # print(single.sequence)
                        context.move_to(motif_start,line_value)
                        context.line_to(motif_end,line_value)
                        context.set_source_rgb(color1, color2, color3)
                        context.stroke()
                        context.set_font_size(13)
                        context.move_to(location1,location2)
                        context.show_text(single.sequence)
                        context.stroke()
                        #print(start)
                        #whats the length of the segment?
                        #is it upper or lower case?
context.set_source_rgb (0, 0, 0)
context.set_line_width(4)
context.rectangle(30, 645, 100, 125)
#context.rectangle(x coordinate of left side, y coordinate of top side, width, height)
context.stroke()
context.set_font_size(25)
context.move_to(15, 635)
context.show_text("Figure Legend")
context.stroke()
surface.write_to_png(f"{output_filename}.png")
surface.finish()
# for single in motif_classes:
#     print(single.color)
#print(motif_dict)