from randomDNA import random_DNA
from DNA_Codons import *
import sys
import re


class Bioinfo():

    def __init__(self, sequence):
        ''' Programm accepts .fasta and .txt files as files. Both will be flatted to contain only characters 
        and no further information about the sequence'''
        
        self.sequence = sequence

        if not isinstance(self.sequence, str):
            print("Only Strings are accepted")
            sys.exit(1)

        if self.sequence.endswith(".fasta"):
            with open(sequence, "r") as seq_file:
                sequence_plain = "".join(seq_file.readlines()[1::]).replace("\n","")
                self.sequence = sequence_plain

        elif self.sequence.endswith(".txt"):
            with open(sequence, "r") as seq_file:
                sequence_plain = "".join(seq_file.readlines()).replace("\n","")
                self.sequence = sequence_plain

        else:
            print("No valid format or path: programm accepts .fasta or .txt (typo?)")
            sys.exit(1)

    def sequence_identifer(self):
        ''' Checks if string is 1) empty, 2) includes lowercases 3) is DNA 4) or is Protein '''

        regex_dna = (r"[^ATGCN]")
        regex_aa = (r"[^GPAVLIMCFYWHKRQNEDST_]")


        if not self.sequence:
            print("empty file")
            sys.exit(1)

        elif self.sequence.isupper() == False:
            print("sequence with lower charaters")
            sys.exit(1)

        else:
            if re.search(regex_dna, self.sequence) == None:
                return "DNA detected"
            elif re.search(regex_aa, self.sequence) == None:
                return "Protein detected"
            else:
                print("No biological sequence")
                return False 

    def count_characters(self, *character):
        ''' count elements and adds the occurency into  dict, with the searched element as key '''
        count_dict = {}

        for element in character:
            counter = self.sequence.count(element)
            count_dict[element] = counter

        return count_dict

    def change_characters(self, template, translate_dict):
        ''' takes a sequence and changes wanted elements, elements must be given in a dict as as a key,
        while the replacements are the value'''

        changed_list = []

        if not isinstance(translate_dict, dict):
            sys.exit(1)

        if not (isinstance(template, str) or isinstance(template, list) or isinstance(template, tuple)):
            sys.exit(1) 

        for element in template:
            if element in translate_dict:
                changed_list.append(translate_dict[element])
            else:
                changed_list.append(element)

        return "".join(changed_list)
            
    def splice_string_by_position(self, slice_lenght, steps=1):
        ''' Gives slices of a sequence of wanted size, usefull for codon generation from a DNA template'''
        
        fragments = []

        if slice_lenght > len(self.sequence):
            print("slice bigger than sequence")
            sys.exit(1)

        if steps > len(self.sequence):
            print("steps bigger than sequence")
            sys.exit(1)

        for element in range(0, len(self.sequence), steps):
            fragments.append(self.sequence[element:(element+slice_lenght)])

        return fragments

    def splice_string_by_charaters(self, start_position, stop_position):
            '''Gives slices of a sequence with a specific start and end , 
            usefull for identifing genes in a aminoacoid sequence''' 

            strings = []

            if start_position not in self.sequence:
                print("start character not in sequence")
                sys.exit(1)
            
            if stop_position not in self.sequence:
                print("stop character not in sequence")
                sys.exit(1)

            for position, element in enumerate(self.sequence):
                if element == start_position:
                    start_pos = position
                    for i in range(start_pos, len(self.sequence)):
                        strings.append(self.sequence[i])
                        if self.sequence[i] == stop_position:
                            break

            strings = "".join(strings).split(stop_position)
            if strings [-1] == "":
                strings.remove("")
            return strings


class Bioinfo_DNA(Bioinfo):

    def __init__(self, sequence):
        super().__init__(sequence)
        if self.sequence_identifer() != "DNA detected":
            print("Only DNA allowed" )
            sys.exit(1)
        else:
            pass

    def return_DNA(self):
        return self.sequence

    def return_RNA(self):
        rna_dict = {"T":"U"}
        return self.change_characters(self.sequence, rna_dict)

    def reverse_sequence(self):
        return self.sequence[::-1]

    def complementary_sequence(self):
        complementary_dict = {"A":"T", "T":"A", "G":"C", "C":"G"}
        return self.change_characters(self.sequence, complementary_dict)[::-1]

    def nucleotide_content(self):
            nucleotide = self.count_characters("A", "T", "G", "C")
            return ("Adenin: " + str(nucleotide["A"]) + 
                    " Thymin: " + str(nucleotide["T"]) + 
                    " Guanin: " + str(nucleotide["G"]) +
                    " Cytosin: " + str(nucleotide["C"]))

    def GC_content(self):
            nucleotide = self.count_characters("A", "T", "G", "C")
            GC = (nucleotide["G"] + nucleotide["C"])/(nucleotide["A"] + nucleotide["T"]+ nucleotide["G"] + nucleotide["C"])*100
            GC_round = round(GC,2)
            return ("GC content: " + str(GC_round))

    def count_motif(self, motif):
        return self.count_characters(motif)

    def codons(self):
        return self.splice_string_by_position(3,3)

    def translation(self, codon_list):
        return self.change_characters(codon_list, DNA_Codons)


class Bioinfo_Protein(Bioinfo):

    def __init__(self, sequence):
        super().__init__(sequence)

        if self.sequence_identifer() == "DNA detected":
            codon_list = self.splice_string_by_position(3,3)
            protein = self.change_characters(codon_list, DNA_Codons)
            self.sequence = protein
        
        if self.sequence_identifer() == "Protein detected":
            pass
            
    def return_AA(self):
        return self.sequence

    def find_gene(self):
        return self.splice_string_by_charaters("M","_")






### RandomString Generator ####
#random_String = random_DNA(200)

#sequence = Bioinfo(test_fasta_empty)
#print(sequence.sequence_identifer())
#print(sequence.count_characters("Z"))

### Methods for DNA ###
# sequence = Bioinfo_DNA()
# print(sequence.return_DNA())
# print(sequence.return_RNA())
# print(sequence.reverse_sequence())
# print(sequence.complementary_sequence())
# print(sequence.nucleotide_content())
# print(sequence.GC_content())
# print(sequence.count_motif("TATA"))
# print(sequence.codons())
#print(sequence.translation(sequence.codons()))

### Methods for Proteins ###
#sequence = Bioinfo_Protein(test_fasta_dna)
#print(sequence.return_AA())
#print(sequence.find_gene())




