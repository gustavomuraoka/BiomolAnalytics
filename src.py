# This project aims to provide the user plenty of useful info about biological subjects

# This class seeks to describe as close as needed a Poly nucleotide
class Polynucleotyde:

    def __init__(self, sequence):
        self.sequence = sequence


# This class describe the DNA segment
class DNASegment(Polynucleotyde):

    def __init__(self, sequence):
        super().__init__(sequence)
        # This dictionary helps with DNA -> RNA translation
        self.translation_dict = {"A": "U", "G": "C", "C": "G", "T": "A"}

    # Simple way to translate DNA to RNA
    def DNA_to_RNA(self):

        try:
            translated_seq = ""
            for nucleotide in self.sequence:
                translated_seq += self.translation_dict[nucleotide]
            return translated_seq

        # Provide info when the sequence has a letter that shouldn't be part of a DNA sequence
        except KeyError:
            for nucleotide in range(len(self.sequence)):
                if self.sequence[nucleotide] != "A" and self.sequence[nucleotide] != "T" and self.sequence[nucleotide] != "C" and self.sequence[nucleotide] != "G":
                    print(f'Your sequence seem to have an error, please try again [Position: {nucleotide} ({self.sequence[nucleotide]})]')


# This class describe the RNA segment
class RNASegment(Polynucleotyde):

    def __init__(self, sequence):
        super().__init__(sequence)
        # This dictionary helps with DNA -> RNA translation
        self.translation_dict = {"U": "A", "G": "C", "C": "G", "A": "T"}

    # Simple way to translate RNA to DNA
    def RNA_to_DNA(self):

        try:
            translated_seq = ""
            for nucleotide in self.sequence:
                translated_seq += self.translation_dict[nucleotide]
            return translated_seq

        # Provide info when the sequence has a letter that shouldn't be part of an RNA sequence
        except KeyError:
            for nucleotide in range(len(self.sequence)):
                if self.sequence[nucleotide] != "A" and self.sequence[nucleotide] != "U" and self.sequence[nucleotide] != "C" and self.sequence[nucleotide] != "G":
                    print(f'Your sequence seem to have an error, please try again [Position: {nucleotide} ({self.sequence[nucleotide]})]')


# This dictionary is useful when translating mRNA to Protein, having one-letter and three-letter codes, both being
# common ways to deal with aminoacids referencing.
# Stop codons have 'STOP' code twice to avoid bugs due the 'codon_dict[codon][X] format'
codon_dict = {'UUU': ['Phe', 'F'], 'UUC': ['Phe', 'F'], 'UUA': ['Leu', 'L'], 'UUG': ['Leu', 'L'],
              'CUU': ['Leu', 'L'], 'CUC': ['Leu', 'L'], 'CUA': ['Leu', 'L'], 'CUG': ['Leu', 'L'],
              'AUU': ['Ile', 'I'], 'AUC': ['Ile', 'I'], 'AUA': ['Ile', 'I'], 'AUG': ['Met', 'M'],
              'GUU': ['Val', 'V'], 'GUC': ['Val', 'V'], 'GUA': ['Val', 'V'], 'GUG': ['Val', 'V'],
              'UCU': ['Ser', 'S'], 'UCC': ['Ser', 'S'], 'UCA': ['Ser', 'S'], 'UCG': ['Ser', 'S'],
              'CCU': ['Pro', 'P'], 'CCC': ['Pro', 'P'], 'CCA': ['Pro', 'P'], 'CCG': ['Pro', 'P'],
              'ACU': ['Thr', 'T'], 'ACC': ['Thr', 'T'], 'ACA': ['Thr', 'T'], 'ACG': ['Thr', 'T'],
              'GCU': ['Ala', 'A'], 'GCC': ['Ala', 'A'], 'GCA': ['Ala', 'A'], 'GCG': ['Ala', 'A'],
              'UAU': ['Tyr', 'Y'], 'UAC': ['Tyr', 'Y'], 'UAA': ['STOP', 'STOP'], 'UAG': ['STOP', 'STOP'],
              'CAU': ['His', 'H'], 'CAC': ['His', 'H'], 'CAA': ['Gln', 'Q'], 'CAG': ['Gln', 'Q'],
              'AAU': ['Asn', 'N'], 'AAC': ['Asn', 'N'], 'AAA': ['Lys', 'K'], 'AAG': ['Lys', 'K'],
              'GAU': ['Asp', 'D'], 'GAC': ['Asp', 'D'], 'GAA': ['Glu', 'E'], 'GAG': ['Glu', 'E'],
              'UGU': ['Cys', 'C'], 'UGC': ['Cys', 'C'], 'UGA': ['STOP', 'STOP'], 'UGG': ['Trp', 'W'],
              'CGU': ['Arg', 'R'], 'CGC': ['Arg', 'R'], 'CGA': ['Arg', 'R'], 'CGG': ['Arg', 'R'],
              'AGU': ['Ser', 'S'], 'AGC': ['Ser', 'S'], 'AGA': ['Arg', 'R'], 'AGG': ['Arg', 'R'],
              'GGU': ['Gly', 'G'], 'GGC': ['Gly', 'G'], 'GGA': ['Gly', 'G'], 'GGG': ['Gly', 'G']}


# This class describes proteins
class Protein(DNASegment):
    def __init__(self, sequence, is_DNA):
        super().__init__(sequence)
        if is_DNA is True:
            self.sequence = super().DNA_to_RNA()
        self.translation_dict = {"U": "A", "G": "C", "C": "G", "A": "U"}

    # This function translates the RNA to the protein primary structure, where you are able to set up
    # how you expect to deal with the frames, giving better perspectives on the sequence
    def nucleotide_to_protein(self):
        chosen_frame = int(input("Which frame you want to use? (0 = All): "))
        if chosen_frame > 0:
            self.singleframe_setup(chosen_frame)

    def singleframe_setup(self, chosen_frame):
        if 0 < chosen_frame < 4:
            resultant_seq = self.sequence[chosen_frame-1:]
            print(self.frame_toProtein_translation(resultant_seq))
        else:
            resultantseq = self.reversestrand_sequencer()[chosen_frame - 4:]
            print(self.frame_toProtein_translation(resultantseq))

    @staticmethod  # Remove this if needed anytime soon
    def frame_toProtein_translation(self, sequence):
        protein_seq, codon = "", ""
        for nucleotide in sequence:
            codon += nucleotide
            if len(codon) == 3:
                protein_seq += codon_dict[codon][0] + " - "
                codon = ""
        return protein_seq

    def reversestrand_sequencer(self):
        reverse_strand = ""
        for nucleotide in self.sequence:
            reverse_strand += self.translation_dict[nucleotide]
        return reverse_strand


# This variable is responsible to control the program exec, making it stop whenever the user wants to
loop_control = True

while loop_control:

    # Display to the user the possibilities available
    print(f'1. DNA -> mRNA Translation\n'
          f'2. RNA -> DNA Translation\n'
          f'3. RNA/DNA -> Protein Translation')
    option_selected = int(input("Select: "))

    # !! Atualizar pra uma versão mais inteligente (switch/match) quando possivel Python 3.9 nn dá pra usar match!!
    if option_selected == 1:
        DNA_seq = DNASegment(input("Insert your sequence: ").upper())
        print(f'mRNA output: {DNA_seq.DNA_to_RNA()}')

    elif option_selected == 2:
        RNA_seq = RNASegment(input("Insert your sequence: ").upper())
        print(f'DNA output: {RNA_seq.RNA_to_DNA()}')

    elif option_selected == 3:
        nucleotide_seq = input("Insert the sequence to be translated into protein:").upper()
        if nucleotide_seq.find("T") > -1 and nucleotide_seq.find("U") > -1:
            print("Your sequence contains both Uracil and Thymine, try again!")
        elif nucleotide_seq.find("T") > -1:
            protein = Protein(nucleotide_seq, True)
            protein.nucleotide_to_protein()
        elif nucleotide_seq.find("U") > -1:
            protein = Protein(nucleotide_seq, False)
        else:
            print("Your sequence might have some error, check again")

    loop_control = False
