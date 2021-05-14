"""
Project 01 Template File

CS/BIOS 112 - Project 01

   In this project, I will write a simple gene finder with a function that calls other functions, and then use it to find open reading frames in a 51,000 bp segment of human chromosome 9. 
   

@author:    <Adithya Jose>
Due Date:   <10/12/2020>
"""


# The following function should NOT be modified!!
def read_one_seq_fasta(fasta_file):
    """Read a FASTA file that contains one sequence."""
    seq = ''
    with open(fasta_file, 'r') as f:
        f.readline()  # by pass description in first line of data file
        for line in f.readlines():
            seq = seq + line[:-1]
    return seq
# The above function should NOT be modified!!!


def gc_content(seq):
    ''' Takes one argument, seq. It counts the number of G and C characters in the sequence and adds these counts and divides the GC count total by the length of the sequence and returns this fraction.  '''
    G_C = seq.count('G') + seq.count('C')
    divide = G_C / len(seq)
    return float(divide)


# Tests for gc_content function. Should print True in all cases.
print('\ngc_content Tests')
print(gc_content('ATGTGAA') == 0.2857142857142857)
print(gc_content('ATGAGATAAG') == 0.3)


def get_orf(seq):
    ''' Takes one argument, seq. The function assumes that the given DNA sequence begins with a start codon ATG. Finds the first in-frame stop codon, and returns the sequence from the start to that stop codon. The sequence that is returned includes the start codon but not the stop codon. If there is no in-frame stop codon, get_orf assumes that the reading frame extends through the end of the sequence and simply returns the entire sequence. '''
    es = ''
    if seq[:3] != 'ATG':
        return es
    else:
        for i in range (0, len (seq), 3):
            if(seq[i:i+3]) == 'TAG' or seq[i:i+3] == 'TGA' or seq[i:i+3] == 'TAA':
                return seq[0:i]
        return seq


# Tests for get_orf function. Should print True in all cases.
print('\nget_orf Tests')
print(get_orf('ATGTGAA') == 'ATG')
print(get_orf('ATGAGATAAG') == 'ATGAGA')
print(get_orf('ATGAGATAGG') == 'ATGAGA')
print(get_orf('ATGAGATGAGGGTAA') == 'ATGAGA')
print(get_orf('ATGAAATT') == 'ATGAAATT')


def one_frame(seq):
    ''' Takes one argument, seq. The argument (seq) is a string containing bases for a DNA sequence. Returns a list of ORFs contained in the DNA sequence. Finds the first in-frame start codon: 'ATG'. Calls get_orf( ) with the DNA sequence beginning at that start codon. Adds the ORF returned by get_orf( ) to a list. Looks for the next start codon beginning where the ORF ends. '''
    i = -3
    of_list = []
    while i <= len(seq):
        i = i + 3
        if seq[i: i + 3] == 'ATG':
            of_list.append(get_orf(seq[i:]))
            i = i + len(get_orf(seq[i:]))
    return of_list


# Tests for one_frame function. Should print True in all cases.
print('\none_frame')
print(one_frame('ATGTGAA') == ['ATG'])
print(one_frame('ATGAGATAAG') == ['ATGAGA'])
print(one_frame('ATGAGATAGG') == ['ATGAGA'])
print(one_frame('ATGAGATGAGGGTAA') == ['ATGAGA'])
print(one_frame('ATGAAATT') == ['ATGAAATT'])
print(one_frame('ATGAGATGAACCATGGGGTAA') == ['ATGAGA', 'ATGGGG'])


def forward_frames(seq):
    ''' Takes one argument, seq. The argument (seq) is a string containing bases for a DNA sequence. Returns a list of ORFs contained in the DNA sequence at ALL FRAMES. '''
    ff_list = [] 
    s = 0
    while s < 3:
        ff_list.extend(one_frame(seq[s:]))
        s = s + 1
    return ff_list


# Tests for forward_frames function. Should print True in all cases.
print('\nforward_frames')
print(forward_frames('ATGAGATAAG') == ['ATGAGA'])
print(forward_frames('ATGAGATGAGGGTAA') == ['ATGAGA', 'ATGAGGGTAA'])
print(forward_frames('ATGAAATT') == ['ATGAAATT'])
print(forward_frames('ATGAGATGACACCATGGGGTAA') == ['ATGAGA', 'ATGGGG', 'ATGACACCATGGGGTAA'])


def gene_finder(fasta_file, min_len, min_gc):
    ''' Take 3 arguments: First Argument: name of a data file, Second Argument: minimum length of ORFs to be found, Third argument: minimum GC content of ORF to be found
    orfs = gene_finder('gene_finder_test.fasta', 6, 0.45). Build a “list of lists” of ORFs contained in the data file that: that exceed the length requirement and exceed the GC content requirement. '''
    seq = read_one_seq_fasta(fasta_file)
    orf_list = []
    candidates = forward_frames(seq)
    for orf in candidates:
        if len(orf) >= min_len and gc_content(orf) >= min_gc:
            orf_list.append([len(orf), gc_content(orf), orf])     
    return orf_list
    

# Tests for gene_finder function. Should print True.
print('\ngene_finder')
calculated_result = gene_finder('gene_finder_test.fasta', 6, 0.45)
desired_result = [[6, 0.6666666666666666, 'ATGCCC'], [9, 0.7777777777777778, 'ATGCCCCGG']]
print( calculated_result == desired_result)

orf_list = gene_finder('human_chr9_segment.fasta', 550, 0.45)
print (len(orf_list) == 5)

# viewing the results of the gene_finder( ) calculations
orf_list = gene_finder('gene_finder_test.fasta', 6, 0.45)
print (orf_list)

"""

Identify ORFs in the provided X73525.fasta file

< The gene is the protein coding gene for the salmonella strain. The majority of the coding sequences for the genes is for proteins. Using gene_finder() for the X73525.fasta file, I used the following arguments for fasta_file = X73525.fasta, min_len = 550, min_gc = 0.45. My end result output or data was [[909, 0.528052805280528, 'ATGTCATTGCGTGTGAGACAGATTGATCGTCGCGAATGGCTATTGGCGCAAACCGCGACAGAATGCCAGCGCCATGGCCGGGAAGCGACGCTGGAATATCCGACGCGACAGGGAATGTGGGTTCGGTTGAGCGATGCAGAAAAACGGTGGTCGGCCTGGATTAAACCTGGGGACTGGCTTGAGCATGTCTCTCCCGCTCTGGCTGGGGCGGCGGTTTCTGCTGGCGCTGAGCACCTGGTCGTTCCCTGGCTTGCTGCAACAGAGCGACCGTTTGAGTTGCCCGTGCCGCATTTGTCCTGTCGGCGTTTATGCGTAGAGAACCCCGTACCGGGAAGCGCGCTGCCGGAAGGGAAATTGTTGCACATTATGAGCGATCGGGGCGGCCTGTGGTTTGAGCATCTTCCTGAACTGCCTGCAGTCGGGGGCGGCAGGCCGAAAATGCTGCGTTGGCCGTTGCGCTTTGTAATCGGTAGCAGTGATACGCAGCGTTCGTTGCTGGGCCGAATCGGGATCGGAGATGTACTCCTGATTCGTACTTCCCGTGCGGAAGTTTATTGCTACGCGAAAAAGTTAGGTCATTTCAACCGTGTTGAAGGGGGAATTATTGTGGAAACGTTAGATATTCAACATATCGAAGAAGAAAATAATACAACTGAAACTGCAGAAACTCTGCCTGGCTTGAATCAATTGCCCGTCAAACTGGAATTTGTTTTGTATCGTAAGAACGTTACCCTCGCCGAACTCGAAGCCATGGGGCAGCAACAGCTATTATCACTGCCGACCAATGCTGAACTTAACGTTGAAATTATGGCGAATGGTGTTTTGCTGGGTAATGGCGAACTGGTACAGATGAATGACACCTTAGGCGTTGAGATCCATGAATGGCTGAGCGAGTCTGGTAATGGGGAA'], [831, 0.5342960288808665, 'ATGGGCATTTTTGCCTCCGCAGGATGCGGTAAGACCATGCTGATGCATATGCTGATCGAGCAAACGGAGGCGGATGTCTTTGTTATCGGTCTTATCGGTGAACGAGGCCGTGAGGTCACTGAATTCGTGGATATGTTGCGCGCTTCGCATAAGAAAGAAAAATGCGTGCTGGTTTTTGCCACTTCCGATTTCCCCTCGGTCGATCGCTGCAATGCGGCGCAACTGGCGACAACCGTAGCGGAATATTTTCGCGACCAGGGAAAACGGGTCGTGCTTTTTATCGATTCCATGACCCGTTATGCGCGTGCTTTGCGAGACGTGGCACTGGCGTCGGGAGAGCGTCCGGCTCGTCGAGGTTATCCCGCCTCCGTATTCGATAATTTGCCCCGCTTGCTGGAACGCCCAGGGGCGACCAGCGAGGGAAGCATTACTGCCTTTTATACGGTACTGCTGGAAAGCGAGGAAGAGGCGGACCCGATGGCGGATGAAATTCGCTCTATCCTTGACGGTCACCTGTATCTGAGCAGAAAGCTGGCCGGGCAGGGACATTACCCGGCAATCGATGTACTGAAAAGCGTAAGCCGCGTTTTTGGACAAGTCACGACGCCGACACATGCTGAACAGGCATCTGCCGTGCGTAAATTAATGACGCGTTTGGAAGAGCTCCAGCTTTTCATTGACTTGGGAGAATATCGTCCTGGCGAAAATATCGATAACGATCGGGCGATGCAGATGCGGGATAGCCTGAAAGCCTGGTTATGCCAGCCGGTAGCGCAGTATTCATCCTTTGATGACACGTTGAGCGGTATGAATGCATTCGCTGACCAGAAT'], [1008, 0.4851190476190476, 'ATGGGCGATGTGTCAGCTGTCAGTTCATCCGGGAACATTTTACTGCCGCAGCAGGATGAGGTTGGCGGTTTATCAGAAGCATTAAAAAAAGCGGTGGAAAAACATAAGACAGAATATTCCGGTGATAAAAAAGATCGCGACTATGGCGATGCTTTCGTAATGCATAAAGAAACGGCTTTACCGTTATTACTGGCGGCATGGCGACATGGCGCGCCAGCGAAATCAGAACATCACAATGGCAACGTTTCTGGTCTGCATCATAACGGAAAAAGCGAACTCAGGATTGCTGAAAAACTGTTGAAAGTCACTGCTGAAAAATCTGTCGGTTTGATCTCTGCGGAGGCCAAAGTAGATAAATCCGCAGCGTTGCTATCGTCTAAAAATAGGCCGTTAGAAAGCGTAAGCGGTAAAAAATTATCTGCTGATTTAAAAGCTGTGGAATCCGTTAGTGAAGTAACCGATAACGCCACGGGAATCTCTGACGATAATATCAAGGCATTGCCTGGGGATAATAAAGCCATCGCGGGCGAAGGCGTTCGTAAAGAGGGCGCGCCGCTGGCGCGGGATGTCGCACCTGCCCGAATGGCCGCAGCCAATACCGGTAAGCCTGAAGATAAAGATCATAAAAAGGTTAAAGATGTTTCTCAGCTTCCGCTGCAACCAACCACTATCGCCGATCTTAGCCAATTAACCGGCGGCGATGAAAAAATGCCTTTAGCGGCGCAATCAAAGCCGATGATGACTATTTTTCCCACTGCCGATGGCGTGAAAGGAGAGGATAGCTCGCTGACTTACCGTTTTCAGCGCTGGGGAAATGACTATTCCGTCAATATTCAGGCGCGGCAAGCAGGGGAGTTTTCGTTAATACCGTCAAATACGCAGGTTGAACATCGTTTGCATGATCAATGGCAAAACGGTAATCCCCAGCGCTGGCACCTGACGCGAGACGATCAACAAAATCCGCAGCAGCAACAGCACAGACAGCAATCTGGCGAGGAGGATGACGCC'], [789, 0.4955640050697085, 'ATGTTTTACGCGTTGTACTTTGAAATTCATCACCTGGTTGCGTCTGCGGCGCTAGGGTTTGCTCGCGTGGCGCCGATTTTTTTCTTCCTGCCGTTTTTGAATAGCGGGGTATTAAGCGGTGCGCCGAGAAACGCCATTATCATCCTGGTGGCATTGGGAGTATGGCCGCATGCATTGAACGAGGCGCCGCCGTTTTTATCGGTGGCGATGATCCCGTTAGTTCTGCAAGAAGCGGCGGTAGGCGTCATGCTGGGCTGTCTGCTGTCATGGCCTTTTTGGGTTATGCATGCGCTGGGTTGTATTATCGATAACCAGCGAGGGGCAACGCTAAGTAGTAGTATCGATCCGGCAAACGGTATTGATACCTCGGAAATGGCTAATTTCCTGAATATGTTTGCCGCTGTCGTTTATTTACAAAACGGCGGTCTGGTCACGATGGTTGACGTGTTAAATAAAAGCTATCAGCTATGCGATCCGATGAACGAGTGCACGCCTTCATTACCGCCGCTATTAACGTTTATTAATCAGGTGGCTCAAAACGCCTTGGTTCTGGCCAGTCCGGTGGTATTAGTGCTGTTGCTGTCAGAAGTATTCCTGGGTTTATTGTCGCGCTTTGCTCCGCAAATGAACGCTTTTGCGATTTCACTGACGGTAAAAAGCGGTATTGCCGTTTTAATTATGCTGCTTTATTTCTCTCCGGTACTACCGGACAATGTACTGCGACTCTCTTTCCAGGCCACAGGGTTAAGCAGTTGGTTTTACGAGCGAGGGGCGACGCATGTCCTCGAA']].  >

"""
