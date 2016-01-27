# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Coleman Ellis

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###

#Define a dictionary for easy complement switching
complements = {'A':'T','T':'A','C':'G','G':'C'}
def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

        Added two doctests, so all 4 possible nucleotides are tested.
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    return complements[nucleotide]


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string of only the characters 'A','T','C', and 'G'
        returns: the reverse complementary DNA sequence represented as a string

        Added a couple doctests for short strings
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("C")
    'G'
    >>> get_reverse_complement("AG")
    'CT'
    """

    #First pass:
    '''#initialize an empty string
    reverse_complement = ''

    #for each character (nucleotide) in the reversed dna string
    for nucleotide in dna[::-1]:
        #add the complement to the reverse_complement string
        reverse_complement += get_complement(nucleotide)

    return reverse_complement'''

    #List comprehension:
    return ''.join([get_complement(nucleotide) for nucleotide in dna[::-1]])


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        Added a doctest for a DNA strand without a stop codon that isn't
        some multiple of three bases, so that it must be added on
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAG")
    'ATGAG'
    """
    #initialize an empty string to store the ORF
    ORF = ''

    #break dna into a list of length 3 strings (codons)
    #this is taken directly from a stackoverflow response
    for codon in [dna[i:i+3] for i in range(0, len(dna), 3)]:

        #if a codon isn't complete, add it on and break the loop
        if len(codon) < 3:
            ORF += codon
            break

        #if a codon is a stop codon, break the loop
        elif aa_table[codon] == '|':
            break

        #otherwise, add the codon to the ORF
        else:
            ORF += codon

    return ORF


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        Added a doctest that tests both nested start codons and a strand
        that isn't divisible by 3
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGATGATGATGATGTAGCC")
    ['ATGATGATGATGATG']
    """
    #initialize empty list to store ORFs
    ORFs = []

    #use a while loop so we can manipulate i
    i = 0
    while i < len(dna):

        #if there's not a whole codon left, break the loop
        if len(dna[i:])<3:
            break

        #if the codon at index i is a start codon
        if aa_table[dna[i:i+3]] == 'M':

            #find the rest of the ORF, passing in the whole dna
            #string from the start codon, and add the result to ORFs
            ORF = rest_of_ORF(dna[i:])
            ORFs.append(ORF)

            #advance i the length of the ORF, to avoid overlapping ORFs
            i += len(ORF)

        #advance i by one codon
        i += 3

    return ORFs



def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        This doctest should be sufficient, it has ORFs in all three reference
        frames. This shouldn't need to test for length (not divisible by 3),
        because find_all_ORFs_oneframe() does that.
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """

    #original
    '''
    #initialize an empty list for ORFs
    ORFs = []

    #find all the ORFs for all three reference frames
    ORFs += find_all_ORFs_oneframe(dna)
    ORFs += find_all_ORFs_oneframe(dna[1:])
    ORFs += find_all_ORFs_oneframe(dna[2:])
    return ORFs
    '''

    #list comprehensions
    ORFs = [ORF for i in range(0,3) for ORF in find_all_ORFs_oneframe(dna[i:])]
    return ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    #original
    '''
    ORFs = []

    ORFs += find_all_ORFs(dna)
    ORFs += find_all_ORFs(get_reverse_complement(dna))

    return ORFs
    '''

    #list comprehensions
    dnas = [dna,get_reverse_complement(dna)]
    ORFs = [ORF for temp_dna in dnas for ORF in find_all_ORFs(temp_dna)]
    return ORFs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
