import re
import sys

# This program does not check it, but the minimum length between primers should be 150 base pairs
# or, at most 10 kilo base pairs
def main():
    user_input_fp = input("Forward Primer: ").upper()
    user_input_rp = input("Reverse Primer: ").upper()
    valid_fp = validate(user_input_fp)
    valid_rp = validate(user_input_rp)
    check_islengthfp = is_length(valid_fp)
    check_islengthrp = is_length(valid_rp)
    fp_temp = int(input('What is the melting tempature of the Forward primer in Celsius? '))
    rp_temp = int(input('What is the melting tempature of the Reverse primer in Celsius? '))
    check_diff = temp_diff(fp_temp, rp_temp)
    check_gccontentfp = gc_content(valid_fp)
    check_gccontentrp = gc_content(valid_rp)
    check_runsrepeatsfp = runs_repeats(valid_fp)
    check_runsrepeatsrp = runs_repeats(valid_rp)
    check_primerdimer = primer_dimer(valid_fp, valid_rp)
    error_counter = 0
    if check_islengthfp == 1:
        print("Forward primer is longer than 30 nucleotides")
        error_counter += 1
    if check_islengthfp == 2:
        print('Forward primer is shorter than 18 nucleotides')
        error_counter += 1
    if check_islengthrp == 1:
        print('Reverse primer is longer than 30 nucleotides')
        error_counter += 1
    if check_islengthrp == 2:
        print('Reverse primer is shorter than 18 nucleotides')
        error_counter += 1
    if check_diff == 1:
        print('The difference of the melting tempatures of the primers is greater than 5 degrees')
        error_counter += 1
    if check_gccontentfp == 1:
        print('The GC content in the Forward primer is out of the min/max bounds')
        error_counter += 1
    if check_gccontentrp == 1:
        print('The GC content in the Reverse primer is out of the min/max bounds')
        error_counter += 1
    if check_runsrepeatsfp == 1:
        print('The Forward primer has runs and repeats that are greater than the current allowed limit')
        error_counter += 1
    if check_runsrepeatsrp == 1:
        print('The Reverse primer has runs and repeats that are greater than the current allowed limit')
        error_counter += 1
    if check_primerdimer == 1:
        print('The current forward and reverse primers have high complimentary sequences together, a primer dimer will likely form')
        error_counter += 1
    if error_counter == 0:
        print('Both forword and reverse primers are suitable sequences')
    
    # print error counter
    if error_counter == 0:
        print('Congratz your Primers have 0 errors')
    print(f"Error counter: {error_counter}")

# This function makes sure that the user only inputs a valid sequence
def validate(s):    
    pattern = r'^[ATGC]+$|^[AUGC]+$'
    match = re.search(pattern, s)
    if match:
        return s
    else:
        sys.exit("Only enter a valid DNA sequence") 


# This function only checks to see if given primer is within minimum, maximum length requirements
def is_length(seq):
    if len(seq) > 30:
        return 1
    elif len(seq) < 18:
        return 2
    else:
        return 0

# This function checks the melting tempature differential to make sure it less then or equal to 5
def temp_diff(x, y):

    absolute_value = abs(x - y)
    if absolute_value <= 5:
        return 0
    else:
        return 1

# This function checks the G,C content of the given primer
def gc_content(s1):
    # Get length of primer
    s1_length = int(len(s1))
    # Create lists corresponding to each nuclotide
    flsa = []
    flst = []
    flsg = []
    flsc = []
    # Iterating over each nucleotide in sequence, assigning nucleotide to appropriate list
    for nucleotide in s1:
        if nucleotide == 'A':
            flsa.append(nucleotide)
        elif nucleotide == 'T':
            flst.append(nucleotide)
        elif nucleotide == 'G':
            flsg.append(nucleotide)
        elif nucleotide == 'C':
            flsc.append(nucleotide)
    # Checking GC content is within bounds
    len_g = len(flsg)
    len_c = len(flsc)
    len_gc = len_g + len_c
    if 40 <= ((len_gc/s1_length) * 100) <= 60:
        return 0
    else:
        return 1
    
# This function checks for runs and repeats
def runs_repeats(sequence):
    ls_seqat_counter = 0
    ls_seqgc_counter = 0
    for nucleotide in range(len(sequence) - 1):
        if sequence[nucleotide] == 'A' and sequence[nucleotide+1] == 'T':
            ls_seqat_counter += 1
        if sequence[nucleotide] == 'G' and sequence[nucleotide+1] == 'C':
            ls_seqgc_counter += 1
    if ls_seqat_counter >= 3 or ls_seqgc_counter >= 3:
        return 1
    else:
        return 0
    


# This function checks both primers for any high complimentary base pair presence
def primer_dimer(x, y):
    match_occurences_AT = []
    consecutive_AT = 0
    max_consecutive_AT = 0
    for index, (nuc_a, nuc_t) in enumerate(zip(x, y)):
        if nuc_a == 'A' and nuc_t == 'T':
            match_occurences_AT.append(index)
            if len(match_occurences_AT) > 1 and match_occurences_AT[-1] == match_occurences_AT[-2] + 1:
                consecutive_AT += 1
            else:
                consecutive_AT = 1
            max_consecutive_AT = max(max_consecutive_AT, consecutive_AT)


    match_occurences_GC = []
    consecutive_GC = 0
    max_consecutive_GC = 0
    for index, (nuc_a, nuc_t) in enumerate(zip(x, y)):
        if nuc_a == 'G' and nuc_t == 'C':
            match_occurences_GC.append(index)
            if len(match_occurences_GC) > 1 and match_occurences_GC[-1] == match_occurences_GC[-2] + 1:
                consecutive_GC += 1
            else:
                consecutive_GC = 1
            max_consecutive_GC = max(max_consecutive_GC, consecutive_GC)
    #setting max consecutive complimentary base pairs to 3, total complimentary base paris to 4    
    if max_consecutive_AT < 4 and len(match_occurences_AT) < 5 and max_consecutive_GC < 4 and len(match_occurences_GC) < 5:
        return 0
    else:
        return 1



    





if __name__ == "__main__":
    main()