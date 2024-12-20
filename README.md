# Primer-Design
#### Description:
    Just to give quick context of what PCR is and how it might be usefull, PCR tests are used to amplify and certain segments of DNA rapidly to create many copies.

    It has three steps which are:
     -Denaturation (unwinding the DNA using heat),
     -Annealing (raising the tempature of the DNA just enough to facilitate bonding of the primers to the main DNA strand, but not enough to destroy the DNA), and
     -Elongation (making many copies of resulting DNA)
    A quick example everyone should be familiar with would be the covid-19 PCR tests which detected specific DNA sequences found in the corona virus (they were known to be the most accurate tests)

    I will breifly explain what each custom function accomplishes

    main
        Gathers user input, passes that input through other functions, manages return values and determines which errors are present
        - the return values of my custom functions are either 0 or 1 (or 2 in some cases). a return value of 0 generally means the given argument to the function was successful and that argument has a favorable effect on primer design process
    validate
        takes 1 argument, uses regex to validate a user inputed DNA sequence, uses sys.exit if invalid
        - uses regular expressions (re.search) to match user input
        - currently accepts either a DNA strand or RNA strand (in futre may do something with RNA)
    is_length
        makes sure that a DNA strand is 18 - 30 nucleotides long
        - just using conditionals
    temp_diff
        determines if the two melting tempatures are withing 5 degrees celcius of each other
        - uses abs function
    gc_content
        calculates the amount of guanine and cytosine nucleotides in given strand as a percentage, should be within 40 to 60 percent
        - creates list correspong to each nucleotide, uses for loop to append each nucleotide to corresponding list, get length of G and C lists, converts it to a percentage
    runs_repeats
        loops over the whole string and counts how many time adenine is followed by thymine or guanine is followed by cytosine
        - uses a for loop to check current nucleotide and the previous one at the same time
    primer_dimer
        checks both primers for any high complimentary base pair presence
        -Probably most complex function in the program. I use for loop to compare two DNA strands at once (forward and reverse primer)

    What lessons I've learned: I've learned from making this program that I now know why software engineers are considered engineers, beacuse they have to build something that is used by a large amount of people, good design is important, and you have to think practically

    A technical lesson I've learned from making this project would be making for loops. You almost have to picture it in your mind before you implement it.
    Also using Psuedocode really helps.