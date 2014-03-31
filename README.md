FilterDNASequences
==================

This package reads through a fasta + quali file pair and produces a new, trimmed down fasta file. The new file contains only those DNA sequences from the original which have passed all tests. 

The program is initiated by calling `Main()` in `FilterDNASequences/Main.py`.


#####General Usage Instructions
- Check that `FilterDNASequences/ConstantsAndStructures.py` contains the correct primer and flanking sequences for you. If incorrect, alter this file accordingly.
- You must supply a *tissue tag* csv file for `Main()`. Here is an example of such a file (the header line *is* required):
```
ORIGIN,LABEL,TAG
testpatient1,sample1,ATT
testpatient1,sample2,ACTA
testpatient1,sample3,ATA
testpatient1,sample4,ATC
testpatient2,sample5,AATC
testpatient2,sample6,ACA
```
- Your supplied fasta and quali file pair should be in the regular format.
