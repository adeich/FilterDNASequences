from collections import namedtuple


# Contains all information about a given sequence.
tSequenceReport = namedtuple('tSequenceReport', 
		['bPasses_filters', 
		'printready_output_string',
		'input_id',
		'input_seq',
		'output_id',
		'output_seq',
		'seq_is_reversed',
		'start_pos_forward_primer',
		'end_pos_forward_primer',
		'start_pos_ending_seq',
		'end_pos_ending_seq'])


sFORWARD_PRIMER = 'CGCAATTCCTTTAGTTGTTC'
sFORWARD_PRIMER_REVERSED = 'GAACAACTAAAGGAATTGCG'
sREVERSE_PRIMER = 'AACCTCATACAGAAAATTCA'
sREVERSE_PRIMER_REVERSED = 'TGAATTTTCTGTATGAGGTTTTGC'
sBEGINNING_FLANKING_SEQUENCE = ''



def ParseTissueTagsFile():
	pass

dTissueTagsDictionary = ParseTissueTagsFile()
