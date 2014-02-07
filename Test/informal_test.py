# Add other directories to the PYTHONPATH. 
import sys, os
testdir = os.path.dirname(__file__)
srcdir = '../FilterDNASequences'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))

import SequenceFileZipper as SFZ
import Main
import FilterBank as FB
import CustomLog


b = Main.Main(
 	sOutputSequenceFileName='test_output_dump/test_created_generated_output_seq.txt',
	sFastaFileName='test_data/300LinesTest.fna',
	sQualiFileName='test_data/300LinesTest.qual',
	sLogFileName='test_output_dump/log_from_informal_test.txt',
	sTagFileName='test_data/test_tags.csv')


with open('test_data/test_tags.csv', 'r') as oTagFile, open('test_output_dump/log2.txt', 'w') as oLogFile:

	oLog = CustomLog.Log(oLogFile)

	a = FB.FilterBank(oTagFile, oLog)

	print str(a.oLoadedTags.ReturnTagList())
	print str(a.oIsATissueTag.dAllPossibleSubtags)

