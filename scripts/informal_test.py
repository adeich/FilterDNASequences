import SequenceFileZipper as SFZ
import Main
import FilterBank as FB
import CustomLog



b = Main.Main(
 	sOutputSequenceFileName='../testing_dir/generated_output_seq.txt',
	sFastaFileName='../test_data/300LinesTest.fna',
	sQualiFileName='../test_data/300LinesTest.qual',
	sLogFileName='../testing_dir/log.txt',
	sTagFileName='../test_data/test_tags.csv')


with open('../test_data/test_tags.csv', 'r') as oTagFile, open('../testing_dir/test2.txt', 'w') as oLogFile:

	oLog = CustomLog.Log(oLogFile)

	a = FB.FilterBank(oTagFile, oLog)
	a.TestAllClasses()

	print str(a.oLoadedTags.ReturnTagList())
	print str(a.ClassesDict['IsATissueTag'].dAllPossibleSubtags)

