import SequenceFileZipper as SFZ
import Main
import FilterBank as FB
import CustomLog



b = Main.Main(
 	sOutputSequenceFileName='testing_dir/generated_output_seq.txt',
	sFastaFileName='300LinesTest.fna',
	sQualiFileName='300LinesTest.qual',
	sLogFileName='testing_dir/log.txt',
	sTagFileName='Human6-Tags.txt')


with open('testing_dir/test1.txt', 'w') as oTagFile, open('testing_dir/test2.txt', 'w') as oLogFile:

	oLog = CustomLog.Log(oLogFile)

	a = FB.FilterBank(oTagFile, oLog)
	a.TestAllClasses()


