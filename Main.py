from collections import namedtuple
import SequenceFileZipper as SFZ


def Main(sFastaFileName, sQualiFileName, sOutputSequenceFileName, sLogFileName, sTagFileName):

	# set filenames
	sFastaFileName = ''
	sQualiFileName = ''
	sOutputSequenceFileName = 'test.fna'
	sLogFileName = 'test.log'
	sTagFileName = 'tag.txt' # tag file is deliberately omitted from github repo.


	with open(sOutputSequenceFileName, 'w') as oOuputFile, 
		open(sFastaFileName, 'r') as oFastaFile,
		open(sQualiFileName, 'r') as oQualiFile,
		open(sLogFileName, 'w') as oLogFile,
		open(sTagFileName, 'r') as oTagFile:

		# load sequence files into sequence generator.
		oSeqIterator = SFZ.SeqNamedTupleGenerator(oFastaFile, oQualiFile)
		
		# create filtering object
		oFilterBank = FilteringPipeline(oTagFile)

		# create the log object
		oLog = CustomLog(oLogFile)

		# iterate through sequences, checking each if it passes the filters.
		for sequence, id_string in seqIterator:
			tReport = oFilterBank.ProcessAndReturnReport(sequence, id_string)
			oLog.IngestReportAndLog(tReport)
			if tReport.bPasses_filters:
				oOuputFile.writeline(tReport.output_id)
				oOuputFile.writeline(tReport.output_seq)


# Checks individual sequences against several pattern matchings.
# Returns a namedtuple report on whether the sequence is fit for inclusion.
class FilteringPipeline:
	def __init__(self, oTagFile):
		self.oTagFile = oTagFile 

	def ProcessAndReturnReport(self, sSequence, id_string):
		pass

class CustomLog:
	def __init__(self, oLogFile):
		self.oLogFile = oLogFile

	def IngestReportAndLog(tReport):
		pass
