from collections import namedtuple
import SequenceFileZipper
import FilterBank
import CustomLog


def Main(sFastaFileName, sQualiFileName, sOutputSequenceFileName, sLogFileName, sTagFileName):

	with open(sOutputSequenceFileName, 'w') as oOuputFile, open(sFastaFileName, 'r') as oFastaFile, open(sQualiFileName, 'r') as oQualiFile, open(sLogFileName, 'w') as oLogFile, open(sTagFileName, 'r') as oTagFile:

		# create the log object
		oLog = CustomLog.Log(oLogFile)

		# load sequence files into sequence generator.
		oSeqIterator = SequenceFileZipper.SeqNamedTupleGenerator(oFastaFile, oQualiFile)
		
		# create filtering object
		oFilterBank = FilterBank.FilterBank(oTagFile, oLog)

		# iterate through sequences, checking each if it passes the filters.
		for sIDString, sSequence, sQualiSequence in oSeqIterator:
			tReport = oFilterBank.RunCompositeAnalysisOnSequence(sIDString, sSequence, sQualiSequence)
			oLog.IngestReportAndLog(tReport)
			if tReport.bPasses_filters:
				oOuputFile.write(''.join([tReport.output_id, '\n']))
				oOuputFile.write(''.join([tReport.output_seq, '\n']))

