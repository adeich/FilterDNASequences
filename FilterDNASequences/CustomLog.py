import time
import ConstantsAndStructures

class Log:
	
	class StatisticalAnalysisRecord:

		def __init__(self):
			self.iTotalSequencesSeen = 0
			self.iTotalSuccessfulSequences = 0

		def RecordSequenceSeen(self, bSuccess):
			self.iTotalSequencesSeen += 1
			if bSuccess:
				self.iTotalSuccessfulSequences += 1

		def GetTotalAndSuccesses(self):
			return self.iTotalSequencesSeen, self.iTotalSuccessfulSequences


	def __init__(self, oLogFile):
		self.oLogFile = oLogFile
		self.AnalysisRecord = self.StatisticalAnalysisRecord()

	def Log(self, sMessage):
		self.oLogFile.writeline(sMessage)

	def IngestReportAndLog(self, tReport):
		if type(tReport) == str:
			print tReport
		self.oLogFile.write(str(tReport))
		self.oLogFile.write('\n')
		self.AnalysisRecord.RecordSequenceSeen(tReport.bPasses_filters)
		
	def AddAnalysisReport(self):
		iTotal, iSuccesses = self.AnalysisRecord.GetTotalAndSuccesses()
		self.oLogFile.write('There {} successful sequences out of {} total seen.\n'.format(
		 	iSuccesses, iTotal))
