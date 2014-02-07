import time
import ConstantsAndStructures

class Log:
	
	class StatisticalAnalysisRecord:

		def __init__(self):
			self.iTotalSequencesSeen = 0
			self.iTotalSuccessfulSequences = 0
			self.ErrorsSeen = {}

		def RecordSequenceSeen(self, bSuccess, lErrors):
			self.iTotalSequencesSeen += 1
			if bSuccess:
				self.iTotalSuccessfulSequences += 1
			for sError in lErrors:
				if sError not in self.ErrorsSeen:
					self.ErrorsSeen[sError] = 0
				else:
					self.ErrorsSeen[sError] += 1

		def GetTotalAndSuccesses(self):
			return self.iTotalSequencesSeen, self.iTotalSuccessfulSequences

		def GetErrorsDict(self):
			return self.ErrorsSeen


	def __init__(self, oLogFile):
		self.oLogFile = oLogFile
		self.AnalysisRecord = self.StatisticalAnalysisRecord()

	def Log(self, sMessage):
		self.oLogFile.writeline(sMessage)

	def IngestReportAndLog(self, tReport):
		if not tReport.bPasses_filters:
			self.oLogFile.write('Bad sequence:\n')
			self.oLogFile.write('\tID: {}\n'.format(tReport.input_id));
			self.oLogFile.write('\tSeq: {}\n'.format(tReport.input_seq));
			self.oLogFile.write('\tQuali: {}\n'.format(tReport.input_quali_seq));
			self.oLogFile.write('\tErrors: {}\n\n'.format(', '.join(tReport.errors_within_sequence))); 
		self.AnalysisRecord.RecordSequenceSeen(tReport.bPasses_filters, 
			tReport.errors_within_sequence)
		
	def AddAnalysisReport(self):
		iTotal, iSuccesses = self.AnalysisRecord.GetTotalAndSuccesses()
		self.oLogFile.write('There {} successful sequences out of {} total seen.\n'.format(
		 	iSuccesses, iTotal))
		self.oLogFile.write('Total sequence errors:\n')
		for sError in self.AnalysisRecord.GetErrorsDict():
			self.oLogFile.write('\t"{}":\n\t\t occurred {} times.\n'.format(sError,
				self.AnalysisRecord.GetErrorsDict()[sError]))
