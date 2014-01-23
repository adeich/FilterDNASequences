import re
from collections import namedtuple
import Exceptions
import ConstantsAndStructures
import csv
import SequenceAnalysisClasses 


class FilterBank:

	class LoadedTissueTags:

		def __init__(self, oTissueTagFile):
			self.tTissueTag = namedtuple('tTissueTag', ['origin', 'label', 'tag'])
			self.lLoadedTissueTagTuples = self.LoadTissueTagsFromFile(oTissueTagFile)
		
		def LoadTissueTagsFromFile(self, oTagFile):
			lLoadedTagList = []
			oCSVReader = csv.reader(oTagFile, delimiter = ',')
			for lRow in oCSVReader:
				lLoadedTagList.append(self.tTissueTag(*lRow))
			return lLoadedTagList
	
		def ReloadTagFile(self, oTagFile):
			self.lLoadedTissueTagTuples = self.LoadTissueTagsFromFile(oTagFile)

		def ReturnTagList(self):
			return self.lLoadedTissueTagTuples

	def __init__(self, oTagFile, oLog):
		self.oLog = oLog
		self.oTagFile = oTagFile
		self.oLoadedTags = self.LoadedTissueTags(oTagFile)
		self.ClassesDict = {
			'ContainsForwardAndReversePrimers': SequenceAnalysisClasses.ContainsForwardAndReversePrimers(),
			'ContainsForwardAndReversePrimers_Complement': SequenceAnalysisClasses.ContainsForwardAndReversePrimers_Complement(),
			'IsATissueTag': SequenceAnalysisClasses.IsATissueTag(self.oLoadedTags.ReturnTagList()),
			'ContainsBothFlankingSequences': SequenceAnalysisClasses.ContainsBothFlankingSequences(),
			'ContainsBothFlankingSequences_Complement': SequenceAnalysisClasses.ContainsBothFlankingSequences_Complement()
		}
		

	def TestAllClasses(self):
		for oClass in self.ClassesDict:
			try:
				self.ClassesDict[oClass].TestSelf()
			except AttributeError as e:
				oLog.Log("Looks like a testing class is missing the TestSelf() method")
				raise	
				

	# Returns a ConstantsAndStructures.tSequenceReport namedtuple.
	def RunCompositeAnalysisOnSequence(self, sIDString, sCompleteSequence, sQualiSequence):

		# a simple class for keeping failure diagnoses externally non-reversable.
		class StatusFailureRecorder:
			def __init__(self):
				self.bSequenceSucceeds = True 
			def SetFailureStatus(self):
				self.bSequenceSucceeds = False
			def ReturnSuccessStatus(self):
				return self.bSequenceSucceeds

		oFailureRecorder = StatusFailureRecorder()
		bSequenceIsReversed = None
		lErrorsWithSequence = []

		# is sequence in forward direction?
		if self.ClassesDict['ContainsForwardAndReversePrimers'].Run(sCompleteSequence):
			bSequenceIsReversed = False

			# does sequence contain proper forward and reverse primers?
			if not self.ClassesDict['ContainsBothFlankingSequences'].Run(sCompleteSequence):
				oFailureRecorder.SetFailureStatus()
				lErrorsWithSequence.append('Flanking sequences aren\'t right.')

			# does sequence begin with a proper tissue tag?
			sPossibleTissueTag = self.ClassesDict[
				'ContainsForwardAndReversePrimers'].ReturnSequencePrependingForwardPrimer(
					sCompleteSequence)
			if not self.ClassesDict['IsATissueTag'].Run(sPossibleTissueTag, bMatchEntireTagOnly=False):
				oFailureRecorder.SetFailureStatus()
				lErrorsWithSequence.append('Tissue tag does not match tests.')

			# does the insert sequence pass its tests?
			tInsSeqBegEndPos = self.ClassesDict[
				'ContainsBothFlankingSequences'].ReturnInsSeqBegEndPos(sCompleteSequence)


		# otherwise, is sequence in reverse direction?
		elif self.ClassesDict['ContainsForwardAndReversePrimers_Complement'].Run(sCompleteSequence):
			bSequenceIsReversed = True
			if self.ClassesDict['ContainsBothFlankingSequences_Complement'].Run(sCompleteSequence):
				pass
				# if begins with tissue tag
				# if insert sequence passes tests

		# otherwise, the sequence is bad.
		else:
			oFailureRecorder.SetFailureStatus()
			lErrorsWithSequence.append('Contains no forward-reverse primer pairs, normal or complement.')
			

		tOutputReport = ConstantsAndStructures.tSequenceReport(
			bPasses_filters = oFailureRecorder.ReturnSuccessStatus(),
			printready_output_string = 'hello',
			input_id = sIDString,
			input_seq = sCompleteSequence,
			output_id = sIDString,
			output_seq = sCompleteSequence,
			seq_is_reversed = bSequenceIsReversed,
			start_pos_forward_primer = 0,
			end_pos_forward_primer = 0,
			start_pos_ending_seq = 0,
			end_pos_ending_seq = 0
			)

		return tOutputReport
