import re
from collections import namedtuple
import Exceptions
import ConstantsAndStructures
import csv
import SequenceAnalysisClasses 


class FilterBank:

	# Called by FilterBank __init__.
	# This class parses the csv file of tissue tags and can return them as a 
	# list of namedtuples.
	class TissueTagLoader:

		def __init__(self, oTissueTagFile):
			self.tTissueTag = namedtuple('tTissueTag', ['origin', 'label', 'tag'])
			self.lTissueTagLoaderTuples = self.LoadTissueTagsFromFile(oTissueTagFile)
		
		def LoadTissueTagsFromFile(self, oTagFile):
			lLoadedTagList = []
			oCSVReader = csv.reader(oTagFile, delimiter = ',')
			bFirstRowSeen = False # boolean for skipping the first row (which is just a header).
			for lRow in oCSVReader:
				if bFirstRowSeen:
					lLoadedTagList.append(self.tTissueTag(*lRow))
				else:
					bFirstRowSeen = True
			return lLoadedTagList
	
		def ReloadTagFile(self, oTagFile):
			self.lTissueTagLoaderTuples = self.LoadTissueTagsFromFile(oTagFile)

		def ReturnTagList(self):
			return self.lTissueTagLoaderTuples

	# a simple class for keeping failure diagnoses externally non-reversable.
	class StatusFailureRecorder:
		def __init__(self):
			self.bSequenceSucceeds = True 
		def SetFailureStatus(self):
			self.bSequenceSucceeds = False
		def ReturnSuccessStatus(self):
			return self.bSequenceSucceeds



	def __init__(self, oTagFile, oLog):
		self.oLog = oLog
		self.oTagFile = oTagFile
		self.oLoadedTags = self.TissueTagLoader(oTagFile)
		self.oIsATissueTag = SequenceAnalysisClasses.IsATissueTag(self.oLoadedTags.ReturnTagList())
		self.oContainsForwardAndReversePrimers = SequenceAnalysisClasses.ContainsForwardAndReversePrimers()
		self.oContainsForwardAndReversePrimers_Complement = SequenceAnalysisClasses.ContainsForwardAndReversePrimers_Complement()
		self.oContainsBothFlankingSequences = SequenceAnalysisClasses.ContainsBothFlankingSequences()
		self.oContainsBothFlankingSequences_Complement = SequenceAnalysisClasses.ContainsBothFlankingSequences_Complement()


	# Returns a ConstantsAndStructures.tSequenceReport namedtuple.
	def RunCompositeAnalysisOnSequence(self, sIDString, sCompleteSequence, sQualiSequence):

		oFailureRecorder = self.StatusFailureRecorder()
		bSequenceIsReversed = None
		lErrorsWithinSequence = []

		# is sequence in forward direction?
		if self.oContainsForwardAndReversePrimers.Ask(sCompleteSequence):
			bSequenceIsReversed = False

			# does sequence contain proper forward and reverse primers?
			if not self.oContainsBothFlankingSequences.Ask(sCompleteSequence):
				oFailureRecorder.SetFailureStatus()
				lErrorsWithinSequence.append('Flanking sequences aren\'t right.')

			# does sequence begin with a proper tissue tag?
			sPossibleTissueTag = self.oContainsForwardAndReversePrimers.ReturnSequencePrependingForwardPrimer(
					sCompleteSequence)
			if not self.oIsATissueTag.Ask(sPossibleTissueTag, bMatchEntireTagOnly=False):
				oFailureRecorder.SetFailureStatus()
				lErrorsWithinSequence.append('Tissue tag does not match tests.')

			# does the insert sequence pass its tests?
			tInsSeqBegEndPos = self.oContainsBothFlankingSequences.ReturnInsSeqBegEndPos(sCompleteSequence)


		# otherwise, is sequence in reverse direction?
		elif self.oContainsForwardAndReversePrimers_Complement.Ask(sCompleteSequence):
			bSequenceIsReversed = True
			if self.oContainsBothFlankingSequences_Complement.Ask(sCompleteSequence):
				pass
				# if begins with tissue tag
				# if insert sequence passes tests

		# otherwise, the sequence is bad.
		else:
			oFailureRecorder.SetFailureStatus()
			lErrorsWithinSequence.append('Contains no forward-reverse primer pairs, normal or complement.')
			

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
			end_pos_ending_seq = 0,
			errors_within_sequence = lErrorsWithinSequence
			)

		return tOutputReport
