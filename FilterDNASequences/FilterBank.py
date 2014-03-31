import re
from collections import namedtuple
import Exceptions
import ConstantsAndStructures
import csv
import SequenceAnalysisClasses 
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


class FilterBank:

	class TissueTagLoader:
		"""This subclass parses the csv file of tissue tags and can return them as a 
		list of namedtuples.
		Called by FilterBank __init__().
		"""

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

	class StatusFailureRecorder:
		""" a simple subclass for keeping failure diagnoses externally non-reversable. """
		def __init__(self):
			self.bSequenceSucceeds = True 
			self.lErrorsWithinSequence = []
		def AddFailureReason(self, sErrorMessage=None):
			self.bSequenceSucceeds = False
			if sErrorMessage:
				self.lErrorsWithinSequence.append(sErrorMessage)
		def ReturnSuccessStatus(self):
			return self.bSequenceSucceeds
		def ReturnErrorList(self):
			return self.lErrorsWithinSequence


	class SequenceConstantsFacade:
		""" Stores sequence constants such as Forward Primer and Forward Primer Complement.
		This object on init loads in half its sequences from the constants file
		and then adds the sequences' complements, generated on the fly with the Biopython library.  
		"""
		
		def __init__(self, 
			sFORWARD_PRIMER, 
			sREVERSE_PRIMER, 
			sBEGINNING_FLANKING_SEQUENCE, 	
			sENDING_FLANKING_SEQUENCE):

			# Add constant sequences.
			self.sFORWARD_PRIMER = sFORWARD_PRIMER
			self.sREVERSE_PRIMER = sREVERSE_PRIMER
			self.sBEGINNING_FLANKING_SEQUENCE = sBEGINNING_FLANKING_SEQUENCE
			self.sENDING_FLANKING_SEQUENCE = sENDING_FLANKING_SEQUENCE

			# Generate and add complement sequences.
			self.sFORWARD_PRIMER_COMPLEMENT = str(Seq(
				self.sFORWARD_PRIMER, IUPAC.unambiguous_dna).reverse_complement())
			self.sREVERSE_PRIMER_COMPLEMENT = str(Seq(
				self.sREVERSE_PRIMER, IUPAC.unambiguous_dna).reverse_complement())
			self.sBEGINNING_FLANKING_SEQUENCE_COMPLEMENT = str(Seq(
				self.sBEGINNING_FLANKING_SEQUENCE, IUPAC.unambiguous_dna).reverse_complement())
			self.sENDING_FLANKING_SEQUENCE_COMPLEMENT = str(Seq(
				self.sENDING_FLANKING_SEQUENCE, IUPAC.unambiguous_dna).reverse_complement())



	def __init__(self, oTagFile, oLog):
		self.SequenceConstants = self.SequenceConstantsFacade(
			sFORWARD_PRIMER = ConstantsAndStructures.sFORWARD_PRIMER,
			sREVERSE_PRIMER = ConstantsAndStructures.sREVERSE_PRIMER,
			sBEGINNING_FLANKING_SEQUENCE = ConstantsAndStructures.sBEGINNING_FLANKING_SEQUENCE,
			sENDING_FLANKING_SEQUENCE = ConstantsAndStructures.sENDING_FLANKING_SEQUENCE)

		self.oLog = oLog
		self.oTagFile = oTagFile
		self.oLoadedTags = self.TissueTagLoader(oTagFile)
		self.oIsATissueTag = SequenceAnalysisClasses.IsATissueTag(self.oLoadedTags.ReturnTagList())
		self.oContainsForwardAndReversePrimers = SequenceAnalysisClasses.ContainsForwardAndReversePrimers(
			self.SequenceConstants)
		self.oContainsForwardAndReversePrimers_Complement = SequenceAnalysisClasses.ContainsForwardAndReversePrimers_Complement(
			self.SequenceConstants)
		self.oContainsBothFlankingSequences = SequenceAnalysisClasses.ContainsBothFlankingSequences(
			self.SequenceConstants)
		self.oContainsBothFlankingSequences_Complement = SequenceAnalysisClasses.ContainsBothFlankingSequences_Complement(
			self.SequenceConstants)
		self.oIsATissueTag_Complement = SequenceAnalysisClasses.IsATissueTag_Complement(
			self.oLoadedTags.ReturnTagList())
		self.oQualiInsertSequenceAllAboveThreshold = SequenceAnalysisClasses.QualiInsertSequenceAllAboveThreshold()


	def RunCompositeAnalysisOnSequence(self, sIDString, sCompleteSequence, 
		sQualiSequence, bSuppressQualiChecks=False):
		""" Runs all checks on a single sequence and IDstring combo. Called by Main.py.
		Returns a ConstantsAndStructures.tSequenceReport namedtuple.
		"""

		oFailureRecorder = self.StatusFailureRecorder()
		bSequenceIsReversed = None


		### SEQUENCE IS IN FORWARD DIRECTION WITH CORRECT FORWARD AND REVERSE PRIMERS. ### 
		if self.oContainsForwardAndReversePrimers.Ask(sCompleteSequence):
			bSequenceIsReversed = False

			# does sequence contain proper flanking sequences?
			if not self.oContainsBothFlankingSequences.Ask(sCompleteSequence):
				oFailureRecorder.AddFailureReason(sErrorMessage='Flanking sequences aren\'t correct.')

			# does sequence begin with a proper tissue tag?
			sPossibleTissueTag = self.oContainsForwardAndReversePrimers.ReturnSequencePrependingForwardPrimer(
					sCompleteSequence)
			if not self.oIsATissueTag.Ask(sPossibleTissueTag, bMatchEntireTagOnly=False):
				oFailureRecorder.AddFailureReason(sErrorMessage='Tissue tag does not match tests.')

			# does the insert sequence pass its tests?
			if not bSuppressQualiChecks:
				tInsSeqBegEndPos = self.oContainsBothFlankingSequences.ReturnInsSeqBegEndPos(sCompleteSequence)
				if not self.oQualiInsertSequenceAllAboveThreshold.Ask(sQualiSequence, tInsSeqBegEndPos):
					oFailureRecorder.AddFailureReason('Quali sequence does not pass tests.')


		### SEQUENCE IS IN REVERSE DIRECTION AND WITH CORRECT COMPLEMENT PRIMERS. ###
		elif self.oContainsForwardAndReversePrimers_Complement.Ask(sCompleteSequence):
			bSequenceIsReversed = True
				
			# does sequence contain proper flanking sequences?
			if not self.oContainsBothFlankingSequences_Complement.Ask(sCompleteSequence):
				oFailureRecorder.AddFailureReason(sErrorMessage='Flanking sequences (complement) aren\'t correct.')

			# does sequence *end* with a proper (complement) tissue tag?
			sPossibleTissueTag = self.oContainsForwardAndReversePrimers_Complement.ReturnSequenceAppendingForwardPrimer(
				sCompleteSequence)
			if not self.oIsATissueTag_Complement.Ask(sPossibleTissueTag, bMatchEntireTagOnly=False):
				oFailureRecorder.AddFailureReason(sErrorMessage='Tissue tag (complement) does not match tests.')

			# does the insert sequence pass its tests?
			if not bSuppressQualiChecks:
				tInsSeqBegEndPos = self.oContainsBothFlankingSequences_Complement.ReturnInsSeqBegEndPos(sCompleteSequence)
				if not self.oQualiInsertSequenceAllAboveThreshold.Ask(sQualiSequence, tInsSeqBegEndPos):
					oFailureRecorder.AddFailureReason('Quali sequence (complement) does not pass tests.')



		### ELSE SEQUENCE IS AUTOMATICALLY INCORRECT. ###
		else:
			oFailureRecorder.AddFailureReason('Contains no forward-reverse primer pairs, normal or complement.')
			

		tOutputReport = ConstantsAndStructures.tSequenceReport(
			bPasses_filters = oFailureRecorder.ReturnSuccessStatus(),
			printready_output_string = 'hello',
			input_id = sIDString,
			input_seq = sCompleteSequence,
			input_quali_seq = sQualiSequence,
			output_id = sIDString,
			output_seq = sCompleteSequence,
			seq_is_reversed = bSequenceIsReversed,
			start_pos_forward_primer = None,
			end_pos_forward_primer = None,
			start_pos_ending_seq = 0,
			end_pos_ending_seq = 0,
			errors_within_sequence = oFailureRecorder.ReturnErrorList()
			)

		return tOutputReport
