import re
from collections import namedtuple
import Exceptions
import ConstantsAndStructures
import csv

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


	class SequenceSubstringCheck:
		"""A basic template class to be subclassed, below.

		The Run() method should always return a boolean.
		"""
		def __init__(self):
			self.sRegex = None
			self.oCompiledRegex = None
		def Run(self, sDNAString):
			pass
		def TestSelf(self):
			pass


	class ContainsForwardAndReversePrimers(SequenceSubstringCheck):
		def __init__(self):
			self.sRegex = '(.*){}.+{}.*'.format(ConstantsAndStructures.sFORWARD_PRIMER,
				ConstantsAndStructures.sREVERSE_PRIMER)
			self.oCompiledRegex = re.compile(self.sRegex)
		def ReturnSequencePrependingForwardPrimer(self, sCompleteSequence):
			return self.oCompiledRegex.search(sCompleteSequence).group(1)
		def Run(self, sDNAString):
			return bool(self.oCompiledRegex.search(sDNAString))
		def TestSelf(self):
			# sequences which should succeed:
			CSdU = ConstantsAndStructures.dUnitTestSequences
			for sUnitTestSequence in [CSdU['everything good # 1'], CSdU['everything good # 2']]:
				if self.Run(sUnitTestSequence) == False:
					raise Exceptions.SequenceCheckUnitTestFail('ContainsForwardAndReversePrimers1')
			# sequences which should fail:
			for sUnitTestSequence in [CSdU['everything good # 3 (reversed and with complement)'],
				CSdU['forward primer has random char inserted']]:
				if self.Run(sUnitTestSequence) == True:
					raise Exceptions.SequenceCheckUnitTestFail('ContainsForwardAndReversePrimers2')
			print 'Tested self from ContainsForwardAndReversePrimers'
					

	class ContainsForwardAndReversePrimers_Complement(SequenceSubstringCheck):
		def __init__(self):
			self.sRegex = '(.*){}.+{}.*'.format(ConstantsAndStructures.sREVERSE_PRIMER_COMPLEMENT,
				ConstantsAndStructures.sFORWARD_PRIMER_COMPLEMENT,)
			self.oCompiledRegex = re.compile(self.sRegex)
		def Run(self, sCompleteSequence):
			return bool(self.oCompiledRegex.search(sCompleteSequence))
		# This is for getting the sequence which might be the tissue tag.
		def TestSelf(self):
			# sequences which should succeed:
			CSdU = ConstantsAndStructures.dUnitTestSequences
			for sUnitTestSequence in [CSdU['everything good # 3 (reversed and with complement)'],
				CSdU['everything good # 4 (reversed and with complement)']]:
				if self.Run(sUnitTestSequence) == False:
					raise Exceptions.SequenceCheckUnitTestFail('ContainsForwardAndReversePrimers_Complement1')
			# sequences which should fail:
			for sUnitTestSequence in [CSdU['everything good # 1'],
				CSdU['forward primer has random char inserted']]:
				if self.Run(sUnitTestSequence) == True:
					raise Exceptions.SequenceCheckUnitTestFail('ContainsForwardAndReversePrimers_Complement2')
			print	'tested self from ContainsForwardAndReversePrimers_Complement()'


	class IsATissueTag(SequenceSubstringCheck):

		def __init__(self, lTissueTagTupleList):
			self.oCompiledRegex = None
			self.dAllPossibleSubtags = self.FillOutAllPossibleSubstrings(lTissueTagTupleList)
			self.lAllCompleteTags = set([sTuple.tag for sTuple in lTissueTagTupleList])

		# this is mostly meant to be run on the class init.
		def FillOutAllPossibleSubstrings(self, lTissueTagTupleList):
			dAllPossibleSubtags = {}
			for sTag in [sTuple.tag for sTuple in lTissueTagTupleList]:
				# Only include subtags down to length 3, because we're not accepting any shorter.
				for iStartPos in range(0, len(sTag) - 2):
					sSubTag = sTag[iStartPos:]
					if sSubTag in dAllPossibleSubtags:
						dAllPossibleSubtags[sSubTag] += 1
					else:
						dAllPossibleSubtags[sSubTag] = 0
			return dAllPossibleSubtags

		def Run(self, sSequence, bMatchEntireTagOnly=False):
			bReturnValue = False
			if bMatchEntireTagOnly:
				if sSequence in self.lAllCompleteTags:
					bReturnValue = True
			else:
				# Second condition ensures that we avoid ambiguous positives:
				# 'AT' could have come from both 'GAT' and 'AAT', so don't count it.
				if sSequence in self.dAllPossibleSubtags and self.dAllPossibleSubtags[sSequence] < 1:
					bReturnValue = True
			return bReturnValue

		def TestSelf(self):
			pass	


	class ContainsBothFlankingSequences(SequenceSubstringCheck):

		def __init__(self):
			self.sRegex = '.+{}(.+){}.+'.format(
				ConstantsAndStructures.sBEGINNING_FLANKING_SEQUENCE, 
				ConstantsAndStructures.sENDING_FLANKING_SEQUENCE)
			self.oCompiledRegex = re.compile(self.sRegex)

		def ReturnInsertSequence(self, sCompleteSequence):
			sInsertSequence = None
			try:
				sInsertSequence = self.oCompiledRegex.search(sCompleteSequence).group(1)
			except AttributeError as e:
				sInsertSequence = 'no insert sequence pattern matched.'
			return sInsertSequence

		def ReturnInsSeqBegEndPos(self, sCompleteSequence):
			tInsSeqBegEndPos = None
			try:
				tInsSeqBegEndPos = self.oCompiledRegex.search(sCompleteSequence).span(1)
			except AttributeError as e:
				pass
			return tInsSeqBegEndPos

		def Run(self, sCompleteSequence):
			return bool(self.oCompiledRegex.search(sCompleteSequence))

		def TestSelf(self):
			# sequences which should succeed:
			CSdU = ConstantsAndStructures.dUnitTestSequences
			for sUnitTestSequence in [CSdU['everything good # 1'], CSdU['everything good # 2']]:
				if self.Run(sUnitTestSequence) == False:
					raise Exceptions.SequenceCheckUnitTestFail('ContainsBothFlankingSequences 1')
			# sequences which should fail:
			for sUnitTestSequence in [CSdU['everything good # 3 (reversed and with complement)'],
				CSdU['contains incorrect beginning flanking sequence']]:
				if self.Run(sUnitTestSequence) == True:
					raise Exceptions.SequenceCheckUnitTestFail('ContainsBothFlankingSequences 2')
			print 'Tested self from ContainsBothFlankingSequences'
					

	class ContainsBothFlankingSequences_Complement(SequenceSubstringCheck):
		def __init__(self):
			self.sRegex = '.+{}(.+){}.+'.format(
				ConstantsAndStructures.sENDING_FLANKING_SEQUENCE_COMPLEMENT,
				ConstantsAndStructures.sBEGINNING_FLANKING_SEQUENCE_COMPLEMENT)
			self.oCompiledRegex = re.compile(self.sRegex)
		def ReturnInsertSequence(self, sCompleteSequence):
			sInsertSequence = None
			try:
				sInsertSequence = self.oCompiledRegex.search(sCompleteSequence).group(1)
			except AttributeError as e:
				sInsertSequence = 'no insert sequence pattern matched.'
			return sInsertSequence
		def Run(self, sCompleteSequence):
			return bool(self.oCompiledRegex.search(sCompleteSequence))
		def TestSelf(self):
			CSdU = ConstantsAndStructures.dUnitTestSequences
			# sequences which should succeed:
			for sUnitTestSequence in [CSdU['everything good # 3 (reversed and with complement)'],
				CSdU['everything good # 4 (reversed and with complement)']]:
				if self.Run(sUnitTestSequence) == False:
					raise Exceptions.SequenceCheckUnitTestFail('ContainsBothFlankingSequences_Complement 1')
			# sequences which should fail:
			for sUnitTestSequence in [CSdU['everything good # 1'],
				CSdU['forward primer has random char inserted']]:
				if self.Run(sUnitTestSequence) == True:
					raise Exceptions.SequenceCheckUnitTestFail('ContainsBothFlankingSequences_Complement 2')
	
			print 'Tested self from ContainsBothFlankingSequences_Complement'
					
	
	class InsertSequencePassesTests(SequenceSubstringCheck):
		def __init__(self):
			pass
		def Run(self, sCompleteSequence, sQualiSequence, tInsSeqBegEndPos):
			pass

	def __init__(self, oTagFile, oLog):
		self.oLog = oLog
		self.oTagFile = oTagFile
		self.oLoadedTags = self.LoadedTissueTags(oTagFile)
		self.ClassesDict = {
			'ContainsForwardAndReversePrimers': self.ContainsForwardAndReversePrimers(),
			'ContainsForwardAndReversePrimers_Complement': self.ContainsForwardAndReversePrimers_Complement(),
			'IsATissueTag': self.IsATissueTag(self.oLoadedTags.ReturnTagList()),
			'ContainsBothFlankingSequences': self.ContainsBothFlankingSequences(),
			'ContainsBothFlankingSequences_Complement': self.ContainsBothFlankingSequences_Complement()
		}
		

	def TestAllClasses(self):
		for oClass in self.ClassesDict:
			try:
				self.ClassesDict[oClass].TestSelf()
			except AttributeError as e:
				oLog.Log("Looks like a testing class is missing the TestSelf() method")
				raise	
				

	# Returns a ConstantsAndStructures.tSequenceReport namedtuple.
	def RunCompositeTestOnSequence(self, sIDString, sCompleteSequence, sQualiSequence):

		# a simple class for keeping failure diagnoses externally non-reversable.
		class StatusFailureRecorder:
			def __init__(self):
				self.bSequenceSucceeds = True 
			def SetFailureStatus(self):
				self.bSequenceSucceeds = False
			def ReturnSuccessStatus(self):
				return self.bSequenceSucceeds

		oFailureRecorder = StatusFailureRecorder()
		lErrorsWithSequence = []

		# is sequence in forward direction?
		if self.ClassesDict['ContainsForwardAndReversePrimers'].Run(sCompleteSequence):

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
			if self.ClassesDict['ContainsBothFlankingSequences_Complement'].Run(sCompleteSequence):
				pass
				# if begins with tissue tag
				# if insert sequence passes tests

		# otherwise, the sequence is bad.
		else:
			oFailureRecorder.SetFailureStatus()
			lErrorsWithSequence.append('Contains no forward-reverse primer pairs, normal or complement.')
			

		tOutputReport = ConstantsAndStructures.tSequenceReport(
			bPasses_filters = True,
			printready_output_string = 'hello',
			input_id = sIDString,
			input_seq = sCompleteSequence,
			output_id = sIDString,
			output_seq = sCompleteSequence,
			seq_is_reversed = False,
			start_pos_forward_primer = 6,
			end_pos_forward_primer = 6,
			start_pos_ending_seq = 6,
			end_pos_ending_seq = 6
			)

		return tOutputReport
