import re 
from collections import namedtuple
import Exceptions
import ConstantsAndStructures


class SequenceSubstringCheck:
	"""A basic template class to be subclassed, below.
	
	By convention, these class names generally describe a yes-or-no question;
	accordingly, the Ask() method should always return a boolean.
	"""
	def __init__(self):
		self.sRegex = None
		self.oCompiledRegex = None
	def Ask(self, sDNAString):
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
	def Ask(self, sDNAString):
		return bool(self.oCompiledRegex.search(sDNAString))
#	def TestSelf(self):
#		# sequences which should succeed:
#		CSdU = ConstantsAndStructures.dUnitTestSequences
#		for sUnitTestSequence in [CSdU['everything good # 1'], CSdU['everything good # 2']]:
#			if self.Ask(sUnitTestSequence) == False:
#				raise Exceptions.SequenceCheckUnitTestFail('ContainsForwardAndReversePrimers1')
#		# sequences which should fail:
#		for sUnitTestSequence in [CSdU['everything good # 3 (reversed and with complement)'],
#			CSdU['forward primer has random char inserted']]:
#			if self.Ask(sUnitTestSequence) == True:
#				raise Exceptions.SequenceCheckUnitTestFail('ContainsForwardAndReversePrimers2')
#		print 'Tested self from ContainsForwardAndReversePrimers'
#				

class ContainsForwardAndReversePrimers_Complement(SequenceSubstringCheck):
	def __init__(self):
		self.sRegex = '(.*){}.+{}.*'.format(ConstantsAndStructures.sREVERSE_PRIMER_COMPLEMENT,
			ConstantsAndStructures.sFORWARD_PRIMER_COMPLEMENT,)
		self.oCompiledRegex = re.compile(self.sRegex)
	def Ask(self, sCompleteSequence):
		return bool(self.oCompiledRegex.search(sCompleteSequence))
	# This is for getting the sequence which might be the tissue tag.
#	def TestSelf(self):
#		# sequences which should succeed:
#		CSdU = ConstantsAndStructures.dUnitTestSequences
#		for sUnitTestSequence in [CSdU['everything good # 3 (reversed and with complement)'],
#			CSdU['everything good # 4 (reversed and with complement)']]:
#			if self.Ask(sUnitTestSequence) == False:
#				raise Exceptions.SequenceCheckUnitTestFail('ContainsForwardAndReversePrimers_Complement1')
#		# sequences which should fail:
#		for sUnitTestSequence in [CSdU['everything good # 1'],
#			CSdU['forward primer has random char inserted']]:
#			if self.Ask(sUnitTestSequence) == True:
#				raise Exceptions.SequenceCheckUnitTestFail('ContainsForwardAndReversePrimers_Complement2')
#		print	'tested self from ContainsForwardAndReversePrimers_Complement()'
#

class IsATissueTag(SequenceSubstringCheck):

	def __init__(self, lTissueTagTupleList):
		self.oCompiledRegex = None
		self.dAllPossibleSubtags = self.FillOutAllPossibleSubstrings(lTissueTagTupleList)
		self.lAllCompleteTags = set([sTuple.tag for sTuple in lTissueTagTupleList])

	# this is run on class init, but can be run at other times for testing.
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

	def Ask(self, sPossibleTag, bMatchEntireTagOnly=False):
		bIsTissueTag = False
		if bMatchEntireTagOnly:
			if sPossibleTag in self.lAllCompleteTags:
				bIsTissueTag = True
		else:
			# Second condition, below, ensures that we avoid ambiguous positives:
			# 'AT' could have come from both 'GAT' and 'AAT', so don't count it.
			if sPossibleTag in self.dAllPossibleSubtags and self.dAllPossibleSubtags[sPossibleTag] < 1:
				bIsTissueTag = True
		return bIsTissueTag

	# used for debugging
	def GetTagsAndPartialTags(self):
		return {'tags': self.lAllCompleteTags, 'partial_tags': self.dAllPossibleSubtags}

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

	def Ask(self, sCompleteSequence):
		return bool(self.oCompiledRegex.search(sCompleteSequence))

	def TestSelf(self):
		# sequences which should succeed:
		CSdU = ConstantsAndStructures.dUnitTestSequences
		for sUnitTestSequence in [CSdU['everything good # 1'], CSdU['everything good # 2']]:
			if self.Ask(sUnitTestSequence) == False:
				raise Exceptions.SequenceCheckUnitTestFail('ContainsBothFlankingSequences 1')
		# sequences which should fail:
		for sUnitTestSequence in [CSdU['everything good # 3 (reversed and with complement)'],
			CSdU['contains incorrect beginning flanking sequence']]:
			if self.Ask(sUnitTestSequence) == True:
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
	def Ask(self, sCompleteSequence):
		return bool(self.oCompiledRegex.search(sCompleteSequence))
	def TestSelf(self):
		CSdU = ConstantsAndStructures.dUnitTestSequences
		# sequences which should succeed:
		for sUnitTestSequence in [CSdU['everything good # 3 (reversed and with complement)'],
			CSdU['everything good # 4 (reversed and with complement)']]:
			if self.Ask(sUnitTestSequence) == False:
				raise Exceptions.SequenceCheckUnitTestFail('ContainsBothFlankingSequences_Complement 1')
		# sequences which should fail:
		for sUnitTestSequence in [CSdU['everything good # 1'],
			CSdU['forward primer has random char inserted']]:
			if self.Ask(sUnitTestSequence) == True:
				raise Exceptions.SequenceCheckUnitTestFail('ContainsBothFlankingSequences_Complement 2')

		print 'Tested self from ContainsBothFlankingSequences_Complement'
				

class InsertSequencePassesTests(SequenceSubstringCheck):
	def __init__(self):
		pass
	def Ask(self, sCompleteSequence):
		# Is the length of the string a multiple of 3?
		return (len(sCompleteSequence) % 3 == 0)

class QualiSequenceIsValid(SequenceSubstringCheck):
	def __init__(self):
		pass	
	def QualiSequenceIsValid(self, sQualiSequence):
		# doesn't have any characters besides integers and spaces.
		# doesn't have any doubled spaces.		
		return None
	def Ask(self, sQualiSequence):
		pass


class QualiInsertSequenceAllAboveThreshold(SequenceSubstringCheck):

	def __init__(self):
		pass

	def GetQualiInsertSubsequence(sQualiSequence, tInsSeqBegEndPos):
		QualiSequenceIsValid(sQualiSequence) # look for errors
		sSplitString = [int(sDigits) for sDigits in sQualiSequence.split(' ')]
		return sSplitString[tInsSeqBegEndPos[0]:tInsSeqBegEndPos[1]]		

	def Ask(self, sQualiSequence, tInsSeqBegEndPos):
		sQualiInsertSequence = GetQualiInsertSubsequence(sQualiSequence, tInsSeqBegEndPos)
		for iInteger in sInsertSequence:
			if iInteger < 20:
				return False
		return True


