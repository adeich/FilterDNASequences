import re
from collections import namedtuple
import Exceptions
import ConstantsAndStructures

class FilterBank:

	class SequenceSubstringCheck:
		"""A basic template class to be subclassed, below.

		The Run() method should always either return a boolean.
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
			self.sRegex = '.*{}.+{}.*'.format(ConstantsAndStructures.sFORWARD_PRIMER,
				ConstantsAndStructures.sREVERSE_PRIMER)
			self.oCompiledRegex = re.compile(self.sRegex)
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
			self.sRegex = '.*{}.+{}.*'.format(ConstantsAndStructures.sREVERSE_PRIMER_COMPLEMENT,
				ConstantsAndStructures.sFORWARD_PRIMER_COMPLEMENT,)
			self.oCompiledRegex = re.compile(self.sRegex)
		def Run(self, sCompleteSequence):
			return bool(self.oCompiledRegex.search(sCompleteSequence))
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

	class BeginsWithTissueTag(SequenceSubstringCheck):
		def __init__(self, oTissueTagFile=None):
			self.oCompiledRegex = None
			self.tTissueTag = namedtuple('tTissueTag', ['origin', 'label', 'tag'])
			self.lTissueTagList = [self.tTissueTag('Hsa-Pat6-Aut', 'huh?', ''),
				self.tTissueTag('Hsa-Pat6-Aut', 'huh?', 'ATTACC'),
				self.tTissueTag('Hsa-Pat6-Aut', 'huh?', 'TGTGTT'),
				self.tTissueTag('Hsa-Pat6-Aut', 'huh?', 'CATATT'),
				self.tTissueTag('Hsa-Pat6-Aut', 'huh?', 'ATACTT'),
				self.tTissueTag('Hsa-Pat6-Aut', 'huh?', 'GTTACA'),
				self.tTissueTag('Hsa-Pat6-Aut', 'huh?', 'CTGTAA'),
				self.tTissueTag('Hsa-Pat6-Aut', 'huh?', 'GAAGAT'),
				self.tTissueTag('Hsa-Pat6-Aut', 'huh?', 'CATACA'),
				self.tTissueTag('Hsa-Pat6-Aut', 'huh?', 'CTGATA'), 
				]
		def Run(self, sCompleteSequence):
			for sTag in self.lTissueTagList:
				if sCompleteSequence[:len(sTag.tag)] == sTag.tag:
					return True
			return False
		def TestSelf(self):
			pass	

	class ContainsBothFlankingSequences(SequenceSubstringCheck):
		def __init__(self):
			self.sRegex = '.+{}.+{}.+'.format(
				ConstantsAndStructures.sBEGINNING_FLANKING_SEQUENCE, 
				ConstantsAndStructures.sENDING_FLANKING_SEQUENCE)
			self.oCompiledRegex = re.compile(self.sRegex)
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
			self.sRegex = '.+{}.+{}.+'.format(
				ConstantsAndStructures.sENDING_FLANKING_SEQUENCE_COMPLEMENT,
				ConstantsAndStructures.sBEGINNING_FLANKING_SEQUENCE_COMPLEMENT)
			self.oCompiledRegex = re.compile(self.sRegex)
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
					



	class ExtractInsSequence(SequenceSubstringCheck):
		def __init__(self):
			self.oCompiledRegex = None
		def Run(self, sCompleteSequence):
			pass
		def TestSelf(self):
			print 'tested self from extract ins sequence'

	class ExtractInsSequence_Complement(SequenceSubstringCheck):
		def __init__(self):
			self.oCompiledRegex = None
		def Run(self, sCompleteSequence):
			pass
		def TestSelf(self):
			pass	



	def __init__(self, oTagFile, oLog):
		self.ClassesDict = {
			'ContainsForwardAndReversePrimers': self.ContainsForwardAndReversePrimers(),
			'ContainsForwardAndReversePrimers_Complement': self.ContainsForwardAndReversePrimers_Complement(),
			'BeginsWithTissueTag': self.BeginsWithTissueTag(),
			'ContainsBothFlankingSequences': self.ContainsBothFlankingSequences(),
			'ContainsBothFlankingSequences_Complement': self.ContainsBothFlankingSequences_Complement(),
			'ExtractInsSequence': self.ExtractInsSequence(),
			'ExtractInsSequence_Complement': self.ExtractInsSequence_Complement()
		}
		self.oLog = oLog
		self.oTagFile = oTagFile
		

	def TestAllClasses(self):
		for oClass in self.ClassesDict:
			try:
				self.ClassesDict[oClass].TestSelf()
			except AttributeError as e:
				oLog.Log("Looks like a testing class is missing the TestSelf() method")
				raise	
				

	# Returns a ConstantsAndStructures.tSequenceReport namedtuple.
	def RunCompositeTestOnSequence(self, sIDString, sCompleteSequence, sQualiSequence):
		tOutputReport = ConstantsAndStructures.tSequenceReport(
			bPasses_filters = True,
			printready_output_string = 'hello',
			input_id = 'hello',
			input_seq = 'hello',
			output_id = 'hello',
			output_seq = 'hello', 
			seq_is_reversed = False,
			start_pos_forward_primer = 6,
			end_pos_forward_primer = 6,
			start_pos_ending_seq = 6,
			end_pos_ending_seq = 6
			)

		return tOutputReport
