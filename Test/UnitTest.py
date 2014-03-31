import sys, os
testdir = os.path.dirname(__file__)
srcdir = '../FilterDNASequences'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))

import unittest
import FilterBank
import CustomLog
import Main
import ConstantsAndStructures as CAS
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class AnalysisClassTestCase(unittest.TestCase):

	def setUp(self):

		# Create a test tag file using a local list of tags.
		def CreateTestTagFile(sFilePath, lTagList):
			with open(sFilePath, 'w') as oFile:
				oFile.write('ORIGIN,LABEL,TAG\n')
				for sTag in lTagList:
					oFile.write('UnitTest.py,no_label,{}\n'.format(sTag))

		self.lTestTags = ['ABCDE', 'BCDEF', 'CDEFG', 'DDEFG']
		self.sTestTagFileName = 'test_input_data/test_tags_from_UnitTest.csv'
		CreateTestTagFile(self.sTestTagFileName, self.lTestTags)

		self.oTagFile = open(self.sTestTagFileName, 'r')
		self.oLogFile = open('test_output_dump/UnitTest.py-created-log1.txt', 'a') 
		self.oLog = CustomLog.Log(self.oLogFile)
		
		# Initialize the FilterBank object using the freshly-created 
		# test tag file.
		self.oFilterBank = FilterBank.FilterBank(self.oTagFile, self.oLog)
		self.CSdU = CAS.dUnitTestSequences

	def tearDown(self):
		self.oTagFile.close()
		self.oLogFile.close()

	def test_ContainsForwardAndReversePrimers(self):
		Ask = self.oFilterBank.oContainsForwardAndReversePrimers.Ask
		self.assertTrue(Ask(self.CSdU['correct_normal_form_1']), '')
		self.assertTrue(Ask(self.CSdU['correct_normal_form_2']), '')
		self.assertFalse(Ask(self.CSdU['correct_reversed_complement_3']), '')
		self.assertFalse(Ask(self.CSdU['correct_reversed_complement_4']), '')
		self.assertFalse(Ask(self.CSdU['forward_primer_has_random_char_inserted']), '')

	def test_ContainsForwardAndReversePrimers_Complement(self):
		Ask = self.oFilterBank.oContainsForwardAndReversePrimers_Complement.Ask
		self.assertFalse(Ask(self.CSdU['correct_normal_form_1']), '')
		self.assertFalse(Ask(self.CSdU['correct_normal_form_2']), '')
		self.assertTrue(Ask(self.CSdU['correct_reversed_complement_3']), '')
		self.assertTrue(Ask(self.CSdU['correct_reversed_complement_4']), '')

	def test_IsATissueTag(self):
		Ask = self.oFilterBank.oIsATissueTag.Ask	

		# Generate some partial tags 
		lSomePartialTags = [sTag[1:] for sTag in self.lTestTags]
		lSomePartialTags += [sTag[2:] for sTag in self.lTestTags]

		# Test that *full* tags are recognized.
		for sTissueTag in self.lTestTags:
			self.assertTrue(Ask(sTissueTag, bMatchEntireTagOnly = True), 
				'The sequence \'{}\' should be a member of the tags {}'.format(sTissueTag, 
					self.oFilterBank.oIsATissueTag.GetTagsAndPartialTags()['tags']))
			self.assertTrue(Ask(sTissueTag, bMatchEntireTagOnly = False)) 
		for sPartialTag in lSomePartialTags:
			self.assertFalse(Ask(sPartialTag, bMatchEntireTagOnly = True))
			
		# Test that *partial* tags are recognized.
		for sPartialTag in lSomePartialTags:
			if self.oFilterBank.oIsATissueTag.GetTagsAndPartialTags()['partial_tags'][sPartialTag] < 1:
				self.assertTrue(Ask(sPartialTag, bMatchEntireTagOnly = False),
					'The tag \'{}\' should be a member of the tags {}'.format(sPartialTag, 
						self.oFilterBank.oIsATissueTag.GetTagsAndPartialTags()['partial_tags']))
			else:
				self.assertFalse(Ask(sPartialTag, bMatchEntireTagOnly = False),
					'The tag \'{}\' should be a member of the tags {}'.format(sPartialTag, 
						self.oFilterBank.oIsATissueTag.GetTagsAndPartialTags()['partial_tags']))
			self.assertFalse(Ask(sPartialTag, bMatchEntireTagOnly = True))

		# Test the AddTag function.
		sRandomTagName = 'Whoa!TotallyRandomTag!'
		self.oFilterBank.oIsATissueTag.AddTag(sRandomTagName)
		self.assertTrue(Ask(sRandomTagName, bMatchEntireTagOnly = True))

	
	def test_IsATissueTag_Complement(self):
		Ask = self.oFilterBank.oIsATissueTag_Complement.Ask	

		def ConvertToComplement(sInputSeq):
			return str(Seq(sInputSeq, IUPAC.unambiguous_dna).reverse_complement())

		# Generate some partial tags (of the complement form)
		lSomePartialTags = [ConvertToComplement(sTag)[:3] for sTag in self.lTestTags]
		lSomePartialTags += [ConvertToComplement(sTag)[:4] for sTag in self.lTestTags]

		# Test that *full* tags are recognized.
		for sTissueTag in [ConvertToComplement(sTag) for sTag in self.lTestTags]:
			self.assertTrue(Ask(sTissueTag, bMatchEntireTagOnly = True), 
				'The sequence \'{}\' should be a member of the tags {}'.format(sTissueTag, 
					self.oFilterBank.oIsATissueTag_Complement.GetTagsAndPartialTags()['tags']))
			self.assertTrue(Ask(sTissueTag, bMatchEntireTagOnly = False)) 
		for sPartialTag in lSomePartialTags:
			self.assertFalse(Ask(sPartialTag, bMatchEntireTagOnly = True))
			
		# Test that *partial* tags are recognized.
		for sPartialTag in lSomePartialTags:
			if self.oFilterBank.oIsATissueTag_Complement.GetTagsAndPartialTags()['partial_tags'][sPartialTag] < 1:
				self.assertTrue(Ask(sPartialTag, bMatchEntireTagOnly = False),
					'The tag \'{}\' should be a member of the tags {}'.format(sPartialTag, 
						self.oFilterBank.oIsATissueTag_Complement.GetTagsAndPartialTags()['partial_tags']))
			else:
				self.assertFalse(Ask(sPartialTag, bMatchEntireTagOnly = False),
					'The tag \'{}\' should be a member of the tags {}'.format(sPartialTag, 
						self.oFilterBank.oIsATissueTag_Complement.GetTagsAndPartialTags()['partial_tags']))
			self.assertFalse(Ask(sPartialTag, bMatchEntireTagOnly = True))

		# Test the AddTag function.
		sRandomTagName = 'Whoa!TotallyRandomTag!'
		self.oFilterBank.oIsATissueTag_Complement.AddTag(sRandomTagName)
		self.assertTrue(Ask(sRandomTagName, bMatchEntireTagOnly = True))
		
	def test_ContainsBothFlankingSequences(self):
		sCorrectExample = 'XXX{}XXX{}XXX'.format(
			self.oFilterBank.SequenceConstants.sBEGINNING_FLANKING_SEQUENCE,
			self.oFilterBank.SequenceConstants.sENDING_FLANKING_SEQUENCE)
				
		# lacks beginning flanking sequence.
		sIncorrectExample = 'XXXXXX{}XXX'.format(
			self.oFilterBank.SequenceConstants.sENDING_FLANKING_SEQUENCE)

		Ask = self.oFilterBank.oContainsBothFlankingSequences.Ask

		self.assertTrue(Ask(sCorrectExample))
		self.assertFalse(Ask(sIncorrectExample))
		

	def test_ContainsBothFlankingSequences_Complement(self):
		sCorrectExample = 'XXX{}XXX{}XXX'.format(
			self.oFilterBank.SequenceConstants.sENDING_FLANKING_SEQUENCE_COMPLEMENT,
			self.oFilterBank.SequenceConstants.sBEGINNING_FLANKING_SEQUENCE_COMPLEMENT)
				
		# lacks beginning flanking sequence.
		sIncorrectExample = 'XXXXXX{}XXX'.format(
			self.oFilterBank.SequenceConstants.sENDING_FLANKING_SEQUENCE)

		Ask = self.oFilterBank.oContainsBothFlankingSequences_Complement.Ask

		self.assertTrue(Ask(sCorrectExample))
		self.assertFalse(Ask(sIncorrectExample))

	def test_QualiInsertSequenceAllAboveThreshold(self):
		Ask = self.oFilterBank.oQualiInsertSequenceAllAboveThreshold.Ask

		sOuterStuff = 'XXX'
		sInsertSequence = 'INSERTSEQUENCE'
		sCorrectExample = '{stuff1}{beginning_flank}{insert}{ending_flank}{stuff2}'.format(
			stuff1 = sOuterStuff,
			beginning_flank = self.oFilterBank.SequenceConstants.sBEGINNING_FLANKING_SEQUENCE,
			insert = sInsertSequence,
			ending_flank = self.oFilterBank.SequenceConstants.sENDING_FLANKING_SEQUENCE,
			stuff2 = sOuterStuff
			)

		tInsSeqBegEndPos = self.oFilterBank.oContainsBothFlankingSequences.ReturnInsSeqBegEndPos(sCorrectExample)

		sCorrectQualiSequence = '{stuff1} {beginning_flank} {insert} {ending_flank} {stuff2}'.format(
			stuff1 = ' '.join(['10' for i in range(len(sOuterStuff))]),
			beginning_flank = ' '.join(['05' for i in range(len(self.oFilterBank.SequenceConstants.sBEGINNING_FLANKING_SEQUENCE))]),
			insert = ' '.join(['25' for i in range(len(sInsertSequence))]),
			ending_flank = ' '.join(['05' for i in range(len(self.oFilterBank.SequenceConstants.sENDING_FLANKING_SEQUENCE))]),
			stuff2 = ' '.join(['10' for i in range(len(sOuterStuff))])
			)

		self.assertTrue(Ask(sCorrectQualiSequence, tInsSeqBegEndPos))

		#                    0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
		sQualiSnippetGood = '10 11 12 20 20 20 25 25 25 20 20 20 25 25 25 20 20 12 11 10'
		sQualiSnippetBad =  '10 11 12 20 20 20 19 00 25 20 20 20 25 25 25 20 20 12 11 10'

		self.assertTrue(Ask(sQualiSnippetGood, (3, 17))) 
		# False because includes an integer lower than 20 (12)
		self.assertFalse(Ask(sQualiSnippetGood, (2, 17))) 
		self.assertFalse(Ask(sQualiSnippetBad, (3, 17)))

	def test_QualiInsertSequenceAllAboveThreshold_Complement(self):
		Ask = self.oFilterBank.oQualiInsertSequenceAllAboveThreshold.Ask

		sOuterStuff = 'XXX'
		sInsertSequence = 'INSERTSEQUENCE'
		sCorrectExample = '{stuff1}{ending_flank_complement}{insert}{beginning_flank_complement}{stuff2}'.format(
			stuff1 = sOuterStuff,
			ending_flank_complement = self.oFilterBank.SequenceConstants.sENDING_FLANKING_SEQUENCE_COMPLEMENT,
			insert = sInsertSequence,
			beginning_flank_complement = self.oFilterBank.SequenceConstants.sBEGINNING_FLANKING_SEQUENCE_COMPLEMENT,
			stuff2 = sOuterStuff
			)

		tInsSeqBegEndPos = self.oFilterBank.oContainsBothFlankingSequences_Complement.ReturnInsSeqBegEndPos(sCorrectExample)

		sCorrectQualiSequence = '{stuff1} {beginning_flank} {insert} {ending_flank} {stuff2}'.format(
			stuff1 = ' '.join(['10' for i in range(len(sOuterStuff))]),
			beginning_flank = ' '.join(['05' for i in range(len(self.oFilterBank.SequenceConstants.sENDING_FLANKING_SEQUENCE_COMPLEMENT))]),
			insert = ' '.join(['25' for i in range(len(sInsertSequence))]),
			ending_flank = ' '.join(['05' for i in range(len(self.oFilterBank.SequenceConstants.sBEGINNING_FLANKING_SEQUENCE_COMPLEMENT))]),
			stuff2 = ' '.join(['10' for i in range(len(sOuterStuff))])
			)

		self.assertTrue(Ask(sCorrectQualiSequence, tInsSeqBegEndPos))

		#                    0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
		sQualiSnippetGood = '10 11 12 20 20 20 25 25 25 20 20 20 25 25 25 20 20 12 11 10'
		sQualiSnippetBad =  '10 11 12 20 20 20 19 00 25 20 20 20 25 25 25 20 20 12 11 10'

		self.assertTrue(Ask(sQualiSnippetGood, (3, 17))) 
		# False because includes an integer lower than 20 (12)
		self.assertFalse(Ask(sQualiSnippetGood, (2, 17))) 
		self.assertFalse(Ask(sQualiSnippetBad, (3, 17)))



class FilterBankTestCase(unittest.TestCase):
	
	def setUp(self):
		self.oTagFile = open('test_input_data/test_tags.csv', 'r') 
		self.oLogFile = open('test_output_dump/UnitTest.py-created-log2.txt', 'a') 
		self.oLog = CustomLog.Log(self.oLogFile)
		self.oFilterBank = FilterBank.FilterBank(self.oTagFile, self.oLog)

	def tearDown(self):
		self.oTagFile.close()
		self.oLog.AddAnalysisReport()
		self.oLogFile.close()

	def	test_CorrectSequencesSucceed(self):
		self.oFilterBank.oIsATissueTag.AddTag('AGTCAT')
		tReport1 = self.oFilterBank.RunCompositeAnalysisOnSequence(
			sIDString = None,
			sCompleteSequence = CAS.dUnitTestSequences['correct_normal_form_1'],
			sQualiSequence = None,
			bSuppressQualiChecks = True)

		self.oLog.IngestReportAndLog(tReport1)

		self.assertTrue(tReport1.bPasses_filters, 'custom errors: {}'.format(
			tReport1.errors_within_sequence))


class MainProgramTestCase(unittest.TestCase):

	# Not much precise to test here except that no errors are thrown.
	# This test is also useful in manually inspecting the log it creates.
	def test_Main(self):
		Main.Main(
 			sOutputSequenceFileName='test_output_dump/test_created_generated_output_seq.txt',
			sFastaFileName='test_input_data/300LinesTest.fna',
			sQualiFileName='test_input_data/300LinesTest.qual',
			sLogFileName='test_output_dump/UnitTest.py-created-log3.txt',
			sTagFileName='test_input_data/test_tags.csv')



if __name__ == '__main__':
	unittest.main()
