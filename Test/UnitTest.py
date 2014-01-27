import sys, os
testdir = os.path.dirname(__file__)
srcdir = '../FilterDNASequences'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))

import unittest
import FilterBank
import CustomLog
import ConstantsAndStructures as CAS

class AnalysisClassTestCase(unittest.TestCase):

	def setUp(self):

		# Create a test tag file using a local list of tags.
		def CreateTestTagFile(sFilePath, lTagList):
			with open(sFilePath, 'w') as oFile:
				oFile.write('ORIGIN,LABEL,TAG\n')
				for sTag in lTagList:
					oFile.write('UnitTest.py,no_label,{}\n'.format(sTag))

		self.lTestTags = ['ABCDE', 'BCDEF', 'CDEFG']
		self.sTestTagFileName = 'test_data/test_tags.csv'
		CreateTestTagFile(self.sTestTagFileName, self.lTestTags)

		self.oTagFile = open(self.sTestTagFileName, 'r')
		self.oLogFile = open('test_output_dump/log.txt', 'w') 
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
				'The sequence \'{}\' should be succeeding for the tags {}'.format(sTissueTag, 
					self.oFilterBank.oIsATissueTag.GetTagsAndPartialTags()['tags']))
			self.assertTrue(Ask(sTissueTag, bMatchEntireTagOnly = False)) 
		for sPartialTag in lSomePartialTags:
			self.assertFalse(Ask(sPartialTag, bMatchEntireTagOnly = True))
			
		# Test that *partial* tags are recognized.
		for sPartialTag in lSomePartialTags:
			self.assertTrue(Ask(sPartialTag, bMatchEntireTagOnly = False),
				'The tag \'{}\' should be succeeding for the tags {}'.format(sPartialTag, 
					self.oFilterBank.oIsATissueTag.GetTagsAndPartialTags()['partial_tags']))
			self.assertFalse(Ask(sPartialTag, bMatchEntireTagOnly = True))
		
	def test_ContainsBothFlankingSequences(self):
		pass		
		

	def test_ContainsBothFlankingSequences_Complement(self):
		pass		





if __name__ == '__main__':
	unittest.main()
