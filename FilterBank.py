import re
from collections import namedtuple

tTissueTag = namedtuple('tTissueTag', ['origin', 'label', 'tag'])
lTissueTagList = [tTissueTag('Hsa-Pat6-Aut', 'huh?', ''),
	tTissueTag('Hsa-Pat6-Aut', 'huh?', 'ATTACC'),
	tTissueTag('Hsa-Pat6-Aut', 'huh?', 'TGTGTT'),
	tTissueTag('Hsa-Pat6-Aut', 'huh?', 'CATATT'),
	tTissueTag('Hsa-Pat6-Aut', 'huh?', 'ATACTT'),
	tTissueTag('Hsa-Pat6-Aut', 'huh?', 'GTTACA'),
	tTissueTag('Hsa-Pat6-Aut', 'huh?', 'CTGTAA'),
	tTissueTag('Hsa-Pat6-Aut', 'huh?', 'GAAGAT'),
	tTissueTag('Hsa-Pat6-Aut', 'huh?', 'CATACA'),
	tTissueTag('Hsa-Pat6-Aut', 'huh?', 'CTGATA'), 
]

# compile all regular expressions to be used.
# We're keeping them here, in the global scope, so that
# they only need to be compiled once.
dRegexDict = {}

def ContainsForwardPrimer(sCompleteSequence):
	pass

def ContainsForwardPrimer_Complement(sCompleteSequence):
	pass

# sequence begins with a snippet of tissue tag.
def BeginsWithTissueTag(sCompleteSequence):
	for sTag in lTissueTagList:
		if sCompleteSequence[:len(sTag.tag)] == sTag.tag:
			return True
	return False

def ExtractInsSequence(sCompleteSequence):
	pass

def ExtractInsSequence_Complement(sCompleteSequence):
	pass

# sequence is in the forward direction.

# sequence is in the reverse direction.

# contains a forward primer.

# contains a forward primer, but reversed.

# contains a beginning flanking sequence OR its complement.


#class SequenceSubstringCheck():
#"""A basic template class to be subclassed, below.
#
#"""
#
#	def __init__(self):
#		self.oCompiledRegex = None
#
#	def CheckString(self, sDNAString):
#		pass
#
#	def UnitTest(self):
#		pass
#
#class ContainsForwardPrimer(SequenceSubstringCheck):
#
#	def __init__(self):
#		self.oCompiledRegex = re.match('')
#
#	def CheckString(self, sDNAString):
#		a

