import re
from collections import namedtuple
import ConstantsAndStructures as CS


# This regex looks for the sequence ID line which should be preceding every set of lines 
# that represent a single sequence. Here's an example such ID line:
# >000017_0846_3280 length=72 uaccno=FHVE77H01CCMLE 
oSeqIDSimpleRegex = re.compile('^>.+')
oSeqIDComplexRegex = re.compile('^>([\d_]+)(\slength=(\d+))?(\suaccno=(\w+))?')

def GetNextSeq_Generator(oSeqFile, iStartPosition=0):
	"""Yields next DNA sequence from file.

	Works on both fasta and quali files. Yields all content between two lines
	which begin with '>'.
	"""

	if iStartPosition:
		oSeqFile.seek(iStartPosition)
	
	# first line read must be a sequence ID.
	sCurrentIDLine = oSeqFile.readline()
	if not oSeqIDSimpleRegex.match(sCurrentIDLine):
		raise BaseException('first line should be a sequence ID')

	# get ensuing lines until next sequence ID line.
	sNextLine = oSeqFile.readline()
	lSeqContentLines = []
	while sNextLine:
		if oSeqIDSimpleRegex.match(sNextLine):
			yield sCurrentIDLine, lSeqContentLines
			sCurrentIDLine = sNextLine
			# empty the list of sequence content strings.
			del(lSeqContentLines[:])
		else:
			lSeqContentLines.append(sNextLine)
	
		# read next line from file.
		sNextLine = oSeqFile.readline()

		
	

def SeqNamedTupleGenerator(oFastaFile, oQualiFile):
	"""Yields a tuple containing ID line, DNA sequence, and quali sequence.

	Uses GetNextSeq_Generator, from above.
	Takes in a fasta file object and quali file object.
	"""

	oFastaFileIterator = GetNextSeq_Generator(oFastaFile)
	oQualiFileIterator = GetNextSeq_Generator(oQualiFile)

	for sFastaSeqIDLine, lFastaContentLines in oFastaFileIterator:
		sQualiSeqIDLine, lQualiContentLines = oQualiFileIterator.next()

		sFastaSeqID = oSeqIDComplexRegex.match(sFastaSeqIDLine).group(1)
		sQualiSeqID = oSeqIDComplexRegex.match(sQualiSeqIDLine).group(1)

	#	sInputSeq = 'hello'
	#	sInputID = sFastaSeqID
		sEntireIDInputLine = sFastaSeqIDLine
		sFastaSequence = str(lFastaContentLines)
		sQualiSequence = str(lQualiContentLines)

		# Should yield input_id, input_id_metadata, fasta_sequence, quali_sequence.
		yield sEntireIDInputLine, sFastaSequence, sQualiSequence
			
