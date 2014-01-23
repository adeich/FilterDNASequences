class Error(Exception):
	pass


class FastaAndQualiDontMatch(Error):
	pass

class SequenceIsStupid(Error):
	pass

class SequenceCheckFailsSelfTest(Error):
	pass

class SequenceCheckUnitTestFail(Error):
	pass
