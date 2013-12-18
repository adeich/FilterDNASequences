import SequenceFileZipper as SFZ

with open('300LinesTest.fna', 'r') as f, open('300LinesTest.qual', 'r') as q:
	for seq_idF, seq_idQ, lFContents, lQContents in SFZ.SeqNamedTupleGenerator(f, q):
		print 'fasta id: ' + seq_idF
		print 'quali id: ' + seq_idQ
