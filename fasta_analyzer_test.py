from nose.tools import *
from fasta_analyzer import fasta_analyzer
import os

name = ">Sekvensnamn\n"
seq = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA"


filename = 'setup_fasta.fa'
fasta_file = open(filename, 'a')


#def setup():
fasta_file.write('>contig1\n')
fasta_file.write('GATTACAGATTACAGATTACAGATTACA\nGATTACAGATTACAGATTACAGATTACA\nGATTACAGATTACAGATTACAGATTACA\n')
fasta_file.write('>contig2\n')
fasta_file.write('GGATTACAGGATTACAGGATTACAGGATTACA\nGGATTACAGGATTACAGGATTACAGGATTACA\nGGATTACAGGATTACAGGATTACAGGATTACA\n')
fasta_file.write('>contig3\n')
fasta_file.write('GGATTACCAGGATTACCAGGATTACCAGGATTACCA\nGGATTACCAGGATTACCAGGATTACCAGGATTACCA\nGGATTACCAGGATTACCAGGATTACCAGGATTACCA\n')



#def test_read_file():
#read_file = fasta_analyzer.read_file('setup_fasta.fa')
#assert_equal(



#def test_class():
fasta = fasta_analyzer.Fasta(name, seq)
assert_equal(fasta.name, name)
assert_equal(fasta.seq, seq)
assert_equal(fasta.header(), "Sekvensnamn")
assert_equal(fasta.sequence(), fasta.seq)
assert_equal(fasta.length(), 77)
assert_equal(fasta.gccount(), float(2)/7 * 100)
fasta.cov = "1234"
assert_equal(fasta.coverage(), "1234")


#def cleanup():
fasta_file.close()
os.remove('setup_fasta.fa')


