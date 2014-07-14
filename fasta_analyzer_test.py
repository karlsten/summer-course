from nose.tools import *
from fasta_analyzer import fasta_analyzer
import os


# Name and seq used for testing the class.
name = ">Sekvensnamn\n"
seq = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA"


fasta_filename = 'setup_fasta.fa'
fasta_file = open(fasta_filename, 'a+')
coverage_filename = 'setup_cov.txt'
coverage_file = open(coverage_filename, 'a+')

# Create temporary fasta and coverage files that will be used in testing.
def setup():
    fasta_file.write('>contig1\n')
    fasta_file.write('GATTACAGATTACAGATTACAGATTACA\nGATTACAGATTACAGATTACAGATTACA\nGATTACAGATTACAGATTACAGATTACA\n')
    fasta_file.write('>contig2\n')
    fasta_file.write('GGATTACAGGATTACAGGATTACAGGATTACA\nGGATTACAGGATTACAGGATTACAGGATTACA\nGGATTACAGGATTACAGGATTACAGGATTACA\n')
    fasta_file.write('>contig3\n')
    fasta_file.write('GGATTACCAGGATTACCAGGATTACCAGGATTACCA\nGGATTACCAGGATTACCAGGATTACCAGGATTACCA\nGGATTACCAGGATTACCAGGATTACCAGGATTACCA\n')
    coverage_file.write('contig1\t1000')
    coverage_file.write('contig2\t1001')
    coverage_file.write('contig3\t1002')


# Test if the fasta infile is read and saved in dictionary as expected. 
def test_read_file():
    fasta_file.seek(0)
    read_file = fasta_analyzer.read_file(fasta_file)
    assert_equal(read_file['contig1'].header(), 'contig1')
    assert_equal(read_file['contig1'].sequence(), 'GATTACAGATTACAGATTACAGATTACA\nGATTACAGATTACAGATTACAGATTACA\nGATTACAGATTACAGATTACAGATTACA\n')


# Test the Fasta class.
def test_class():
    fasta = fasta_analyzer.Fasta(name, seq)
    assert_equal(fasta.name, name)
    assert_equal(fasta.seq, seq)
    assert_equal(fasta.header(), "Sekvensnamn")
    assert_equal(fasta.sequence(), fasta.seq)
    assert_equal(fasta.length(), 77)
    assert_equal(fasta.gccount(), float(2)/7 * 100)
    fasta.cov = "1234"
    assert_equal(fasta.coverage(), "1234")


# Test if the coverage file is read and saved in dictionary as expected. --NOT DONE--
#def test_read_covfile():
#    coverage_file.seek(0)
#    read_covfile = fasta_analyzer.read_covfile(coverage_file)
#    assert_equal(read_covfile['contig1'].cov, '1000')
    

# Close files so that they can be removed.
def cleanup():
    fasta_file.close()
    coverage_file.close()


# Remove both files after all tests are performed.
os.remove('setup_fasta.fa')
os.remove('setup_cov.txt')
