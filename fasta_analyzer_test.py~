from nose.tools import *
import fasta_analyzer
import os





fasta_filename = 'setup_fasta.fa'
fasta_file = open(fasta_filename, 'a+')
coverage_filename = 'setup_cov.txt'
coverage_file = open(coverage_filename, 'a+')





# Create temporary fasta and coverage files that will be used in testing.
def setup():
    fasta_file.write('>contig1\n')
    fasta_file.write('GATTACAGATTACAGATTACAGATTACA\n'\
                     'GATTACAGATTACAGATTACAGATTACA\n'\
                     'GATTACAGATTACAGATTACAGATTACAN\n')
    fasta_file.write('>contig2\n')
    fasta_file.write('GGATTACAGGATTACAGGATTACAGGATTACA\n'\
                     'GGATTACAGGATTACAGGATTACAGGATTACA\n'\
                     'GGATTACAGGATTACAGGATTACAGGATTACANN\n')
    fasta_file.write('>contig3\n')
    fasta_file.write('GGATTACCAGGATTACCAGGATTACCAGGATTACCA\n'\
                     'GGATTACCAGGATTACCAGGATTACCAGGATTACCA\n'\
                     'GGATTACCAGGATTACCAGGATTACCAGGATTACCANNN\n')
    coverage_file.write('contig1\t1000\n')
    coverage_file.write('contig2\t1001\n')
    coverage_file.write('contig3\t1002\n')





# Test the Fasta class.
def test_class():
    name = ">Sequence_name\n"
    seq = 'gattacagattacagattacagattacagattacagattaca'\
          'gattacagattacagattacagattacagattacannn'
    fasta = fasta_analyzer.Fasta(name, seq)
    assert_equal(fasta.name, name[1:].rstrip())
    assert_equal(fasta.seq, seq)
    assert_equal(fasta.header(), "Sequence_name")
    assert_equal(fasta.sequence(), fasta.seq)
    assert_equal(fasta.length(), 80)
    assert_equal(fasta.gccount(), round(float(2)/7 * 100, 3))
    fasta.cov = "1234"
    assert_equal(fasta.coverage(), "1234")
    assert_equal(fasta.ncontent(), 3)





# Test if the fasta infile is read and saved in dictionary as expected. 
def test_read_file():
    fasta_file.seek(0)
    read_file = fasta_analyzer.read_file(fasta_file)
    assert_equal(read_file['contig1'].header(), 'contig1')
    assert_equal(read_file['contig1'].sequence(), 'gattacagattaca'\
                                                  'gattacagattaca\n'\
                                                  'gattacagattaca'\
                                                  'gattacagattaca\n'\
                                                  'gattacagattaca'\
                                                  'gattacagattacan\n')
    assert_equal(read_file['contig2'].header(), 'contig2')
    assert_equal(read_file['contig2'].sequence(), 'ggattacaggattaca'\
                                                  'ggattacaggattaca\n'\
                                                  'ggattacaggattaca'\
                                                  'ggattacaggattaca\n'\
                                                  'ggattacaggattaca'\
                                                  'ggattacaggattacann\n')
    assert_equal(read_file['contig3'].header(), 'contig3')
    assert_equal(read_file['contig3'].sequence(), 'ggattaccaggattacca'\
                                                  'ggattaccaggattacca\n'\
                                                  'ggattaccaggattacca'\
                                                  'ggattaccaggattacca\n'\
                                                  'ggattaccaggattacca'\
                                                  'ggattaccaggattaccannn\n')





# Test if the coverage file is read and saved in dictionary as expected.
def test_read_covfile():
    coverage_file.seek(0)
    read_covfile = fasta_analyzer.read_covfile(coverage_file)
    assert_equal(read_covfile['contig1'].cov, int(1000))
    assert_equal(read_covfile['contig2'].cov, int(1001))
    assert_equal(read_covfile['contig3'].cov, int(1002))





# Close files so that they can be removed.
def cleanup():
    fasta_file.close()
    coverage_file.close()





# Remove both files after all tests are performed.
os.remove('setup_fasta.fa')
os.remove('setup_cov.txt')
