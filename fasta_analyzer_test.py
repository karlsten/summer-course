from nose.tools import *
import fasta_analyzer
import os
from mock import patch





fasta_filename = 'setup_fasta.fa'
fasta_file = open(fasta_filename, 'a+')
coverage_filename = 'setup_cov.txt'
coverage_file = open(coverage_filename, 'a+')

test_dict = {
'contig1' : fasta_analyzer.Fasta('>contig1', 'GATTACA'), 
'contig2' : fasta_analyzer.Fasta('>contig2', 'GATTACAGATTACAGC'), 
'contig3' : fasta_analyzer.Fasta('>contig3', 'GATTACAGATTACAGATTACAGCGC')
}

dict_longcontigs = {
'contig1' : fasta_analyzer.Fasta('>contig1', 'GATTACA'), 
'contig2' : fasta_analyzer.Fasta('>contig2', "GATTACAGC"*2000), 
'contig3' : fasta_analyzer.Fasta('>contig3', "GATTACAGCGCG"*10000)
}

dict_longcontigs['contig1'].cov = 1000
dict_longcontigs['contig2'].cov = 1001
dict_longcontigs['contig3'].cov = 1002





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
    assert_equal(fasta.gccount(), round(float(2)/7 * 100, 1))
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
    read_covfile = fasta_analyzer.read_covfile(coverage_file, test_dict)
    assert_equal(read_covfile['contig1'].cov, int(1000))
    assert_equal(read_covfile['contig2'].cov, int(1001))
    assert_equal(read_covfile['contig3'].cov, int(1002))





# Lines containing a %% must be commented out to enable matplot.pyplot.show.
#def test_lengcplot(): # This line must be uncommented to enable 
# matplot.pyplot.show.
@patch("matplotlib.pyplot.show") # %%
def test_lengcplot(mock_pyplot_show): # %%
    mock_pyplot_show.return_value = None # %%
    xlist, ylist, namelist = fasta_analyzer.lengcplot(test_dict)
    assert_equal(xlist, [25,16,7])
    assert_equal(ylist, [40.0, 37.5, 28.6])
    assert_equal(namelist, ['contig3', 'contig2', 'contig1'])





# Lines containing a %% must be commented out to enable matplot.pyplot.show.
#def test_covgcplot(): # This line must be uncommented to enable 
# matplot.pyplot.show.
@patch("matplotlib.pyplot.show") # %%
def test_covgcplot(mock_pyplot_show): # %%
    mock_pyplot_show.return_value = None # %%
    dict_longcontigs['contig1'].cov = 1000
    dict_longcontigs['contig2'].cov = 1001
    dict_longcontigs['contig3'].cov = 1002
    xlist, xlist1, xlist2, xlist3, ylist, ylist1, ylist2, ylist3, \
    namelist = fasta_analyzer.covgcplot(dict_longcontigs)
    assert_equal(xlist, [1002, 1001, 1000])
    assert_equal(xlist1, [1000])
    assert_equal(xlist2, [1001])
    assert_equal(xlist3, [1002])
    assert_equal(ylist, [58.3, 44.4, 28.6])
    assert_equal(ylist1, [28.6])
    assert_equal(ylist2, [44.4])
    assert_equal(ylist3, [58.3])
    assert_equal(namelist, ['contig3', 'contig2', 'contig1'])





# Lines containing a %% must be commented out to enable matplot.pyplot.show.
#def test_covlenplot(): # This line must be uncommented to enable 
# matplot.pyplot.show.
@patch("matplotlib.pyplot.show") # %%
def test_covlenplot(mock_pyplot_show): # %%
    mock_pyplot_show.return_value = None # %%
    xlist, xlist1, xlist2, xlist3, ylist, ylist1, ylist2, ylist3, \
    namelist = fasta_analyzer.covlenplot(dict_longcontigs)
    assert_equal(xlist, [1002, 1001, 1000])
    assert_equal(xlist1, [1000])
    assert_equal(xlist2, [1001])
    assert_equal(xlist3, [1002])
    assert_equal(ylist, [120000, 18000, 7])
    assert_equal(ylist1, [7])
    assert_equal(ylist2, [18000])
    assert_equal(ylist3, [120000])
    assert_equal(namelist, ['contig3', 'contig2', 'contig1'])





# Lines containing a %% must be commented out to enable matplot.pyplot.show.
#def test_covhistogram(): # This line must be uncommented to enable 
# matplot.pyplot.show.
@patch("matplotlib.pyplot.show") # %%
def test_covhistogram(mock_pyplot_show): # %%
    mock_pyplot_show.return_value = None # %%
    histlist = fasta_analyzer.covhistogram(dict_longcontigs)
    assert_equal(histlist, [1002, 1001, 1000])





# Lines containing a %% must be commented out to enable matplot.pyplot.show.
#def test_lenhistogram(): # This line must be uncommented to enable 
# matplot.pyplot.show.
@patch("matplotlib.pyplot.show") # %%
def test_lenhistogram(mock_pyplot_show): # %%
    mock_pyplot_show.return_value = None # %%
    histlist = fasta_analyzer.lenhistogram(dict_longcontigs)
    assert_equal(histlist, [120000, 18000, 7])





# Close files so that they can be removed.
def cleanup():
    fasta_file.close()
    coverage_file.close()





# Remove both files after all tests are performed.
os.remove('setup_fasta.fa')
os.remove('setup_cov.txt')
