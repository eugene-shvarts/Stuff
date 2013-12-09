import fnmatch
import os

"""
HW 0, Problem 3

This code includes two functions, scissors and glue. scissors disassembles a string into its characters, and stores
each character in a separate text file. glue takes a directory of such files, and reconstructs the original string.

We are given a string, sentence, and perform glue(scissors(sentence,'letters')). The expectation is that we are returned
the string, sentence.
"""


def scissors(string, dirname):
    """
    Creates a directory dirname, and populates it with files 'out_*.txt', the nth of which contains the nth
    character of the provided string.

    @param string: a string which is to be disassembled.
    @param dirname: the directory containing text files, in turn containing characters from the disassembled string.
    @return: returns dirname, for convenience of later applying the glue function.
    """

    l = list(string)
    n = len(l)

    current = os.open(os.getcwd(), os.O_DIRECTORY)    # keep track of current directory
    os.mkdir(dirname)                               # put letters in their own directory
    os.chdir(dirname)

    for i in range(1, n+1):
        if i < 10:                      # for consistent filenames
            fnum = repr(0)+repr(i)
        else:
            fnum = repr(i)
        fname = 'out_'+fnum+'.txt'
        f = open(fname, 'w')
        f.write(l[i-1])          # write a character to its file
        f.close()

    os.fchdir(current)
    return dirname


def glue(ldir):
    """
    Finds all files in ldir matching 'out_*.txt', i.e., those created by the scissors function. Grabs their content
    (each contains one character) into a list, then joins the list to return the string. By virtue of the linear
    naming system used (i.e. the nth character is sent to the nth file by scissors), the string is reconstructed in
    the right order.

    @param ldir: the directory containing text files, in turn containing characters from the disassembled string.
    @return: returns the reassembled, or un-disassembled, string.
    """

    current = os.open(os.getcwd(), os.O_DIRECTORY)    # keep track of current directory
    os.chdir(ldir)
    filelist = [file for file in os.listdir() if fnmatch.fnmatch(file, 'out_*.txt')]

    lstring = []
    for fname in filelist:
        f = open(fname)
        lstring += f.read()     # assemble a list of characters
        f.close()

    os.fchdir(current)
    return ''.join(lstring)     # the reconstructed string

sentence = 'Hello, my name is Bob. I am a statistician. I like statistics very much.'

print(glue(scissors(sentence,'letters')))