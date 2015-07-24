__author__ = 'Shashank'
#
# -- basicCodeBlock.py
#
from pymol import cmd, stored

def yourFunction( arg1, arg2 ):
    '''
DESCRIPTION

    Brief description what this function does goes here
    '''
    #
    # Your code goes here
    #
    print "Hello, PyMOLers"
    print "You passed in %s and %s" % (arg1, arg2)
    print "I will return them to you in a list.  Here you go."
    return (arg1, arg2)

cmd.extend( "yourFunction", yourFunction("Hello","World") );
cmd.load("../PDB/1JS1.pdb")
