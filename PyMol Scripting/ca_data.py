# Import PyMOL's stored module.  This will allow us with a 
# way to pull out the PyMOL data and modify it in our script.
# See below.
from pymol import cmd, stored
 
def data( userSelection ):
    # this array will be used to hold the coordinates.  It
    # has access to PyMOL objects and, we have access to it.
    stored.alphaCarbons = []
 
    # let's just get the alpha carbons, so make the
    # selection just for them
    userSelection = userSelection + " and n. CA"
 
    # iterate over state 1, or the userSelection -- this just means
    # for each item in the selection do what the next parameter says.
    # And, that is to append the (x,y,z) coordinates to the stored.alphaCarbon
    # array.
    cmd.iterate_state(1, selector.process(userSelection), "stored.alphaCarbons.append([x,y,z])")
    return stored.alphaCarbons
#userSelection = []
#cmd.extend( "data", data(userSelection) );
 
    # stored.alphaCarbons now has the data you want.