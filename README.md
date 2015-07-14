SlipKnotFind
============

- KnotFind
------------
The code accepts a **.pdb** file as a input parameter and searches the entire chain of *Alpha Carbons* for knots by trimming residues one at a time and ends up with either no knots or the smallest knot core that is found.

-  SlipKnot Detector
--------------------
	A slipknot is formed when the remainder of the knotted protien double backs on itself, making the traditional knotfind algorithms not detect it as a knot. This algorithm here searches for the k1 atom with which the knot can be removed by a simple pull. If the protein has a slipknot, it returns a **.pdb** file with core of the slipknot.

-  Python: PyMol Visualization 
------------------------------
A simple python script using the PyMol API to connect all the *Alpha Carbons* using the *cmd.bond()* method and then coloring it using *cmd.spectrum()*.