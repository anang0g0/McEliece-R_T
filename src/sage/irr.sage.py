
# This file was *autogenerated* from the file irr.sage
from sage.all_cmdline import *   # import sage library

_sage_const_128 = Integer(128); _sage_const_1 = Integer(1); _sage_const_16 = Integer(16)
BIN = GF(_sage_const_16 )['X']; (X,) = BIN._first_ngens(1)
while _sage_const_1 :
      poly = BIN.random_element(_sage_const_128 );
      if poly.is_irreducible():
      	 break;
print(poly)
 

