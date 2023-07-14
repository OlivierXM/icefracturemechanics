'''
    Viscoelastic Compressive Damage Tensor

'''

import sympy as sym
sym.init_printing()
# Declare Variables
## - d: Damage (0 - 1) // (damaged - undamaged)
## - V: Total anelastic compliance coefficient
## - E: Elastic Modulus
d, V, E = sym.symbols('d V E')

# Elastic Compliance Matrix
A = sym.Matrix([[(2*d+1)/(E*d*(d+2)), (d-1)/(E*d*(d+2)), 0/E], [(d-1)/(E*d*(d+2)), (2*d+1)/(E*d*(d+2)), 0/E],[ 0/E, 0/E, 1/(E*d)]])
# Anelastic Compliance Matrix
B = sym.Matrix([[(2*V)/(3*d), (-V)/(3*d), 0*V], [(-V)/(3*d), (2*V)/(3*d), 0*V],[ 0*V, 0*V, V/d]])

# Combine and invert compliance matrix
D = (A+B).inv()

# Simplify the stiffness matrix
D_s = sym.simplify(D)

# Highlight various test cases
## No Damage, no anelastic mechanisms
## E * I
print(D_s.subs([(d, 1), (V, 0.0)]))
## Fully Damaged, no anelastic mechanisms
## E/3 * [1 1 0
##        1 1 0
##        0 0 0]
print(D_s.subs([(V, 0.0), (d, 1e-9)]))

# Write to File
# out_file = open("latex.txt","w")
# out_file.write(sym.latex(D_s))
# out_file.close()

# Explore incremental previous stress coefficient
# d1, d2 = sym.symbols('d1 d2')
# print(sym.simplify(D_s.subs([(d, 1e-9),(V, 0)])))
# print( sym.simplify((D_s.subs([(d, d1), (V, 0.0)])*(D_s.subs([(d, d2), (V, 0.0)]).inv()))) )

