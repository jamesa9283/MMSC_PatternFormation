import sympy as sym
x, t, a, L = sym.symbols('x t a L')
u = (x**2 - L**2)*t
def pde(u):
    return sym.diff(u, t) - sym.diff(u, x, x)
f = sym.simplify(pde(u))
print(f)