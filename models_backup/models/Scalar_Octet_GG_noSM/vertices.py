# This file was automatically created by FeynRules $Revision: 436 $
# Mathematica version: 7.0 for Linux x86 (32-bit) (February 18, 2009)
# Date: Thu 19 Apr 2012 11:28:04


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.G, P.G, P.S8 ],
             color = [ 'd(1,2,3)' ],
             lorentz = [ L.VVS1 ],
             couplings = {(0,0):C.GC_3})

V_2 = Vertex(name = 'V_2',
             particles = [ P.G, P.G, P.G, P.S8 ],
             color = [ 'd(1,4,-1)*f(2,3,-1)', 'd(2,4,-1)*f(1,3,-1)', 'd(3,4,-1)*f(1,2,-1)' ],
             lorentz = [ L.VVVS1, L.VVVS2, L.VVVS3 ],
             couplings = {(0,0):C.GC_4,(1,1):C.GC_4,(2,2):C.GC_4})

V_3 = Vertex(name = 'V_3',
             particles = [ P.G, P.G, P.G, P.G, P.S8 ],
             color = [ 'd(5,-1,-2)*f(1,2,-1)*f(3,4,-2)', 'd(5,-1,-2)*f(1,2,-2)*f(3,4,-1)', 'd(5,-1,-2)*f(1,3,-1)*f(2,4,-2)', 'd(5,-1,-2)*f(1,3,-2)*f(2,4,-1)', 'd(5,-1,-2)*f(1,4,-1)*f(2,3,-2)', 'd(5,-1,-2)*f(1,4,-2)*f(2,3,-1)' ],
             lorentz = [ L.VVVVS1, L.VVVVS2, L.VVVVS3 ],
             couplings = {(3,1):C.GC_5,(2,1):C.GC_5,(1,0):C.GC_5,(0,0):C.GC_5,(5,2):C.GC_5,(4,2):C.GC_5})

V_4 = Vertex(name = 'V_4',
             particles = [ P.G, P.S8, P.S8 ],
             color = [ 'f(2,1,3)' ],
             lorentz = [ L.VSS1 ],
             couplings = {(0,0):C.GC_1})

V_5 = Vertex(name = 'V_5',
             particles = [ P.G, P.G, P.S8, P.S8 ],
             color = [ 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
             lorentz = [ L.VVSS1 ],
             couplings = {(1,0):C.GC_2,(0,0):C.GC_2})

