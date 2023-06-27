# This file was automatically created by FeynRules 1.7.221
# Mathematica version: 9.0 for Linux x86 (64-bit) (November 20, 2012)
# Date: Tue 28 Apr 2015 17:36:04


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
             color = [ 'd(-1,1,4)*f(-1,2,3)', 'd(-1,2,4)*f(-1,1,3)', 'd(-1,3,4)*f(-1,1,2)' ],
             lorentz = [ L.VVVS1, L.VVVS2, L.VVVS3 ],
             couplings = {(0,0):C.GC_4,(1,1):C.GC_4,(2,2):C.GC_4})

V_3 = Vertex(name = 'V_3',
             particles = [ P.G, P.G, P.G, P.G, P.S8 ],
             color = [ 'd(-1,-2,5)*f(-2,1,2)*f(-1,3,4)', 'd(-1,-2,5)*f(-2,1,3)*f(-1,2,4)', 'd(-1,-2,5)*f(-2,1,4)*f(-1,2,3)', 'd(-1,-2,5)*f(-2,2,3)*f(-1,1,4)', 'd(-1,-2,5)*f(-2,2,4)*f(-1,1,3)', 'd(-1,-2,5)*f(-2,3,4)*f(-1,1,2)' ],
             lorentz = [ L.VVVVS1, L.VVVVS2, L.VVVVS3 ],
             couplings = {(0,0):C.GC_5,(1,1):C.GC_5,(4,1):C.GC_5,(5,0):C.GC_5,(2,2):C.GC_5,(3,2):C.GC_5})

V_4 = Vertex(name = 'V_4',
             particles = [ P.G, P.S8, P.S8 ],
             color = [ 'f(1,2,3)' ],
             lorentz = [ L.VSS1 ],
             couplings = {(0,0):C.GC_1})

V_5 = Vertex(name = 'V_5',
             particles = [ P.G, P.G, P.S8, P.S8 ],
             color = [ 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
             lorentz = [ L.VVSS1 ],
             couplings = {(1,0):C.GC_2,(0,0):C.GC_2})

