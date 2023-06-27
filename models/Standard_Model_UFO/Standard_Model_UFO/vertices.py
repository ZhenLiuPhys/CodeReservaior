# This file was automatically created by FeynRules 2.3.29
# Mathematica version: 11.0.0 for Mac OS X x86 (64-bit) (July 28, 2016)
# Date: Sun 12 Nov 2017 20:21:30


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.g, P.g, P.g ],
             color = [ 'f(1,2,3)' ],
             lorentz = [ L.VVV2 ],
             couplings = {(0,0):C.GC_6})

V_2 = Vertex(name = 'V_2',
             particles = [ P.g, P.g, P.g, P.g ],
             color = [ 'f(-1,1,2)*f(3,4,-1)', 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
             lorentz = [ L.VVVV6, L.VVVV8, L.VVVV9 ],
             couplings = {(1,1):C.GC_8,(0,0):C.GC_8,(2,2):C.GC_8})

V_3 = Vertex(name = 'V_3',
             particles = [ P.a, P.W__minus__, P.W__plus__ ],
             color = [ '1' ],
             lorentz = [ L.VVV2 ],
             couplings = {(0,0):C.GC_4})

V_4 = Vertex(name = 'V_4',
             particles = [ P.a, P.a, P.W__minus__, P.W__plus__ ],
             color = [ '1' ],
             lorentz = [ L.VVVV7 ],
             couplings = {(0,0):C.GC_5})

V_5 = Vertex(name = 'V_5',
             particles = [ P.W__minus__, P.W__plus__, P.Z ],
             color = [ '1' ],
             lorentz = [ L.VVV2 ],
             couplings = {(0,0):C.GC_16})

V_6 = Vertex(name = 'V_6',
             particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__ ],
             color = [ '1' ],
             lorentz = [ L.VVVV7 ],
             couplings = {(0,0):C.GC_11})

V_7 = Vertex(name = 'V_7',
             particles = [ P.a, P.W__minus__, P.W__plus__, P.Z ],
             color = [ '1' ],
             lorentz = [ L.VVVV10 ],
             couplings = {(0,0):C.GC_17})

V_8 = Vertex(name = 'V_8',
             particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z ],
             color = [ '1' ],
             lorentz = [ L.VVVV7 ],
             couplings = {(0,0):C.GC_12})

V_9 = Vertex(name = 'V_9',
             particles = [ P.H, P.H, P.H, P.H ],
             color = [ '1' ],
             lorentz = [ L.SSSS2 ],
             couplings = {(0,0):C.GC_9})

V_10 = Vertex(name = 'V_10',
              particles = [ P.H, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.SSS2 ],
              couplings = {(0,0):C.GC_22})

V_11 = Vertex(name = 'V_11',
              particles = [ P.W__minus__, P.W__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS2 ],
              couplings = {(0,0):C.GC_10})

V_12 = Vertex(name = 'V_12',
              particles = [ P.W__minus__, P.W__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS2 ],
              couplings = {(0,0):C.GC_23})

V_13 = Vertex(name = 'V_13',
              particles = [ P.Z, P.Z, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS2 ],
              couplings = {(0,0):C.GC_21})

V_14 = Vertex(name = 'V_14',
              particles = [ P.Z, P.Z, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS2 ],
              couplings = {(0,0):C.GC_24})

V_15 = Vertex(name = 'V_15',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[lbar[1]]], P.<>CreateObjectParticleName[PartNameMG[lbar[1]]], P.a ],
              color = [ '1' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_3})

V_16 = Vertex(name = 'V_16',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[lbar[2]]], P.<>CreateObjectParticleName[PartNameMG[lbar[2]]], P.a ],
              color = [ '1' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_3})

V_17 = Vertex(name = 'V_17',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[lbar[3]]], P.<>CreateObjectParticleName[PartNameMG[lbar[3]]], P.a ],
              color = [ '1' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_3})

V_18 = Vertex(name = 'V_18',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[1]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[1]]], P.a ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_2})

V_19 = Vertex(name = 'V_19',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[2]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[2]]], P.a ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_2})

V_20 = Vertex(name = 'V_20',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[3]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[3]]], P.a ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_2})

V_21 = Vertex(name = 'V_21',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[1]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[1]]], P.a ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_1})

V_22 = Vertex(name = 'V_22',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[2]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[2]]], P.a ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_1})

V_23 = Vertex(name = 'V_23',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[3]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[3]]], P.a ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_1})

V_24 = Vertex(name = 'V_24',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[1]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[1]]], P.g ],
              color = [ 'T(3,2,1)' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_7})

V_25 = Vertex(name = 'V_25',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[2]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[2]]], P.g ],
              color = [ 'T(3,2,1)' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_7})

V_26 = Vertex(name = 'V_26',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[3]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[3]]], P.g ],
              color = [ 'T(3,2,1)' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_7})

V_27 = Vertex(name = 'V_27',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[1]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[1]]], P.g ],
              color = [ 'T(3,2,1)' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_7})

V_28 = Vertex(name = 'V_28',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[2]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[2]]], P.g ],
              color = [ 'T(3,2,1)' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_7})

V_29 = Vertex(name = 'V_29',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[3]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[3]]], P.g ],
              color = [ 'T(3,2,1)' ],
              lorentz = [ L.FFV4 ],
              couplings = {(0,0):C.GC_7})

V_30 = Vertex(name = 'V_30',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[1]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[1]]], P.W__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_13})

V_31 = Vertex(name = 'V_31',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[2]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[2]]], P.W__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_13})

V_32 = Vertex(name = 'V_32',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[3]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[3]]], P.W__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_13})

V_33 = Vertex(name = 'V_33',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[1]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[1]]], P.W__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_13})

V_34 = Vertex(name = 'V_34',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[2]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[2]]], P.W__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_13})

V_35 = Vertex(name = 'V_35',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[3]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[3]]], P.W__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_13})

V_36 = Vertex(name = 'V_36',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[lbar[1]]], P.<>CreateObjectParticleName[PartNameMG[vlbar[1]]], P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_13})

V_37 = Vertex(name = 'V_37',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[lbar[2]]], P.<>CreateObjectParticleName[PartNameMG[vlbar[2]]], P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_13})

V_38 = Vertex(name = 'V_38',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[lbar[3]]], P.<>CreateObjectParticleName[PartNameMG[vlbar[3]]], P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_13})

V_39 = Vertex(name = 'V_39',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[vlbar[1]]], P.<>CreateObjectParticleName[PartNameMG[lbar[1]]], P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_13})

V_40 = Vertex(name = 'V_40',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[vlbar[2]]], P.<>CreateObjectParticleName[PartNameMG[lbar[2]]], P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_13})

V_41 = Vertex(name = 'V_41',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[vlbar[3]]], P.<>CreateObjectParticleName[PartNameMG[lbar[3]]], P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_13})

V_42 = Vertex(name = 'V_42',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[1]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[1]]], P.Z ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV5, L.FFV8 ],
              couplings = {(0,0):C.GC_15,(0,1):C.GC_18})

V_43 = Vertex(name = 'V_43',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[2]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[2]]], P.Z ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV5, L.FFV8 ],
              couplings = {(0,0):C.GC_15,(0,1):C.GC_18})

V_44 = Vertex(name = 'V_44',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[3]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[3]]], P.Z ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV5, L.FFV8 ],
              couplings = {(0,0):C.GC_15,(0,1):C.GC_18})

V_45 = Vertex(name = 'V_45',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[1]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[1]]], P.Z ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV5, L.FFV6 ],
              couplings = {(0,0):C.GC_14,(0,1):C.GC_18})

V_46 = Vertex(name = 'V_46',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[2]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[2]]], P.Z ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV5, L.FFV6 ],
              couplings = {(0,0):C.GC_14,(0,1):C.GC_18})

V_47 = Vertex(name = 'V_47',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[3]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[3]]], P.Z ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFV5, L.FFV6 ],
              couplings = {(0,0):C.GC_14,(0,1):C.GC_18})

V_48 = Vertex(name = 'V_48',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[vlbar[1]]], P.<>CreateObjectParticleName[PartNameMG[vlbar[1]]], P.Z ],
              color = [ '1' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_20})

V_49 = Vertex(name = 'V_49',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[vlbar[2]]], P.<>CreateObjectParticleName[PartNameMG[vlbar[2]]], P.Z ],
              color = [ '1' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_20})

V_50 = Vertex(name = 'V_50',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[vlbar[3]]], P.<>CreateObjectParticleName[PartNameMG[vlbar[3]]], P.Z ],
              color = [ '1' ],
              lorentz = [ L.FFV5 ],
              couplings = {(0,0):C.GC_20})

V_51 = Vertex(name = 'V_51',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[lbar[1]]], P.<>CreateObjectParticleName[PartNameMG[lbar[1]]], P.Z ],
              color = [ '1' ],
              lorentz = [ L.FFV5, L.FFV7 ],
              couplings = {(0,0):C.GC_14,(0,1):C.GC_19})

V_52 = Vertex(name = 'V_52',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[lbar[2]]], P.<>CreateObjectParticleName[PartNameMG[lbar[2]]], P.Z ],
              color = [ '1' ],
              lorentz = [ L.FFV5, L.FFV7 ],
              couplings = {(0,0):C.GC_14,(0,1):C.GC_19})

V_53 = Vertex(name = 'V_53',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[lbar[3]]], P.<>CreateObjectParticleName[PartNameMG[lbar[3]]], P.Z ],
              color = [ '1' ],
              lorentz = [ L.FFV5, L.FFV7 ],
              couplings = {(0,0):C.GC_14,(0,1):C.GC_19})

V_54 = Vertex(name = 'V_54',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[dqbar[3]]], P.<>CreateObjectParticleName[PartNameMG[dqbar[3]]], P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS3 ],
              couplings = {(0,0):C.GC_25})

V_55 = Vertex(name = 'V_55',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[lbar[3]]], P.<>CreateObjectParticleName[PartNameMG[lbar[3]]], P.H ],
              color = [ '1' ],
              lorentz = [ L.FFS3 ],
              couplings = {(0,0):C.GC_27})

V_56 = Vertex(name = 'V_56',
              particles = [ P.<>CreateObjectParticleName[PartNameMG[uqbar[3]]], P.<>CreateObjectParticleName[PartNameMG[uqbar[3]]], P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS3 ],
              couplings = {(0,0):C.GC_26})

V_57 = Vertex(name = 'V_57',
              particles = [ P.ghG, P.ghG__tilde__, P.g ],
              color = [ 'f(1,2,3)' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_6})

