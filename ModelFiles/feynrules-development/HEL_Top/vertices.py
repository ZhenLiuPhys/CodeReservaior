# This file was automatically created by FeynRules 1.7.221
# Mathematica version: 9.0 for Linux x86 (64-bit) (November 20, 2012)
# Date: Thu 3 Sep 2015 11:43:53


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.G0, P.G0, P.G0, P.G0 ],
             color = [ '1' ],
             lorentz = [ L.SSSS1 ],
             couplings = {(0,0):[ C.GC_170, C.GC_103 ]})

V_2 = Vertex(name = 'V_2',
             particles = [ P.G0, P.G0, P.G__minus__, P.G__plus__ ],
             color = [ '1' ],
             lorentz = [ L.SSSS1 ],
             couplings = {(0,0):[ C.GC_168, C.GC_101 ]})

V_3 = Vertex(name = 'V_3',
             particles = [ P.G__minus__, P.G__minus__, P.G__plus__, P.G__plus__ ],
             color = [ '1' ],
             lorentz = [ L.SSSS1 ],
             couplings = {(0,0):[ C.GC_169, C.GC_102 ]})

V_4 = Vertex(name = 'V_4',
             particles = [ P.G0, P.G0, P.H, P.H ],
             color = [ '1' ],
             lorentz = [ L.SSSS1 ],
             couplings = {(0,0):[ C.GC_168, C.GC_101 ]})

V_5 = Vertex(name = 'V_5',
             particles = [ P.G__minus__, P.G__plus__, P.H, P.H ],
             color = [ '1' ],
             lorentz = [ L.SSSS1 ],
             couplings = {(0,0):[ C.GC_168, C.GC_101 ]})

V_6 = Vertex(name = 'V_6',
             particles = [ P.H, P.H, P.H, P.H ],
             color = [ '1' ],
             lorentz = [ L.SSSS1 ],
             couplings = {(0,0):[ C.GC_170, C.GC_103 ]})

V_7 = Vertex(name = 'V_7',
             particles = [ P.G0, P.G0, P.H ],
             color = [ '1' ],
             lorentz = [ L.SSS1 ],
             couplings = {(0,0):[ C.GC_238, C.GC_114 ]})

V_8 = Vertex(name = 'V_8',
             particles = [ P.G__minus__, P.G__plus__, P.H ],
             color = [ '1' ],
             lorentz = [ L.SSS1 ],
             couplings = {(0,0):[ C.GC_238, C.GC_114 ]})

V_9 = Vertex(name = 'V_9',
             particles = [ P.H, P.H, P.H ],
             color = [ '1' ],
             lorentz = [ L.SSS1 ],
             couplings = {(0,0):[ C.GC_239, C.GC_115 ]})

V_10 = Vertex(name = 'V_10',
              particles = [ P.a, P.a, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS1 ],
              couplings = {(0,0):C.GC_6})

V_11 = Vertex(name = 'V_11',
              particles = [ P.a, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VSS1 ],
              couplings = {(0,0):C.GC_3})

V_12 = Vertex(name = 'V_12',
              particles = [ P.ghA, P.ghWm__tilde__, P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_3})

V_13 = Vertex(name = 'V_13',
              particles = [ P.ghA, P.ghWp__tilde__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_4})

V_14 = Vertex(name = 'V_14',
              particles = [ P.ghWm, P.ghA__tilde__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.UUS1 ],
              couplings = {(0,0):C.GC_276})

V_15 = Vertex(name = 'V_15',
              particles = [ P.ghWm, P.ghA__tilde__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_3})

V_16 = Vertex(name = 'V_16',
              particles = [ P.ghWm, P.ghWm__tilde__, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.UUS1 ],
              couplings = {(0,0):C.GC_271})

V_17 = Vertex(name = 'V_17',
              particles = [ P.ghWm, P.ghWm__tilde__, P.H ],
              color = [ '1' ],
              lorentz = [ L.UUS1 ],
              couplings = {(0,0):C.GC_272})

V_18 = Vertex(name = 'V_18',
              particles = [ P.ghWm, P.ghWm__tilde__, P.a ],
              color = [ '1' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_4})

V_19 = Vertex(name = 'V_19',
              particles = [ P.ghWm, P.ghWm__tilde__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_78})

V_20 = Vertex(name = 'V_20',
              particles = [ P.ghWm, P.ghZ__tilde__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.UUS1 ],
              couplings = {(0,0):C.GC_279})

V_21 = Vertex(name = 'V_21',
              particles = [ P.ghWm, P.ghZ__tilde__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_77})

V_22 = Vertex(name = 'V_22',
              particles = [ P.ghWp, P.ghA__tilde__, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.UUS1 ],
              couplings = {(0,0):C.GC_275})

V_23 = Vertex(name = 'V_23',
              particles = [ P.ghWp, P.ghA__tilde__, P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_4})

V_24 = Vertex(name = 'V_24',
              particles = [ P.ghWp, P.ghWp__tilde__, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.UUS1 ],
              couplings = {(0,0):C.GC_274})

V_25 = Vertex(name = 'V_25',
              particles = [ P.ghWp, P.ghWp__tilde__, P.H ],
              color = [ '1' ],
              lorentz = [ L.UUS1 ],
              couplings = {(0,0):C.GC_272})

V_26 = Vertex(name = 'V_26',
              particles = [ P.ghWp, P.ghWp__tilde__, P.a ],
              color = [ '1' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_3})

V_27 = Vertex(name = 'V_27',
              particles = [ P.ghWp, P.ghWp__tilde__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_77})

V_28 = Vertex(name = 'V_28',
              particles = [ P.ghWp, P.ghZ__tilde__, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.UUS1 ],
              couplings = {(0,0):C.GC_278})

V_29 = Vertex(name = 'V_29',
              particles = [ P.ghWp, P.ghZ__tilde__, P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_78})

V_30 = Vertex(name = 'V_30',
              particles = [ P.ghZ, P.ghWm__tilde__, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.UUS1 ],
              couplings = {(0,0):C.GC_280})

V_31 = Vertex(name = 'V_31',
              particles = [ P.ghZ, P.ghWm__tilde__, P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_77})

V_32 = Vertex(name = 'V_32',
              particles = [ P.ghZ, P.ghWp__tilde__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.UUS1 ],
              couplings = {(0,0):C.GC_277})

V_33 = Vertex(name = 'V_33',
              particles = [ P.ghZ, P.ghWp__tilde__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_78})

V_34 = Vertex(name = 'V_34',
              particles = [ P.ghZ, P.ghZ__tilde__, P.H ],
              color = [ '1' ],
              lorentz = [ L.UUS1 ],
              couplings = {(0,0):C.GC_281})

V_35 = Vertex(name = 'V_35',
              particles = [ P.ghG, P.ghG__tilde__, P.g ],
              color = [ 'f(1,2,3)' ],
              lorentz = [ L.UUV1 ],
              couplings = {(0,0):C.GC_10})

V_36 = Vertex(name = 'V_36',
              particles = [ P.g, P.g, P.g ],
              color = [ 'f(1,2,3)' ],
              lorentz = [ L.VVV1 ],
              couplings = {(0,0):C.GC_10})

V_37 = Vertex(name = 'V_37',
              particles = [ P.g, P.g, P.g, P.g ],
              color = [ 'f(-1,1,2)*f(3,4,-1)', 'f(-1,1,3)*f(2,4,-1)', 'f(-1,1,4)*f(2,3,-1)' ],
              lorentz = [ L.VVVV1, L.VVVV3, L.VVVV4 ],
              couplings = {(1,1):C.GC_12,(0,0):C.GC_12,(2,2):C.GC_12})

V_38 = Vertex(name = 'V_38',
              particles = [ P.u__tilde__, P.d, P.G__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4, L.FFS6 ],
              couplings = {(0,0):C.GC_13,(0,2):C.GC_22,(0,1):C.GC_316,(0,3):C.GC_217})

V_39 = Vertex(name = 'V_39',
              particles = [ P.c__tilde__, P.d, P.G__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4 ],
              couplings = {(0,0):C.GC_14,(0,2):C.GC_23,(0,1):C.GC_355})

V_40 = Vertex(name = 'V_40',
              particles = [ P.t__tilde__, P.d, P.G__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4 ],
              couplings = {(0,0):C.GC_15,(0,2):C.GC_24,(0,1):C.GC_394})

V_41 = Vertex(name = 'V_41',
              particles = [ P.u__tilde__, P.s, P.G__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4 ],
              couplings = {(0,0):C.GC_16,(0,2):C.GC_25,(0,1):C.GC_329})

V_42 = Vertex(name = 'V_42',
              particles = [ P.c__tilde__, P.s, P.G__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4, L.FFS6 ],
              couplings = {(0,0):C.GC_17,(0,2):C.GC_26,(0,1):C.GC_368,(0,3):C.GC_217})

V_43 = Vertex(name = 'V_43',
              particles = [ P.t__tilde__, P.s, P.G__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4 ],
              couplings = {(0,0):C.GC_18,(0,2):C.GC_27,(0,1):C.GC_407})

V_44 = Vertex(name = 'V_44',
              particles = [ P.u__tilde__, P.b, P.G__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4 ],
              couplings = {(0,0):C.GC_19,(0,2):C.GC_28,(0,1):C.GC_342})

V_45 = Vertex(name = 'V_45',
              particles = [ P.c__tilde__, P.b, P.G__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4 ],
              couplings = {(0,0):C.GC_20,(0,2):C.GC_29,(0,1):C.GC_381})

V_46 = Vertex(name = 'V_46',
              particles = [ P.t__tilde__, P.b, P.G__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4, L.FFS6 ],
              couplings = {(0,0):C.GC_21,(0,2):C.GC_30,(0,1):C.GC_420,(0,3):C.GC_217})

V_47 = Vertex(name = 'V_47',
              particles = [ P.d__tilde__, P.d, P.G0 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS2, L.FFS3, L.FFS6 ],
              couplings = {(0,0):C.GC_113,(0,2):C.GC_215,(0,1):C.GC_287})

V_48 = Vertex(name = 'V_48',
              particles = [ P.s__tilde__, P.s, P.G0 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS2, L.FFS3, L.FFS6 ],
              couplings = {(0,0):C.GC_113,(0,2):C.GC_215,(0,1):C.GC_297})

V_49 = Vertex(name = 'V_49',
              particles = [ P.b__tilde__, P.b, P.G0 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS2, L.FFS3, L.FFS6 ],
              couplings = {(0,0):C.GC_113,(0,2):C.GC_215,(0,1):C.GC_283})

V_50 = Vertex(name = 'V_50',
              particles = [ P.d__tilde__, P.d, P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS5 ],
              couplings = {(0,0):C.GC_288})

V_51 = Vertex(name = 'V_51',
              particles = [ P.s__tilde__, P.s, P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS5 ],
              couplings = {(0,0):C.GC_298})

V_52 = Vertex(name = 'V_52',
              particles = [ P.b__tilde__, P.b, P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS5 ],
              couplings = {(0,0):C.GC_284})

V_53 = Vertex(name = 'V_53',
              particles = [ P.ve__tilde__, P.e__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.FFS1 ],
              couplings = {(0,0):C.GC_290})

V_54 = Vertex(name = 'V_54',
              particles = [ P.vm__tilde__, P.mu__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.FFS1 ],
              couplings = {(0,0):C.GC_294})

V_55 = Vertex(name = 'V_55',
              particles = [ P.vt__tilde__, P.ta__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.FFS1 ],
              couplings = {(0,0):C.GC_302})

V_56 = Vertex(name = 'V_56',
              particles = [ P.e__plus__, P.e__minus__, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.FFS3 ],
              couplings = {(0,0):C.GC_291})

V_57 = Vertex(name = 'V_57',
              particles = [ P.mu__plus__, P.mu__minus__, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.FFS3 ],
              couplings = {(0,0):C.GC_295})

V_58 = Vertex(name = 'V_58',
              particles = [ P.ta__plus__, P.ta__minus__, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.FFS3 ],
              couplings = {(0,0):C.GC_303})

V_59 = Vertex(name = 'V_59',
              particles = [ P.e__plus__, P.e__minus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.FFS5 ],
              couplings = {(0,0):C.GC_292})

V_60 = Vertex(name = 'V_60',
              particles = [ P.mu__plus__, P.mu__minus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.FFS5 ],
              couplings = {(0,0):C.GC_296})

V_61 = Vertex(name = 'V_61',
              particles = [ P.ta__plus__, P.ta__minus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.FFS5 ],
              couplings = {(0,0):C.GC_304})

V_62 = Vertex(name = 'V_62',
              particles = [ P.u__tilde__, P.u, P.G0 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS2, L.FFS3, L.FFS6 ],
              couplings = {(0,0):C.GC_112,(0,2):C.GC_216,(0,1):C.GC_306})

V_63 = Vertex(name = 'V_63',
              particles = [ P.c__tilde__, P.c, P.G0 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS2, L.FFS3, L.FFS6 ],
              couplings = {(0,0):C.GC_112,(0,2):C.GC_216,(0,1):C.GC_286})

V_64 = Vertex(name = 'V_64',
              particles = [ P.t__tilde__, P.t, P.G0 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS2, L.FFS3, L.FFS6 ],
              couplings = {(0,0):C.GC_112,(0,2):C.GC_216,(0,1):C.GC_300})

V_65 = Vertex(name = 'V_65',
              particles = [ P.u__tilde__, P.u, P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS5 ],
              couplings = {(0,0):C.GC_305})

V_66 = Vertex(name = 'V_66',
              particles = [ P.c__tilde__, P.c, P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS5 ],
              couplings = {(0,0):C.GC_285})

V_67 = Vertex(name = 'V_67',
              particles = [ P.t__tilde__, P.t, P.H ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS5 ],
              couplings = {(0,0):C.GC_299})

V_68 = Vertex(name = 'V_68',
              particles = [ P.d__tilde__, P.u, P.G__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4, L.FFS6 ],
              couplings = {(0,0):C.GC_31,(0,2):C.GC_40,(0,1):C.GC_218,(0,3):C.GC_217})

V_69 = Vertex(name = 'V_69',
              particles = [ P.s__tilde__, P.u, P.G__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4 ],
              couplings = {(0,0):C.GC_32,(0,2):C.GC_41,(0,1):C.GC_219})

V_70 = Vertex(name = 'V_70',
              particles = [ P.b__tilde__, P.u, P.G__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4 ],
              couplings = {(0,0):C.GC_33,(0,2):C.GC_42,(0,1):C.GC_220})

V_71 = Vertex(name = 'V_71',
              particles = [ P.d__tilde__, P.c, P.G__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4 ],
              couplings = {(0,0):C.GC_34,(0,2):C.GC_43,(0,1):C.GC_221})

V_72 = Vertex(name = 'V_72',
              particles = [ P.s__tilde__, P.c, P.G__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4, L.FFS6 ],
              couplings = {(0,0):C.GC_35,(0,2):C.GC_44,(0,1):C.GC_222,(0,3):C.GC_217})

V_73 = Vertex(name = 'V_73',
              particles = [ P.b__tilde__, P.c, P.G__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4 ],
              couplings = {(0,0):C.GC_36,(0,2):C.GC_45,(0,1):C.GC_223})

V_74 = Vertex(name = 'V_74',
              particles = [ P.d__tilde__, P.t, P.G__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4 ],
              couplings = {(0,0):C.GC_37,(0,2):C.GC_46,(0,1):C.GC_224})

V_75 = Vertex(name = 'V_75',
              particles = [ P.s__tilde__, P.t, P.G__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4 ],
              couplings = {(0,0):C.GC_38,(0,2):C.GC_47,(0,1):C.GC_225})

V_76 = Vertex(name = 'V_76',
              particles = [ P.b__tilde__, P.t, P.G__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.FFS1, L.FFS2, L.FFS4, L.FFS6 ],
              couplings = {(0,0):C.GC_39,(0,2):C.GC_48,(0,1):C.GC_226,(0,3):C.GC_217})

V_77 = Vertex(name = 'V_77',
              particles = [ P.a, P.W__minus__, P.G0, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS1 ],
              couplings = {(0,0):C.GC_80})

V_78 = Vertex(name = 'V_78',
              particles = [ P.a, P.W__minus__, P.G__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS1 ],
              couplings = {(0,0):C.GC_79})

V_79 = Vertex(name = 'V_79',
              particles = [ P.a, P.W__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVS1 ],
              couplings = {(0,0):C.GC_275})

V_80 = Vertex(name = 'V_80',
              particles = [ P.W__minus__, P.G0, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VSS1 ],
              couplings = {(0,0):C.GC_54})

V_81 = Vertex(name = 'V_81',
              particles = [ P.W__minus__, P.G__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VSS1 ],
              couplings = {(0,0):C.GC_52})

V_82 = Vertex(name = 'V_82',
              particles = [ P.a, P.W__minus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVV1 ],
              couplings = {(0,0):C.GC_4})

V_83 = Vertex(name = 'V_83',
              particles = [ P.a, P.W__plus__, P.G0, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS1 ],
              couplings = {(0,0):C.GC_80})

V_84 = Vertex(name = 'V_84',
              particles = [ P.a, P.W__plus__, P.G__minus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS1 ],
              couplings = {(0,0):C.GC_81})

V_85 = Vertex(name = 'V_85',
              particles = [ P.a, P.W__plus__, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.VVS1 ],
              couplings = {(0,0):C.GC_276})

V_86 = Vertex(name = 'V_86',
              particles = [ P.W__plus__, P.G0, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.VSS1 ],
              couplings = {(0,0):C.GC_53})

V_87 = Vertex(name = 'V_87',
              particles = [ P.W__plus__, P.G__minus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VSS1 ],
              couplings = {(0,0):C.GC_52})

V_88 = Vertex(name = 'V_88',
              particles = [ P.W__minus__, P.W__plus__, P.G0, P.G0 ],
              color = [ '1' ],
              lorentz = [ L.VVSS1 ],
              couplings = {(0,0):C.GC_49})

V_89 = Vertex(name = 'V_89',
              particles = [ P.W__minus__, P.W__plus__, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS1 ],
              couplings = {(0,0):C.GC_49})

V_90 = Vertex(name = 'V_90',
              particles = [ P.W__minus__, P.W__plus__, P.H, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVSS1 ],
              couplings = {(0,0):C.GC_49})

V_91 = Vertex(name = 'V_91',
              particles = [ P.W__minus__, P.W__plus__, P.H ],
              color = [ '1' ],
              lorentz = [ L.VVS1 ],
              couplings = {(0,0):C.GC_273})

V_92 = Vertex(name = 'V_92',
              particles = [ P.a, P.a, P.W__minus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVV2 ],
              couplings = {(0,0):C.GC_5})

V_93 = Vertex(name = 'V_93',
              particles = [ P.W__minus__, P.W__plus__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.VVV1 ],
              couplings = {(0,0):C.GC_78})

V_94 = Vertex(name = 'V_94',
              particles = [ P.W__minus__, P.W__minus__, P.W__plus__, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVVV2 ],
              couplings = {(0,0):C.GC_50})

V_95 = Vertex(name = 'V_95',
              particles = [ P.e__plus__, P.ve, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.FFS4 ],
              couplings = {(0,0):C.GC_289})

V_96 = Vertex(name = 'V_96',
              particles = [ P.mu__plus__, P.vm, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.FFS4 ],
              couplings = {(0,0):C.GC_293})

V_97 = Vertex(name = 'V_97',
              particles = [ P.ta__plus__, P.vt, P.G__minus__ ],
              color = [ '1' ],
              lorentz = [ L.FFS4 ],
              couplings = {(0,0):C.GC_301})

V_98 = Vertex(name = 'V_98',
              particles = [ P.a, P.Z, P.G__minus__, P.G__plus__ ],
              color = [ '1' ],
              lorentz = [ L.VVSS1 ],
              couplings = {(0,0):C.GC_92})

V_99 = Vertex(name = 'V_99',
              particles = [ P.Z, P.G0, P.H ],
              color = [ '1' ],
              lorentz = [ L.VSS1 ],
              couplings = {(0,0):C.GC_85})

V_100 = Vertex(name = 'V_100',
               particles = [ P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VSS1 ],
               couplings = {(0,0):C.GC_86})

V_101 = Vertex(name = 'V_101',
               particles = [ P.W__minus__, P.Z, P.G0, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSS1 ],
               couplings = {(0,0):C.GC_8})

V_102 = Vertex(name = 'V_102',
               particles = [ P.W__minus__, P.Z, P.G__plus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSS1 ],
               couplings = {(0,0):C.GC_9})

V_103 = Vertex(name = 'V_103',
               particles = [ P.W__minus__, P.Z, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVS1 ],
               couplings = {(0,0):C.GC_270})

V_104 = Vertex(name = 'V_104',
               particles = [ P.W__plus__, P.Z, P.G0, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSS1 ],
               couplings = {(0,0):C.GC_8})

V_105 = Vertex(name = 'V_105',
               particles = [ P.W__plus__, P.Z, P.G__minus__, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSS1 ],
               couplings = {(0,0):C.GC_7})

V_106 = Vertex(name = 'V_106',
               particles = [ P.W__plus__, P.Z, P.G__minus__ ],
               color = [ '1' ],
               lorentz = [ L.VVS1 ],
               couplings = {(0,0):C.GC_269})

V_107 = Vertex(name = 'V_107',
               particles = [ P.a, P.W__minus__, P.W__plus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVV5 ],
               couplings = {(0,0):C.GC_82})

V_108 = Vertex(name = 'V_108',
               particles = [ P.Z, P.Z, P.G0, P.G0 ],
               color = [ '1' ],
               lorentz = [ L.VVSS1 ],
               couplings = {(0,0):C.GC_94})

V_109 = Vertex(name = 'V_109',
               particles = [ P.Z, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ '1' ],
               lorentz = [ L.VVSS1 ],
               couplings = {(0,0):C.GC_93})

V_110 = Vertex(name = 'V_110',
               particles = [ P.Z, P.Z, P.H, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVSS1 ],
               couplings = {(0,0):C.GC_94})

V_111 = Vertex(name = 'V_111',
               particles = [ P.Z, P.Z, P.H ],
               color = [ '1' ],
               lorentz = [ L.VVS1 ],
               couplings = {(0,0):C.GC_282})

V_112 = Vertex(name = 'V_112',
               particles = [ P.W__minus__, P.W__plus__, P.Z, P.Z ],
               color = [ '1' ],
               lorentz = [ L.VVVV2 ],
               couplings = {(0,0):C.GC_51})

V_113 = Vertex(name = 'V_113',
               particles = [ P.e__plus__, P.ve, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_55})

V_114 = Vertex(name = 'V_114',
               particles = [ P.mu__plus__, P.vm, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_55})

V_115 = Vertex(name = 'V_115',
               particles = [ P.ta__plus__, P.vt, P.W__minus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_55})

V_116 = Vertex(name = 'V_116',
               particles = [ P.ve__tilde__, P.ve, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_87})

V_117 = Vertex(name = 'V_117',
               particles = [ P.vm__tilde__, P.vm, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_87})

V_118 = Vertex(name = 'V_118',
               particles = [ P.vt__tilde__, P.vt, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_87})

V_119 = Vertex(name = 'V_119',
               particles = [ P.e__plus__, P.e__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_3})

V_120 = Vertex(name = 'V_120',
               particles = [ P.mu__plus__, P.mu__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_3})

V_121 = Vertex(name = 'V_121',
               particles = [ P.ta__plus__, P.ta__minus__, P.a ],
               color = [ '1' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_3})

V_122 = Vertex(name = 'V_122',
               particles = [ P.ve__tilde__, P.e__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_55})

V_123 = Vertex(name = 'V_123',
               particles = [ P.vm__tilde__, P.mu__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_55})

V_124 = Vertex(name = 'V_124',
               particles = [ P.vt__tilde__, P.ta__minus__, P.W__plus__ ],
               color = [ '1' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):C.GC_55})

V_125 = Vertex(name = 'V_125',
               particles = [ P.e__plus__, P.e__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV2, L.FFV5 ],
               couplings = {(0,0):C.GC_75,(0,1):C.GC_84})

V_126 = Vertex(name = 'V_126',
               particles = [ P.mu__plus__, P.mu__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV2, L.FFV5 ],
               couplings = {(0,1):C.GC_84,(0,0):C.GC_75})

V_127 = Vertex(name = 'V_127',
               particles = [ P.ta__plus__, P.ta__minus__, P.Z ],
               color = [ '1' ],
               lorentz = [ L.FFV2, L.FFV5 ],
               couplings = {(0,0):C.GC_75,(0,1):C.GC_84})

V_128 = Vertex(name = 'V_128',
               particles = [ P.u__tilde__, P.u, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_96,(0,1):C.GC_122})

V_129 = Vertex(name = 'V_129',
               particles = [ P.c__tilde__, P.c, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_96,(0,1):C.GC_122})

V_130 = Vertex(name = 'V_130',
               particles = [ P.t__tilde__, P.t, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_96,(0,1):C.GC_122})

V_131 = Vertex(name = 'V_131',
               particles = [ P.d__tilde__, P.u, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_128,(0,1):C.GC_126})

V_132 = Vertex(name = 'V_132',
               particles = [ P.s__tilde__, P.u, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_130})

V_133 = Vertex(name = 'V_133',
               particles = [ P.b__tilde__, P.u, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_132})

V_134 = Vertex(name = 'V_134',
               particles = [ P.d__tilde__, P.c, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_134})

V_135 = Vertex(name = 'V_135',
               particles = [ P.s__tilde__, P.c, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_136,(0,1):C.GC_126})

V_136 = Vertex(name = 'V_136',
               particles = [ P.b__tilde__, P.c, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_138})

V_137 = Vertex(name = 'V_137',
               particles = [ P.d__tilde__, P.t, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_140})

V_138 = Vertex(name = 'V_138',
               particles = [ P.s__tilde__, P.t, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_142})

V_139 = Vertex(name = 'V_139',
               particles = [ P.b__tilde__, P.t, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_144,(0,1):C.GC_126})

V_140 = Vertex(name = 'V_140',
               particles = [ P.d__tilde__, P.u, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_127,(0,1):C.GC_124})

V_141 = Vertex(name = 'V_141',
               particles = [ P.s__tilde__, P.u, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_129})

V_142 = Vertex(name = 'V_142',
               particles = [ P.b__tilde__, P.u, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_131})

V_143 = Vertex(name = 'V_143',
               particles = [ P.d__tilde__, P.c, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_133})

V_144 = Vertex(name = 'V_144',
               particles = [ P.s__tilde__, P.c, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_135,(0,1):C.GC_124})

V_145 = Vertex(name = 'V_145',
               particles = [ P.b__tilde__, P.c, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_137})

V_146 = Vertex(name = 'V_146',
               particles = [ P.d__tilde__, P.t, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_139})

V_147 = Vertex(name = 'V_147',
               particles = [ P.s__tilde__, P.t, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_141})

V_148 = Vertex(name = 'V_148',
               particles = [ P.b__tilde__, P.t, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_143,(0,1):C.GC_124})

V_149 = Vertex(name = 'V_149',
               particles = [ P.u__tilde__, P.u, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_95,(0,1):C.GC_123})

V_150 = Vertex(name = 'V_150',
               particles = [ P.c__tilde__, P.c, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_95,(0,1):C.GC_123})

V_151 = Vertex(name = 'V_151',
               particles = [ P.t__tilde__, P.t, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_95,(0,1):C.GC_123})

V_152 = Vertex(name = 'V_152',
               particles = [ P.u__tilde__, P.u, P.a, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_100,(0,1):C.GC_146})

V_153 = Vertex(name = 'V_153',
               particles = [ P.c__tilde__, P.c, P.a, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_100,(0,1):C.GC_146})

V_154 = Vertex(name = 'V_154',
               particles = [ P.t__tilde__, P.t, P.a, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_100,(0,1):C.GC_146})

V_155 = Vertex(name = 'V_155',
               particles = [ P.d__tilde__, P.u, P.a, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_150,(0,1):C.GC_148})

V_156 = Vertex(name = 'V_156',
               particles = [ P.s__tilde__, P.u, P.a, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_152})

V_157 = Vertex(name = 'V_157',
               particles = [ P.b__tilde__, P.u, P.a, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_154})

V_158 = Vertex(name = 'V_158',
               particles = [ P.d__tilde__, P.c, P.a, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_156})

V_159 = Vertex(name = 'V_159',
               particles = [ P.s__tilde__, P.c, P.a, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_158,(0,1):C.GC_148})

V_160 = Vertex(name = 'V_160',
               particles = [ P.b__tilde__, P.c, P.a, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_160})

V_161 = Vertex(name = 'V_161',
               particles = [ P.d__tilde__, P.t, P.a, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_162})

V_162 = Vertex(name = 'V_162',
               particles = [ P.s__tilde__, P.t, P.a, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_164})

V_163 = Vertex(name = 'V_163',
               particles = [ P.b__tilde__, P.t, P.a, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_166,(0,1):C.GC_148})

V_164 = Vertex(name = 'V_164',
               particles = [ P.d__tilde__, P.u, P.a, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_151,(0,1):C.GC_149})

V_165 = Vertex(name = 'V_165',
               particles = [ P.s__tilde__, P.u, P.a, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_153})

V_166 = Vertex(name = 'V_166',
               particles = [ P.b__tilde__, P.u, P.a, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_155})

V_167 = Vertex(name = 'V_167',
               particles = [ P.d__tilde__, P.c, P.a, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_157})

V_168 = Vertex(name = 'V_168',
               particles = [ P.s__tilde__, P.c, P.a, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_159,(0,1):C.GC_149})

V_169 = Vertex(name = 'V_169',
               particles = [ P.b__tilde__, P.c, P.a, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_161})

V_170 = Vertex(name = 'V_170',
               particles = [ P.d__tilde__, P.t, P.a, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_163})

V_171 = Vertex(name = 'V_171',
               particles = [ P.s__tilde__, P.t, P.a, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_165})

V_172 = Vertex(name = 'V_172',
               particles = [ P.b__tilde__, P.t, P.a, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_167,(0,1):C.GC_149})

V_173 = Vertex(name = 'V_173',
               particles = [ P.d__tilde__, P.u, P.a, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_229,(0,1):C.GC_228})

V_174 = Vertex(name = 'V_174',
               particles = [ P.s__tilde__, P.u, P.a, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_230})

V_175 = Vertex(name = 'V_175',
               particles = [ P.b__tilde__, P.u, P.a, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_231})

V_176 = Vertex(name = 'V_176',
               particles = [ P.d__tilde__, P.c, P.a, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_232})

V_177 = Vertex(name = 'V_177',
               particles = [ P.s__tilde__, P.c, P.a, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_233,(0,1):C.GC_228})

V_178 = Vertex(name = 'V_178',
               particles = [ P.b__tilde__, P.c, P.a, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_234})

V_179 = Vertex(name = 'V_179',
               particles = [ P.d__tilde__, P.t, P.a, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_235})

V_180 = Vertex(name = 'V_180',
               particles = [ P.s__tilde__, P.t, P.a, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_236})

V_181 = Vertex(name = 'V_181',
               particles = [ P.b__tilde__, P.t, P.a, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_237,(0,1):C.GC_228})

V_182 = Vertex(name = 'V_182',
               particles = [ P.u__tilde__, P.u, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_2})

V_183 = Vertex(name = 'V_183',
               particles = [ P.c__tilde__, P.c, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_2})

V_184 = Vertex(name = 'V_184',
               particles = [ P.t__tilde__, P.t, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_2})

V_185 = Vertex(name = 'V_185',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_185})

V_186 = Vertex(name = 'V_186',
               particles = [ P.s__tilde__, P.u, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_186})

V_187 = Vertex(name = 'V_187',
               particles = [ P.b__tilde__, P.u, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_187})

V_188 = Vertex(name = 'V_188',
               particles = [ P.d__tilde__, P.c, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_188})

V_189 = Vertex(name = 'V_189',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_189})

V_190 = Vertex(name = 'V_190',
               particles = [ P.b__tilde__, P.c, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_190})

V_191 = Vertex(name = 'V_191',
               particles = [ P.d__tilde__, P.t, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_191})

V_192 = Vertex(name = 'V_192',
               particles = [ P.s__tilde__, P.t, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_192})

V_193 = Vertex(name = 'V_193',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_193})

V_194 = Vertex(name = 'V_194',
               particles = [ P.u__tilde__, P.u, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_175,(0,1):C.GC_178})

V_195 = Vertex(name = 'V_195',
               particles = [ P.c__tilde__, P.c, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_175,(0,1):C.GC_178})

V_196 = Vertex(name = 'V_196',
               particles = [ P.t__tilde__, P.t, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_175,(0,1):C.GC_178})

V_197 = Vertex(name = 'V_197',
               particles = [ P.u__tilde__, P.u, P.W__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_174,(0,1):C.GC_177})

V_198 = Vertex(name = 'V_198',
               particles = [ P.c__tilde__, P.c, P.W__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_174,(0,1):C.GC_177})

V_199 = Vertex(name = 'V_199',
               particles = [ P.t__tilde__, P.t, P.W__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_174,(0,1):C.GC_177})

V_200 = Vertex(name = 'V_200',
               particles = [ P.u__tilde__, P.u, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_242,(0,1):C.GC_244})

V_201 = Vertex(name = 'V_201',
               particles = [ P.c__tilde__, P.c, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_242,(0,1):C.GC_244})

V_202 = Vertex(name = 'V_202',
               particles = [ P.t__tilde__, P.t, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_242,(0,1):C.GC_244})

V_203 = Vertex(name = 'V_203',
               particles = [ P.u__tilde__, P.u, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_175,(0,1):C.GC_178})

V_204 = Vertex(name = 'V_204',
               particles = [ P.c__tilde__, P.c, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_175,(0,1):C.GC_178})

V_205 = Vertex(name = 'V_205',
               particles = [ P.t__tilde__, P.t, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_175,(0,1):C.GC_178})

V_206 = Vertex(name = 'V_206',
               particles = [ P.u__tilde__, P.u, P.W__plus__, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_176,(0,1):C.GC_179})

V_207 = Vertex(name = 'V_207',
               particles = [ P.c__tilde__, P.c, P.W__plus__, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_176,(0,1):C.GC_179})

V_208 = Vertex(name = 'V_208',
               particles = [ P.t__tilde__, P.t, P.W__plus__, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_176,(0,1):C.GC_179})

V_209 = Vertex(name = 'V_209',
               particles = [ P.u__tilde__, P.u, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_243,(0,1):C.GC_245})

V_210 = Vertex(name = 'V_210',
               particles = [ P.c__tilde__, P.c, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_243,(0,1):C.GC_245})

V_211 = Vertex(name = 'V_211',
               particles = [ P.t__tilde__, P.t, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_243,(0,1):C.GC_245})

V_212 = Vertex(name = 'V_212',
               particles = [ P.d__tilde__, P.u, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV3 ],
               couplings = {(0,0):[ C.GC_57, C.GC_66 ],(0,1):C.GC_56})

V_213 = Vertex(name = 'V_213',
               particles = [ P.s__tilde__, P.u, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):[ C.GC_58, C.GC_67 ]})

V_214 = Vertex(name = 'V_214',
               particles = [ P.b__tilde__, P.u, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):[ C.GC_59, C.GC_68 ]})

V_215 = Vertex(name = 'V_215',
               particles = [ P.d__tilde__, P.c, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):[ C.GC_60, C.GC_69 ]})

V_216 = Vertex(name = 'V_216',
               particles = [ P.s__tilde__, P.c, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV3 ],
               couplings = {(0,0):[ C.GC_61, C.GC_70 ],(0,1):C.GC_56})

V_217 = Vertex(name = 'V_217',
               particles = [ P.b__tilde__, P.c, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):[ C.GC_62, C.GC_71 ]})

V_218 = Vertex(name = 'V_218',
               particles = [ P.d__tilde__, P.t, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):[ C.GC_63, C.GC_72 ]})

V_219 = Vertex(name = 'V_219',
               particles = [ P.s__tilde__, P.t, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):[ C.GC_64, C.GC_73 ]})

V_220 = Vertex(name = 'V_220',
               particles = [ P.b__tilde__, P.t, P.W__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV3 ],
               couplings = {(0,0):[ C.GC_65, C.GC_74 ],(0,1):C.GC_56})

V_221 = Vertex(name = 'V_221',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_185,(0,1):C.GC_181})

V_222 = Vertex(name = 'V_222',
               particles = [ P.s__tilde__, P.u, P.W__minus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_186})

V_223 = Vertex(name = 'V_223',
               particles = [ P.b__tilde__, P.u, P.W__minus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_187})

V_224 = Vertex(name = 'V_224',
               particles = [ P.d__tilde__, P.c, P.W__minus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_188})

V_225 = Vertex(name = 'V_225',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_189,(0,1):C.GC_181})

V_226 = Vertex(name = 'V_226',
               particles = [ P.b__tilde__, P.c, P.W__minus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_190})

V_227 = Vertex(name = 'V_227',
               particles = [ P.d__tilde__, P.t, P.W__minus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_191})

V_228 = Vertex(name = 'V_228',
               particles = [ P.s__tilde__, P.t, P.W__minus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_192})

V_229 = Vertex(name = 'V_229',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_193,(0,1):C.GC_181})

V_230 = Vertex(name = 'V_230',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_185,(0,1):C.GC_182})

V_231 = Vertex(name = 'V_231',
               particles = [ P.s__tilde__, P.u, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_186})

V_232 = Vertex(name = 'V_232',
               particles = [ P.b__tilde__, P.u, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_187})

V_233 = Vertex(name = 'V_233',
               particles = [ P.d__tilde__, P.c, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_188})

V_234 = Vertex(name = 'V_234',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_189,(0,1):C.GC_182})

V_235 = Vertex(name = 'V_235',
               particles = [ P.b__tilde__, P.c, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_190})

V_236 = Vertex(name = 'V_236',
               particles = [ P.d__tilde__, P.t, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_191})

V_237 = Vertex(name = 'V_237',
               particles = [ P.s__tilde__, P.t, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_192})

V_238 = Vertex(name = 'V_238',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_193,(0,1):C.GC_182})

V_239 = Vertex(name = 'V_239',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_249,(0,1):C.GC_247})

V_240 = Vertex(name = 'V_240',
               particles = [ P.s__tilde__, P.u, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_250})

V_241 = Vertex(name = 'V_241',
               particles = [ P.b__tilde__, P.u, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_251})

V_242 = Vertex(name = 'V_242',
               particles = [ P.d__tilde__, P.c, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_252})

V_243 = Vertex(name = 'V_243',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_253,(0,1):C.GC_247})

V_244 = Vertex(name = 'V_244',
               particles = [ P.b__tilde__, P.c, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_254})

V_245 = Vertex(name = 'V_245',
               particles = [ P.d__tilde__, P.t, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_255})

V_246 = Vertex(name = 'V_246',
               particles = [ P.s__tilde__, P.t, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_256})

V_247 = Vertex(name = 'V_247',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_257,(0,1):C.GC_247})

V_248 = Vertex(name = 'V_248',
               particles = [ P.u__tilde__, P.u, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_109,(0,1):C.GC_107})

V_249 = Vertex(name = 'V_249',
               particles = [ P.c__tilde__, P.c, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_109,(0,1):C.GC_107})

V_250 = Vertex(name = 'V_250',
               particles = [ P.t__tilde__, P.t, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_109,(0,1):C.GC_107})

V_251 = Vertex(name = 'V_251',
               particles = [ P.u__tilde__, P.u, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV3, L.FFV6 ],
               couplings = {(0,0):[ C.GC_76, C.GC_91 ],(0,2):C.GC_83,(0,1):C.GC_89})

V_252 = Vertex(name = 'V_252',
               particles = [ P.c__tilde__, P.c, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV3, L.FFV6 ],
               couplings = {(0,0):[ C.GC_76, C.GC_91 ],(0,2):C.GC_83,(0,1):C.GC_89})

V_253 = Vertex(name = 'V_253',
               particles = [ P.t__tilde__, P.t, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV3, L.FFV6 ],
               couplings = {(0,0):[ C.GC_76, C.GC_91 ],(0,2):C.GC_83,(0,1):C.GC_89})

V_254 = Vertex(name = 'V_254',
               particles = [ P.u__tilde__, P.u, P.Z, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_111,(0,1):C.GC_106})

V_255 = Vertex(name = 'V_255',
               particles = [ P.c__tilde__, P.c, P.Z, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_111,(0,1):C.GC_106})

V_256 = Vertex(name = 'V_256',
               particles = [ P.t__tilde__, P.t, P.Z, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_111,(0,1):C.GC_106})

V_257 = Vertex(name = 'V_257',
               particles = [ P.u__tilde__, P.u, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_111,(0,1):C.GC_106})

V_258 = Vertex(name = 'V_258',
               particles = [ P.c__tilde__, P.c, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_111,(0,1):C.GC_106})

V_259 = Vertex(name = 'V_259',
               particles = [ P.t__tilde__, P.t, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_111,(0,1):C.GC_106})

V_260 = Vertex(name = 'V_260',
               particles = [ P.u__tilde__, P.u, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_119,(0,1):C.GC_117})

V_261 = Vertex(name = 'V_261',
               particles = [ P.c__tilde__, P.c, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_119,(0,1):C.GC_117})

V_262 = Vertex(name = 'V_262',
               particles = [ P.t__tilde__, P.t, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_119,(0,1):C.GC_117})

V_263 = Vertex(name = 'V_263',
               particles = [ P.d__tilde__, P.u, P.Z, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_198,(0,1):C.GC_195})

V_264 = Vertex(name = 'V_264',
               particles = [ P.s__tilde__, P.u, P.Z, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_200})

V_265 = Vertex(name = 'V_265',
               particles = [ P.b__tilde__, P.u, P.Z, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_202})

V_266 = Vertex(name = 'V_266',
               particles = [ P.d__tilde__, P.c, P.Z, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_204})

V_267 = Vertex(name = 'V_267',
               particles = [ P.s__tilde__, P.c, P.Z, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_206,(0,1):C.GC_195})

V_268 = Vertex(name = 'V_268',
               particles = [ P.b__tilde__, P.c, P.Z, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_208})

V_269 = Vertex(name = 'V_269',
               particles = [ P.d__tilde__, P.t, P.Z, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_210})

V_270 = Vertex(name = 'V_270',
               particles = [ P.s__tilde__, P.t, P.Z, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_212})

V_271 = Vertex(name = 'V_271',
               particles = [ P.b__tilde__, P.t, P.Z, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_214,(0,1):C.GC_195})

V_272 = Vertex(name = 'V_272',
               particles = [ P.d__tilde__, P.u, P.Z, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_197,(0,1):C.GC_196})

V_273 = Vertex(name = 'V_273',
               particles = [ P.s__tilde__, P.u, P.Z, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_199})

V_274 = Vertex(name = 'V_274',
               particles = [ P.b__tilde__, P.u, P.Z, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_201})

V_275 = Vertex(name = 'V_275',
               particles = [ P.d__tilde__, P.c, P.Z, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_203})

V_276 = Vertex(name = 'V_276',
               particles = [ P.s__tilde__, P.c, P.Z, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_205,(0,1):C.GC_196})

V_277 = Vertex(name = 'V_277',
               particles = [ P.b__tilde__, P.c, P.Z, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_207})

V_278 = Vertex(name = 'V_278',
               particles = [ P.d__tilde__, P.t, P.Z, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_209})

V_279 = Vertex(name = 'V_279',
               particles = [ P.s__tilde__, P.t, P.Z, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_211})

V_280 = Vertex(name = 'V_280',
               particles = [ P.b__tilde__, P.t, P.Z, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_213,(0,1):C.GC_196})

V_281 = Vertex(name = 'V_281',
               particles = [ P.d__tilde__, P.u, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_260,(0,1):C.GC_259})

V_282 = Vertex(name = 'V_282',
               particles = [ P.s__tilde__, P.u, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_261})

V_283 = Vertex(name = 'V_283',
               particles = [ P.b__tilde__, P.u, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_262})

V_284 = Vertex(name = 'V_284',
               particles = [ P.d__tilde__, P.c, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_263})

V_285 = Vertex(name = 'V_285',
               particles = [ P.s__tilde__, P.c, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_264,(0,1):C.GC_259})

V_286 = Vertex(name = 'V_286',
               particles = [ P.b__tilde__, P.c, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_265})

V_287 = Vertex(name = 'V_287',
               particles = [ P.d__tilde__, P.t, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_266})

V_288 = Vertex(name = 'V_288',
               particles = [ P.s__tilde__, P.t, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_267})

V_289 = Vertex(name = 'V_289',
               particles = [ P.b__tilde__, P.t, P.Z, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_268,(0,1):C.GC_259})

V_290 = Vertex(name = 'V_290',
               particles = [ P.u__tilde__, P.u, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_11})

V_291 = Vertex(name = 'V_291',
               particles = [ P.c__tilde__, P.c, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_11})

V_292 = Vertex(name = 'V_292',
               particles = [ P.t__tilde__, P.t, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_11})

V_293 = Vertex(name = 'V_293',
               particles = [ P.d__tilde__, P.d, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_97,(0,1):C.GC_120})

V_294 = Vertex(name = 'V_294',
               particles = [ P.s__tilde__, P.s, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_97,(0,1):C.GC_120})

V_295 = Vertex(name = 'V_295',
               particles = [ P.b__tilde__, P.b, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_97,(0,1):C.GC_120})

V_296 = Vertex(name = 'V_296',
               particles = [ P.d__tilde__, P.d, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_98,(0,1):C.GC_121})

V_297 = Vertex(name = 'V_297',
               particles = [ P.s__tilde__, P.s, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_98,(0,1):C.GC_121})

V_298 = Vertex(name = 'V_298',
               particles = [ P.b__tilde__, P.b, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_98,(0,1):C.GC_121})

V_299 = Vertex(name = 'V_299',
               particles = [ P.u__tilde__, P.d, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_310,(0,1):C.GC_125})

V_300 = Vertex(name = 'V_300',
               particles = [ P.c__tilde__, P.d, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_349})

V_301 = Vertex(name = 'V_301',
               particles = [ P.t__tilde__, P.d, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_388})

V_302 = Vertex(name = 'V_302',
               particles = [ P.u__tilde__, P.s, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_323})

V_303 = Vertex(name = 'V_303',
               particles = [ P.c__tilde__, P.s, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_362,(0,1):C.GC_125})

V_304 = Vertex(name = 'V_304',
               particles = [ P.t__tilde__, P.s, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_401})

V_305 = Vertex(name = 'V_305',
               particles = [ P.u__tilde__, P.b, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_336})

V_306 = Vertex(name = 'V_306',
               particles = [ P.c__tilde__, P.b, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_375})

V_307 = Vertex(name = 'V_307',
               particles = [ P.t__tilde__, P.b, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_414,(0,1):C.GC_125})

V_308 = Vertex(name = 'V_308',
               particles = [ P.u__tilde__, P.d, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_309,(0,1):C.GC_124})

V_309 = Vertex(name = 'V_309',
               particles = [ P.c__tilde__, P.d, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_348})

V_310 = Vertex(name = 'V_310',
               particles = [ P.t__tilde__, P.d, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_387})

V_311 = Vertex(name = 'V_311',
               particles = [ P.u__tilde__, P.s, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_322})

V_312 = Vertex(name = 'V_312',
               particles = [ P.c__tilde__, P.s, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_361,(0,1):C.GC_124})

V_313 = Vertex(name = 'V_313',
               particles = [ P.t__tilde__, P.s, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_400})

V_314 = Vertex(name = 'V_314',
               particles = [ P.u__tilde__, P.b, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_335})

V_315 = Vertex(name = 'V_315',
               particles = [ P.c__tilde__, P.b, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1 ],
               couplings = {(0,0):C.GC_374})

V_316 = Vertex(name = 'V_316',
               particles = [ P.t__tilde__, P.b, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFSS1, L.FFSS2 ],
               couplings = {(0,0):C.GC_413,(0,1):C.GC_124})

V_317 = Vertex(name = 'V_317',
               particles = [ P.d__tilde__, P.d, P.a, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_99,(0,1):C.GC_145})

V_318 = Vertex(name = 'V_318',
               particles = [ P.s__tilde__, P.s, P.a, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_99,(0,1):C.GC_145})

V_319 = Vertex(name = 'V_319',
               particles = [ P.b__tilde__, P.b, P.a, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_99,(0,1):C.GC_145})

V_320 = Vertex(name = 'V_320',
               particles = [ P.u__tilde__, P.d, P.a, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_312,(0,1):C.GC_148})

V_321 = Vertex(name = 'V_321',
               particles = [ P.c__tilde__, P.d, P.a, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_351})

V_322 = Vertex(name = 'V_322',
               particles = [ P.t__tilde__, P.d, P.a, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_390})

V_323 = Vertex(name = 'V_323',
               particles = [ P.u__tilde__, P.s, P.a, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_325})

V_324 = Vertex(name = 'V_324',
               particles = [ P.c__tilde__, P.s, P.a, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_364,(0,1):C.GC_148})

V_325 = Vertex(name = 'V_325',
               particles = [ P.t__tilde__, P.s, P.a, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_403})

V_326 = Vertex(name = 'V_326',
               particles = [ P.u__tilde__, P.b, P.a, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_338})

V_327 = Vertex(name = 'V_327',
               particles = [ P.c__tilde__, P.b, P.a, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_377})

V_328 = Vertex(name = 'V_328',
               particles = [ P.t__tilde__, P.b, P.a, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_416,(0,1):C.GC_148})

V_329 = Vertex(name = 'V_329',
               particles = [ P.u__tilde__, P.d, P.a, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_311,(0,1):C.GC_147})

V_330 = Vertex(name = 'V_330',
               particles = [ P.c__tilde__, P.d, P.a, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_350})

V_331 = Vertex(name = 'V_331',
               particles = [ P.t__tilde__, P.d, P.a, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_389})

V_332 = Vertex(name = 'V_332',
               particles = [ P.u__tilde__, P.s, P.a, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_324})

V_333 = Vertex(name = 'V_333',
               particles = [ P.c__tilde__, P.s, P.a, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_363,(0,1):C.GC_147})

V_334 = Vertex(name = 'V_334',
               particles = [ P.t__tilde__, P.s, P.a, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_402})

V_335 = Vertex(name = 'V_335',
               particles = [ P.u__tilde__, P.b, P.a, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_337})

V_336 = Vertex(name = 'V_336',
               particles = [ P.c__tilde__, P.b, P.a, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_376})

V_337 = Vertex(name = 'V_337',
               particles = [ P.t__tilde__, P.b, P.a, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_415,(0,1):C.GC_147})

V_338 = Vertex(name = 'V_338',
               particles = [ P.u__tilde__, P.d, P.a, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_317,(0,1):C.GC_227})

V_339 = Vertex(name = 'V_339',
               particles = [ P.c__tilde__, P.d, P.a, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_356})

V_340 = Vertex(name = 'V_340',
               particles = [ P.t__tilde__, P.d, P.a, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_395})

V_341 = Vertex(name = 'V_341',
               particles = [ P.u__tilde__, P.s, P.a, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_330})

V_342 = Vertex(name = 'V_342',
               particles = [ P.c__tilde__, P.s, P.a, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_369,(0,1):C.GC_227})

V_343 = Vertex(name = 'V_343',
               particles = [ P.t__tilde__, P.s, P.a, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_408})

V_344 = Vertex(name = 'V_344',
               particles = [ P.u__tilde__, P.b, P.a, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_343})

V_345 = Vertex(name = 'V_345',
               particles = [ P.c__tilde__, P.b, P.a, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_382})

V_346 = Vertex(name = 'V_346',
               particles = [ P.t__tilde__, P.b, P.a, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_421,(0,1):C.GC_227})

V_347 = Vertex(name = 'V_347',
               particles = [ P.d__tilde__, P.d, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_1})

V_348 = Vertex(name = 'V_348',
               particles = [ P.s__tilde__, P.s, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_1})

V_349 = Vertex(name = 'V_349',
               particles = [ P.b__tilde__, P.b, P.a ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_1})

V_350 = Vertex(name = 'V_350',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_313})

V_351 = Vertex(name = 'V_351',
               particles = [ P.c__tilde__, P.d, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_352})

V_352 = Vertex(name = 'V_352',
               particles = [ P.t__tilde__, P.d, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_391})

V_353 = Vertex(name = 'V_353',
               particles = [ P.u__tilde__, P.s, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_326})

V_354 = Vertex(name = 'V_354',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_365})

V_355 = Vertex(name = 'V_355',
               particles = [ P.t__tilde__, P.s, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_404})

V_356 = Vertex(name = 'V_356',
               particles = [ P.u__tilde__, P.b, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_339})

V_357 = Vertex(name = 'V_357',
               particles = [ P.c__tilde__, P.b, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_378})

V_358 = Vertex(name = 'V_358',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_417})

V_359 = Vertex(name = 'V_359',
               particles = [ P.d__tilde__, P.d, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_175,(0,1):C.GC_172})

V_360 = Vertex(name = 'V_360',
               particles = [ P.s__tilde__, P.s, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_175,(0,1):C.GC_172})

V_361 = Vertex(name = 'V_361',
               particles = [ P.b__tilde__, P.b, P.W__minus__, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_175,(0,1):C.GC_172})

V_362 = Vertex(name = 'V_362',
               particles = [ P.d__tilde__, P.d, P.W__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_174,(0,1):C.GC_171})

V_363 = Vertex(name = 'V_363',
               particles = [ P.s__tilde__, P.s, P.W__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_174,(0,1):C.GC_171})

V_364 = Vertex(name = 'V_364',
               particles = [ P.b__tilde__, P.b, P.W__minus__, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_174,(0,1):C.GC_171})

V_365 = Vertex(name = 'V_365',
               particles = [ P.d__tilde__, P.d, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_242,(0,1):C.GC_240})

V_366 = Vertex(name = 'V_366',
               particles = [ P.s__tilde__, P.s, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_242,(0,1):C.GC_240})

V_367 = Vertex(name = 'V_367',
               particles = [ P.b__tilde__, P.b, P.W__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_242,(0,1):C.GC_240})

V_368 = Vertex(name = 'V_368',
               particles = [ P.d__tilde__, P.d, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_175,(0,1):C.GC_172})

V_369 = Vertex(name = 'V_369',
               particles = [ P.s__tilde__, P.s, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_175,(0,1):C.GC_172})

V_370 = Vertex(name = 'V_370',
               particles = [ P.b__tilde__, P.b, P.W__plus__, P.G0, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_175,(0,1):C.GC_172})

V_371 = Vertex(name = 'V_371',
               particles = [ P.d__tilde__, P.d, P.W__plus__, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_176,(0,1):C.GC_173})

V_372 = Vertex(name = 'V_372',
               particles = [ P.s__tilde__, P.s, P.W__plus__, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_176,(0,1):C.GC_173})

V_373 = Vertex(name = 'V_373',
               particles = [ P.b__tilde__, P.b, P.W__plus__, P.G__minus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_176,(0,1):C.GC_173})

V_374 = Vertex(name = 'V_374',
               particles = [ P.d__tilde__, P.d, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_243,(0,1):C.GC_241})

V_375 = Vertex(name = 'V_375',
               particles = [ P.s__tilde__, P.s, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_243,(0,1):C.GC_241})

V_376 = Vertex(name = 'V_376',
               particles = [ P.b__tilde__, P.b, P.W__plus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_243,(0,1):C.GC_241})

V_377 = Vertex(name = 'V_377',
               particles = [ P.u__tilde__, P.d, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV3 ],
               couplings = {(0,0):[ C.GC_307, C.GC_308 ],(0,1):C.GC_56})

V_378 = Vertex(name = 'V_378',
               particles = [ P.c__tilde__, P.d, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):[ C.GC_346, C.GC_347 ]})

V_379 = Vertex(name = 'V_379',
               particles = [ P.t__tilde__, P.d, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):[ C.GC_385, C.GC_386 ]})

V_380 = Vertex(name = 'V_380',
               particles = [ P.u__tilde__, P.s, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):[ C.GC_320, C.GC_321 ]})

V_381 = Vertex(name = 'V_381',
               particles = [ P.c__tilde__, P.s, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV3 ],
               couplings = {(0,0):[ C.GC_359, C.GC_360 ],(0,1):C.GC_56})

V_382 = Vertex(name = 'V_382',
               particles = [ P.t__tilde__, P.s, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):[ C.GC_398, C.GC_399 ]})

V_383 = Vertex(name = 'V_383',
               particles = [ P.u__tilde__, P.b, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):[ C.GC_333, C.GC_334 ]})

V_384 = Vertex(name = 'V_384',
               particles = [ P.c__tilde__, P.b, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2 ],
               couplings = {(0,0):[ C.GC_372, C.GC_373 ]})

V_385 = Vertex(name = 'V_385',
               particles = [ P.t__tilde__, P.b, P.W__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV3 ],
               couplings = {(0,0):[ C.GC_411, C.GC_412 ],(0,1):C.GC_56})

V_386 = Vertex(name = 'V_386',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_313,(0,1):C.GC_181})

V_387 = Vertex(name = 'V_387',
               particles = [ P.c__tilde__, P.d, P.W__plus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_352})

V_388 = Vertex(name = 'V_388',
               particles = [ P.t__tilde__, P.d, P.W__plus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_391})

V_389 = Vertex(name = 'V_389',
               particles = [ P.u__tilde__, P.s, P.W__plus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_326})

V_390 = Vertex(name = 'V_390',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_365,(0,1):C.GC_181})

V_391 = Vertex(name = 'V_391',
               particles = [ P.t__tilde__, P.s, P.W__plus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_404})

V_392 = Vertex(name = 'V_392',
               particles = [ P.u__tilde__, P.b, P.W__plus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_339})

V_393 = Vertex(name = 'V_393',
               particles = [ P.c__tilde__, P.b, P.W__plus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_378})

V_394 = Vertex(name = 'V_394',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_417,(0,1):C.GC_181})

V_395 = Vertex(name = 'V_395',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_313,(0,1):C.GC_182})

V_396 = Vertex(name = 'V_396',
               particles = [ P.c__tilde__, P.d, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_352})

V_397 = Vertex(name = 'V_397',
               particles = [ P.t__tilde__, P.d, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_391})

V_398 = Vertex(name = 'V_398',
               particles = [ P.u__tilde__, P.s, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_326})

V_399 = Vertex(name = 'V_399',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_365,(0,1):C.GC_182})

V_400 = Vertex(name = 'V_400',
               particles = [ P.t__tilde__, P.s, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_404})

V_401 = Vertex(name = 'V_401',
               particles = [ P.u__tilde__, P.b, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_339})

V_402 = Vertex(name = 'V_402',
               particles = [ P.c__tilde__, P.b, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_378})

V_403 = Vertex(name = 'V_403',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_417,(0,1):C.GC_182})

V_404 = Vertex(name = 'V_404',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_318,(0,1):C.GC_247})

V_405 = Vertex(name = 'V_405',
               particles = [ P.c__tilde__, P.d, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_357})

V_406 = Vertex(name = 'V_406',
               particles = [ P.t__tilde__, P.d, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_396})

V_407 = Vertex(name = 'V_407',
               particles = [ P.u__tilde__, P.s, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_331})

V_408 = Vertex(name = 'V_408',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_370,(0,1):C.GC_247})

V_409 = Vertex(name = 'V_409',
               particles = [ P.t__tilde__, P.s, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_409})

V_410 = Vertex(name = 'V_410',
               particles = [ P.u__tilde__, P.b, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_344})

V_411 = Vertex(name = 'V_411',
               particles = [ P.c__tilde__, P.b, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_383})

V_412 = Vertex(name = 'V_412',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_422,(0,1):C.GC_247})

V_413 = Vertex(name = 'V_413',
               particles = [ P.d__tilde__, P.d, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_110,(0,1):C.GC_105})

V_414 = Vertex(name = 'V_414',
               particles = [ P.s__tilde__, P.s, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_110,(0,1):C.GC_105})

V_415 = Vertex(name = 'V_415',
               particles = [ P.b__tilde__, P.b, P.Z, P.G__minus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_110,(0,1):C.GC_105})

V_416 = Vertex(name = 'V_416',
               particles = [ P.d__tilde__, P.d, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV3, L.FFV4 ],
               couplings = {(0,0):[ C.GC_75, C.GC_90 ],(0,2):C.GC_83,(0,1):C.GC_88})

V_417 = Vertex(name = 'V_417',
               particles = [ P.s__tilde__, P.s, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV3, L.FFV4 ],
               couplings = {(0,0):[ C.GC_75, C.GC_90 ],(0,2):C.GC_83,(0,1):C.GC_88})

V_418 = Vertex(name = 'V_418',
               particles = [ P.b__tilde__, P.b, P.Z ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFV2, L.FFV3, L.FFV4 ],
               couplings = {(0,0):[ C.GC_75, C.GC_90 ],(0,2):C.GC_83,(0,1):C.GC_88})

V_419 = Vertex(name = 'V_419',
               particles = [ P.d__tilde__, P.d, P.Z, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_108,(0,1):C.GC_104})

V_420 = Vertex(name = 'V_420',
               particles = [ P.s__tilde__, P.s, P.Z, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_108,(0,1):C.GC_104})

V_421 = Vertex(name = 'V_421',
               particles = [ P.b__tilde__, P.b, P.Z, P.G0, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_108,(0,1):C.GC_104})

V_422 = Vertex(name = 'V_422',
               particles = [ P.d__tilde__, P.d, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_108,(0,1):C.GC_104})

V_423 = Vertex(name = 'V_423',
               particles = [ P.s__tilde__, P.s, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_108,(0,1):C.GC_104})

V_424 = Vertex(name = 'V_424',
               particles = [ P.b__tilde__, P.b, P.Z, P.H, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_108,(0,1):C.GC_104})

V_425 = Vertex(name = 'V_425',
               particles = [ P.d__tilde__, P.d, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_118,(0,1):C.GC_116})

V_426 = Vertex(name = 'V_426',
               particles = [ P.s__tilde__, P.s, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_118,(0,1):C.GC_116})

V_427 = Vertex(name = 'V_427',
               particles = [ P.b__tilde__, P.b, P.Z, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_118,(0,1):C.GC_116})

V_428 = Vertex(name = 'V_428',
               particles = [ P.u__tilde__, P.d, P.Z, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_314,(0,1):C.GC_195})

V_429 = Vertex(name = 'V_429',
               particles = [ P.c__tilde__, P.d, P.Z, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_353})

V_430 = Vertex(name = 'V_430',
               particles = [ P.t__tilde__, P.d, P.Z, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_392})

V_431 = Vertex(name = 'V_431',
               particles = [ P.u__tilde__, P.s, P.Z, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_327})

V_432 = Vertex(name = 'V_432',
               particles = [ P.c__tilde__, P.s, P.Z, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_366,(0,1):C.GC_195})

V_433 = Vertex(name = 'V_433',
               particles = [ P.t__tilde__, P.s, P.Z, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_405})

V_434 = Vertex(name = 'V_434',
               particles = [ P.u__tilde__, P.b, P.Z, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_340})

V_435 = Vertex(name = 'V_435',
               particles = [ P.c__tilde__, P.b, P.Z, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_379})

V_436 = Vertex(name = 'V_436',
               particles = [ P.t__tilde__, P.b, P.Z, P.G0, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_418,(0,1):C.GC_195})

V_437 = Vertex(name = 'V_437',
               particles = [ P.u__tilde__, P.d, P.Z, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_315,(0,1):C.GC_194})

V_438 = Vertex(name = 'V_438',
               particles = [ P.c__tilde__, P.d, P.Z, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_354})

V_439 = Vertex(name = 'V_439',
               particles = [ P.t__tilde__, P.d, P.Z, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_393})

V_440 = Vertex(name = 'V_440',
               particles = [ P.u__tilde__, P.s, P.Z, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_328})

V_441 = Vertex(name = 'V_441',
               particles = [ P.c__tilde__, P.s, P.Z, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_367,(0,1):C.GC_194})

V_442 = Vertex(name = 'V_442',
               particles = [ P.t__tilde__, P.s, P.Z, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_406})

V_443 = Vertex(name = 'V_443',
               particles = [ P.u__tilde__, P.b, P.Z, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_341})

V_444 = Vertex(name = 'V_444',
               particles = [ P.c__tilde__, P.b, P.Z, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1 ],
               couplings = {(0,0):C.GC_380})

V_445 = Vertex(name = 'V_445',
               particles = [ P.t__tilde__, P.b, P.Z, P.G__plus__, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS1, L.FFVSS2 ],
               couplings = {(0,0):C.GC_419,(0,1):C.GC_194})

V_446 = Vertex(name = 'V_446',
               particles = [ P.u__tilde__, P.d, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_319,(0,1):C.GC_258})

V_447 = Vertex(name = 'V_447',
               particles = [ P.c__tilde__, P.d, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_358})

V_448 = Vertex(name = 'V_448',
               particles = [ P.t__tilde__, P.d, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_397})

V_449 = Vertex(name = 'V_449',
               particles = [ P.u__tilde__, P.s, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_332})

V_450 = Vertex(name = 'V_450',
               particles = [ P.c__tilde__, P.s, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_371,(0,1):C.GC_258})

V_451 = Vertex(name = 'V_451',
               particles = [ P.t__tilde__, P.s, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_410})

V_452 = Vertex(name = 'V_452',
               particles = [ P.u__tilde__, P.b, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_345})

V_453 = Vertex(name = 'V_453',
               particles = [ P.c__tilde__, P.b, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1 ],
               couplings = {(0,0):C.GC_384})

V_454 = Vertex(name = 'V_454',
               particles = [ P.t__tilde__, P.b, P.Z, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS1, L.FFVS2 ],
               couplings = {(0,0):C.GC_423,(0,1):C.GC_258})

V_455 = Vertex(name = 'V_455',
               particles = [ P.d__tilde__, P.d, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_11})

V_456 = Vertex(name = 'V_456',
               particles = [ P.s__tilde__, P.s, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_11})

V_457 = Vertex(name = 'V_457',
               particles = [ P.b__tilde__, P.b, P.g ],
               color = [ 'T(3,2,1)' ],
               lorentz = [ L.FFV1 ],
               couplings = {(0,0):C.GC_11})

V_458 = Vertex(name = 'V_458',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_180})

V_459 = Vertex(name = 'V_459',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_180})

V_460 = Vertex(name = 'V_460',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_180})

V_461 = Vertex(name = 'V_461',
               particles = [ P.d__tilde__, P.u, P.W__minus__, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2 ],
               couplings = {(0,0):C.GC_246})

V_462 = Vertex(name = 'V_462',
               particles = [ P.s__tilde__, P.c, P.W__minus__, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2 ],
               couplings = {(0,0):C.GC_246})

V_463 = Vertex(name = 'V_463',
               particles = [ P.b__tilde__, P.t, P.W__minus__, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2 ],
               couplings = {(0,0):C.GC_246})

V_464 = Vertex(name = 'V_464',
               particles = [ P.d__tilde__, P.u, P.W__plus__, P.G__minus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_183})

V_465 = Vertex(name = 'V_465',
               particles = [ P.s__tilde__, P.c, P.W__plus__, P.G__minus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_183})

V_466 = Vertex(name = 'V_466',
               particles = [ P.b__tilde__, P.t, P.W__plus__, P.G__minus__, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_183})

V_467 = Vertex(name = 'V_467',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_184})

V_468 = Vertex(name = 'V_468',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_184})

V_469 = Vertex(name = 'V_469',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.G0, P.H ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_184})

V_470 = Vertex(name = 'V_470',
               particles = [ P.u__tilde__, P.d, P.W__plus__, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2 ],
               couplings = {(0,0):C.GC_248})

V_471 = Vertex(name = 'V_471',
               particles = [ P.c__tilde__, P.s, P.W__plus__, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2 ],
               couplings = {(0,0):C.GC_248})

V_472 = Vertex(name = 'V_472',
               particles = [ P.t__tilde__, P.b, P.W__plus__, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVS2 ],
               couplings = {(0,0):C.GC_248})

V_473 = Vertex(name = 'V_473',
               particles = [ P.u__tilde__, P.d, P.W__minus__, P.G__plus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_183})

V_474 = Vertex(name = 'V_474',
               particles = [ P.c__tilde__, P.s, P.W__minus__, P.G__plus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_183})

V_475 = Vertex(name = 'V_475',
               particles = [ P.t__tilde__, P.b, P.W__minus__, P.G__plus__, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFVSS2 ],
               couplings = {(0,0):C.GC_183})

