# This file was automatically created by FeynRules $Revision: 436 $
# Mathematica version: 7.0 for Linux x86 (32-bit) (February 18, 2009)
# Date: Thu 19 Apr 2012 11:28:04


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec



GC_1 = Coupling(name = 'GC_1',
                value = '-G',
                order = {'QCD':1})

GC_2 = Coupling(name = 'GC_2',
                value = 'complex(0,1)*G**2',
                order = {'QCD':2})

GC_3 = Coupling(name = 'GC_3',
                value = '(4*complex(0,1)*G*KapS)/LamNP',
                order = {'QCD':1})

GC_4 = Coupling(name = 'GC_4',
                value = '(4*G**2*KapS)/LamNP',
                order = {'QCD':2})

GC_5 = Coupling(name = 'GC_5',
                value = '(-2*complex(0,1)*G**3*KapS)/LamNP',
                order = {'QCD':3})

