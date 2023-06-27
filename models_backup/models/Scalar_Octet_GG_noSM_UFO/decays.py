# This file was automatically created by FeynRules 1.7.221
# Mathematica version: 9.0 for Linux x86 (64-bit) (November 20, 2012)
# Date: Tue 28 Apr 2015 17:36:04


from object_library import all_decays, Decay
import particles as P


Decay_S8 = Decay(name = 'Decay_S8',
                 particle = P.S8,
                 partial_widths = {(P.G,P.G):'(G**2*KapS**2*MS8**6*dSUN(2,3,1)**2)/(32.*cmath.pi*LamNP**2*abs(MS8)**3)'})

