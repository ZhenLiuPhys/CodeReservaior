# This file was automatically created by FeynRules 2.3.4
# Mathematica version: 10.2.0 for Mac OS X x86 (64-bit) (July 29, 2015)
# Date: Tue 13 Oct 2015 14:49:19



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.9,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

Gf = Parameter(name = 'Gf',
               nature = 'external',
               type = 'real',
               value = 0.0000116637,
               texname = 'G_f',
               lhablock = 'SMINPUTS',
               lhacode = [ 2 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.1184,
               texname = '\\alpha _s',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

ymb = Parameter(name = 'ymb',
                nature = 'external',
                type = 'real',
                value = 4.7,
                texname = '\\text{ymb}',
                lhablock = 'YUKAWA',
                lhacode = [ 5 ])

ymt = Parameter(name = 'ymt',
                nature = 'external',
                type = 'real',
                value = 172,
                texname = '\\text{ymt}',
                lhablock = 'YUKAWA',
                lhacode = [ 6 ])

ymtau = Parameter(name = 'ymtau',
                  nature = 'external',
                  type = 'real',
                  value = 1.777,
                  texname = '\\text{ymtau}',
                  lhablock = 'YUKAWA',
                  lhacode = [ 15 ])

gR = Parameter(name = 'gR',
               nature = 'external',
               type = 'real',
               value = 0.4,
               texname = '\\text{gR}',
               lhablock = 'FRBlock',
               lhacode = [ 1 ])

gYZp = Parameter(name = 'gYZp',
                 nature = 'external',
                 type = 'real',
                 value = 0.363,
                 texname = '\\text{gYZp}',
                 lhablock = 'FRBlock',
                 lhacode = [ 2 ])

gLWp = Parameter(name = 'gLWp',
                 nature = 'external',
                 type = 'real',
                 value = 0.629,
                 texname = '\\text{gLWp}',
                 lhablock = 'FRBlock',
                 lhacode = [ 3 ])

V1e = Parameter(name = 'V1e',
                nature = 'external',
                type = 'real',
                value = 1.,
                texname = '\\text{V1e}',
                lhablock = 'FRBlock',
                lhacode = [ 4 ])

V1mu = Parameter(name = 'V1mu',
                 nature = 'external',
                 type = 'real',
                 value = 1.,
                 texname = '\\text{V1mu}',
                 lhablock = 'FRBlock',
                 lhacode = [ 5 ])

V1tau = Parameter(name = 'V1tau',
                  nature = 'external',
                  type = 'real',
                  value = 1.,
                  texname = '\\text{V1tau}',
                  lhablock = 'FRBlock',
                  lhacode = [ 6 ])

tbeta = Parameter(name = 'tbeta',
                  nature = 'external',
                  type = 'real',
                  value = 1.,
                  texname = '\\text{tbeta}',
                  lhablock = 'FRBlock',
                  lhacode = [ 7 ])

MZ = Parameter(name = 'MZ',
               nature = 'external',
               type = 'real',
               value = 91.1876,
               texname = '\\text{MZ}',
               lhablock = 'MASS',
               lhacode = [ 23 ])

MZp = Parameter(name = 'MZp',
                nature = 'external',
                type = 'real',
                value = 4500.,
                texname = '\\text{MZp}',
                lhablock = 'MASS',
                lhacode = [ 9906663 ])

MWp = Parameter(name = 'MWp',
                nature = 'external',
                type = 'real',
                value = 1900.,
                texname = '\\text{MWp}',
                lhablock = 'MASS',
                lhacode = [ 9916663 ])

MTA = Parameter(name = 'MTA',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{MTA}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 173,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.7,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

MN1 = Parameter(name = 'MN1',
                nature = 'external',
                type = 'real',
                value = 3500.,
                texname = '\\text{MN1}',
                lhablock = 'MASS',
                lhacode = [ 9926662 ])

MND = Parameter(name = 'MND',
                nature = 'external',
                type = 'real',
                value = 1300.,
                texname = '\\text{MND}',
                lhablock = 'MASS',
                lhacode = [ 9936662 ])

MH = Parameter(name = 'MH',
               nature = 'external',
               type = 'real',
               value = 125,
               texname = '\\text{MH}',
               lhablock = 'MASS',
               lhacode = [ 25 ])

MHzero = Parameter(name = 'MHzero',
                   nature = 'external',
                   type = 'real',
                   value = 500.,
                   texname = '\\text{MHzero}',
                   lhablock = 'MASS',
                   lhacode = [ 9936660 ])

MAzero = Parameter(name = 'MAzero',
                   nature = 'external',
                   type = 'real',
                   value = 500.,
                   texname = '\\text{MAzero}',
                   lhablock = 'MASS',
                   lhacode = [ 9946660 ])

MHplus = Parameter(name = 'MHplus',
                   nature = 'external',
                   type = 'real',
                   value = 500.,
                   texname = '\\text{MHplus}',
                   lhablock = 'MASS',
                   lhacode = [ 9956660 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.4414,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.085,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

WZp = Parameter(name = 'WZp',
                nature = 'external',
                type = 'real',
                value = 131.86843257876606,
                texname = '\\text{WZp}',
                lhablock = 'DECAY',
                lhacode = [ 9906663 ])

WWp = Parameter(name = 'WWp',
                nature = 'external',
                type = 'real',
                value = 21.119082230941864,
                texname = '\\text{WWp}',
                lhablock = 'DECAY',
                lhacode = [ 9916663 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.50833649,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

WN1 = Parameter(name = 'WN1',
                nature = 'external',
                type = 'real',
                value = 22.438057020405854,
                texname = '\\text{WN1}',
                lhablock = 'DECAY',
                lhacode = [ 9926662 ])

WND = Parameter(name = 'WND',
                nature = 'external',
                type = 'real',
                value = 0.0006913622648545521,
                texname = '\\text{WND}',
                lhablock = 'DECAY',
                lhacode = [ 9936662 ])

WH = Parameter(name = 'WH',
               nature = 'external',
               type = 'real',
               value = 0.004070000000000001,
               texname = '\\text{WH}',
               lhablock = 'DECAY',
               lhacode = [ 25 ])

WHzero = Parameter(name = 'WHzero',
                   nature = 'external',
                   type = 'real',
                   value = 22.169278164048766,
                   texname = '\\text{WHzero}',
                   lhablock = 'DECAY',
                   lhacode = [ 9936660 ])

WAzero = Parameter(name = 'WAzero',
                   nature = 'external',
                   type = 'real',
                   value = 42.54029305987068,
                   texname = '\\text{WAzero}',
                   lhablock = 'DECAY',
                   lhacode = [ 9946660 ])

WHplus = Parameter(name = 'WHplus',
                   nature = 'external',
                   type = 'real',
                   value = 22.834296721546398,
                   texname = '\\text{WHplus}',
                   lhablock = 'DECAY',
                   lhacode = [ 9956660 ])

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\alpha _{\\text{EW}}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

MW = Parameter(name = 'MW',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2))))',
               texname = 'M_W')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname = 'e')

sw2 = Parameter(name = 'sw2',
                nature = 'internal',
                type = 'real',
                value = '1 - MW**2/MZ**2',
                texname = '\\text{sw2}')

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - sw2)',
               texname = 'c_w')

sw = Parameter(name = 'sw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(sw2)',
               texname = 's_w')

g1 = Parameter(name = 'g1',
               nature = 'internal',
               type = 'real',
               value = 'ee/cw',
               texname = 'g_1')

gw = Parameter(name = 'gw',
               nature = 'internal',
               type = 'real',
               value = 'ee/sw',
               texname = 'g_w')

vev = Parameter(name = 'vev',
                nature = 'internal',
                type = 'real',
                value = '(2*MW*sw)/ee',
                texname = '\\text{vev}')

lam = Parameter(name = 'lam',
                nature = 'internal',
                type = 'real',
                value = 'MH**2/(2.*vev**2)',
                texname = '\\text{lam}')

yb = Parameter(name = 'yb',
               nature = 'internal',
               type = 'real',
               value = '(ymb*cmath.sqrt(2))/vev',
               texname = '\\text{yb}')

yt = Parameter(name = 'yt',
               nature = 'internal',
               type = 'real',
               value = '(ymt*cmath.sqrt(2))/vev',
               texname = '\\text{yt}')

ytau = Parameter(name = 'ytau',
                 nature = 'internal',
                 type = 'real',
                 value = '(ymtau*cmath.sqrt(2))/vev',
                 texname = '\\text{ytau}')

muH = Parameter(name = 'muH',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(lam*vev**2)',
                texname = '\\mu')

