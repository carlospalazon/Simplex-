## Solucions del simplex

### Problema inicial 
```
c=
 -85 -65 -34  92 -83  93  22  76  48 -71  25 -50  79 -99   0   0   0   0   0   0

A=
 -13 -29  23 -27  57 -52 -70  81  -2  64 -54   9  68  10   0   0   0   0   0   0
  -6 -53 -37  36  85 -22  35  87  32  34  65 -64  65 -25   0   0   0   0   0   0
 -13  93  31  16  45  23 -11  83  75 -65  67  10 -92  88   0   0   0   0   0   0
  67  -8 -27 -22  47  44  -3 -23 -16  81 -24 -38 -67  24   0   0   0   0   0   0
  80  57  77  72  52  58  53  84  77  92  57  75  68  88   1   0   0   0   0   0
 -81  21 -14  73  11 -10  -7  45  95 -92  97 -79  89 -92   0  -1   0   0   0   0
  96  -7  85 -63 -24 -95  42 -71 100  42  43 -23  83  60   0   0  -1   0   0   0
 -27  92 -84  77  68  73  28  53  34 -71  59 -17 -36  -1   0   0   0   1   0   0
  -1 -28  96  74  66 -82 -24  91  63  65  64 -29  73  93   0   0   0   0  -1   0
 -60 -34  21 -65  84   4  12  98 -38  -1 -99  25  21  82   0   0   0   0   0  -1

b=
  65 232 350  35 991  55 267 249 520  49

```
#### Output

```
[jh_simplexP] Inici simplex primal amb regla de Bland
[jh_simplexP]   Fase I
[jh_simplexP]     Iteracio    1 : iout =   4, q =   1, B(p) =  24, theta*=  0.522, z = 2813.000
[jh_simplexP]     Iteracio    2 : iout =   8, q =   2, B(p) =  28, theta*=  2.964, z = 2791.060
[jh_simplexP]     Iteracio    3 : iout =   3, q =   3, B(p) =  23, theta*=  0.694, z = 2467.974
[jh_simplexP]     Iteracio    4 : iout =   6, q =   4, B(p) =  26, theta*=  4.106, z = 2256.543
[jh_simplexP]     Iteracio    5 : iout =   7, q =   7, B(p) =  27, theta*=  0.070, z = 1249.582
[jh_simplexP]     Iteracio    6 : iout =  10, q =   5, B(p) =  30, theta*=  2.323, z = 1221.248
[jh_simplexP]     Iteracio    7 : iout =   5, q =   6, B(p) =  25, theta*=  1.050, z =  498.991
[jh_simplexP]     Iteracio    8 : iout =   9, q =   8, B(p) =  29, theta*=  0.049, z =  335.717
[jh_simplexP]     Iteracio    9 : iout =   2, q =   9, B(p) =  22, theta*=  0.227, z =  323.301
[jh_simplexP]     Iteracio   10 : iout =   5, q =  10, B(p) =   6, theta*=  1.846, z =  188.311
[jh_simplexP]     Iteracio   11 : iout =   1, q =  12, B(p) =  21, theta*=  1.456, z =   87.373
[jh_simplexP]     Iteracio   12 : iout =   -, q =   -, B(p) =   -, theta*=     -, z =    0.000
[jh_simplexP]     Solucio basica factible trobada, iteracio 12
[jh_simplexP]   Fase II
[jh_simplexP]     Iteracio   13 : iout =   4, q =  11, B(p) =   1, theta*=  0.180, z = -128.454
[jh_simplexP]     Iteracio   14 : iout =   9, q =   6, B(p) =   8, theta*=  0.258, z = -153.969
[jh_simplexP]     Iteracio   15 : iout =   2, q =  13, B(p) =   9, theta*=  0.820, z = -196.429
[jh_simplexP]     Iteracio   16 : iout =   9, q =  15, B(p) =   6, theta*= 100.881, z = -373.024
[jh_simplexP]     Iteracio   17 : iout =   9, q =  17, B(p) =  15, theta*= 118.527, z = -420.669
[jh_simplexP]     Iteracio   18 : iout =   5, q =   1, B(p) =  10, theta*=  2.047, z = -485.604
[jh_simplexP]     Iteracio   19 : iout =   6, q =  14, B(p) =   4, theta*=  0.713, z = -501.781
[jh_simplexP]     Iteracio   20 : iout =   6, q =  18, B(p) =  14, theta*= 40.301, z = -520.870
[jh_simplexP]     Iteracio   21 : iout =   -, q =   -, B(p) =   -, theta*=     -, z = -533.913
[jh_simplexP]     Solucio optima trobada, iteracio  21, z = -533.912841
[jh_simplexP] Fi simplex primal 

Solucio optima:

vb = 12 13  3 11  1 18  7  2 17  5
xb =    0.9    1.5    2.5    1.4    2.1   40.3    1.6    2.2  289.6    3.2
z = -533.9128
r =    0.0    0.0   -0.0   65.9   -0.0  152.9    0.0  147.7   73.1   23.1   -0.0    0.0    0.0   18.3    0.4    0.3    0.0    0.0    0.2    0.0
```
### Problema 1 conjunt 34

```
c=
   8  -5   6  77 -77  96 -74 -37  91 -100 -10 -59  28  85   0   0   0   0   0   0

A=
 -67  13 -26  24 -32  36 -36  46  11  77 -64 -86  90  42   0   0   0   0   0   0
 -42 -70  30 -74  59  35 -13   6  86 -54  12  91  68  82   0   0   0   0   0   0
  97  80 -70  45  45 -34  22  55 -39 -20  70  31  58 -27   0   0   0   0   0   0
 -45 -33 -35  44  34   6 100 -65  30  56  38 -25  -5 -46   0   0   0   0   0   0
  59  60  92  58  83  75  74  87  70  59  63  97  52  85   1   0   0   0   0   0
 -40  89  -9  79  32   3  44 -36  25 -86  30 -33  24 -84   0   1   0   0   0   0
  82 -24 -100  30 -74   1  95 -57 -49  55 -89  90  45  50   0   0  -1   0   0   0
  92  97 -44  35  35  48  77 -17  -9  74  16  70 -59  29   0   0   0   1   0   0
  70  22   2  69   4  92 -74 -83  56 -69  88  69  32  69   0   0   0   0   1   0
 -42  30  98  94   8 -46  70 -50 -14 -38  30  51  86  25   0   0   0   0   0   1

b=
  28 216 313  54 1015  39  54 445 348 303
```

#### Output

```
[jh_simplexP] Inici simplex primal amb regla de Bland
[jh_simplexP]   Fase I
[jh_simplexP]     Iteracio    1 : iout =   7, q =   1, B(p) =  27, theta*=  0.659, z = 2815.000
[jh_simplexP]     Iteracio    2 : iout =   6, q =   2, B(p) =  26, theta*=  0.845, z = 2707.000
[jh_simplexP]     Iteracio    3 : iout =   3, q =   3, B(p) =  23, theta*=  1.218, z = 2443.242
[jh_simplexP]     Iteracio    4 : iout =   9, q =   4, B(p) =  29, theta*=  1.161, z = 1991.135
[jh_simplexP]     Iteracio    5 : iout =   6, q =   5, B(p) =   2, theta*=  1.134, z = 1624.557
[jh_simplexP]     Iteracio    6 : iout =   3, q =   6, B(p) =   3, theta*=  1.349, z = 1468.862
[jh_simplexP]     Iteracio    7 : iout =   8, q =   7, B(p) =  28, theta*=  0.383, z = 1415.671
[jh_simplexP]     Iteracio    8 : iout =  10, q =   3, B(p) =  30, theta*=  2.740, z = 1222.967
[jh_simplexP]     Iteracio    9 : iout =   5, q =   8, B(p) =  25, theta*=  1.514, z =  957.663
[jh_simplexP]     Iteracio   10 : iout =   3, q =   9, B(p) =   6, theta*=  2.064, z =  735.863
[jh_simplexP]     Iteracio   11 : iout =   4, q =  10, B(p) =  24, theta*=  1.571, z =  682.006
[jh_simplexP]     Iteracio   12 : iout =   1, q =   2, B(p) =  21, theta*=  0.057, z =  474.694
[jh_simplexP]     Iteracio   13 : iout =   1, q =   6, B(p) =   2, theta*=  0.150, z =  463.824
[jh_simplexP]     Iteracio   14 : iout =   8, q =  12, B(p) =   7, theta*=  0.717, z =  463.043
[jh_simplexP]     Iteracio   15 : iout =   3, q =  11, B(p) =   9, theta*=  0.654, z =  364.841
[jh_simplexP]     Iteracio   16 : iout =   7, q =   2, B(p) =   1, theta*=  0.104, z =  187.990
[jh_simplexP]     Iteracio   17 : iout =   6, q =   9, B(p) =   5, theta*=  1.050, z =  157.658
[jh_simplexP]     Iteracio   18 : iout =   2, q =  13, B(p) =  22, theta*=  0.714, z =  156.481
[jh_simplexP]     Iteracio   19 : iout =   -, q =   -, B(p) =   -, theta*=     -, z =    0.000
[jh_simplexP]     Solucio basica factible trobada, iteracio 19
[jh_simplexP]   Fase II
[jh_simplexP]     Iteracio   20 : iout =   6, q =   5, B(p) =   9, theta*=  1.706, z =  -77.097
[jh_simplexP]     Iteracio   21 : iout =   5, q =  15, B(p) =   8, theta*= 179.166, z = -132.510
[jh_simplexP]     Iteracio   22 : iout =   9, q =  16, B(p) =   4, theta*=  9.381, z = -245.214
[jh_simplexP]     Iteracio   23 : iout =   5, q =   1, B(p) =  15, theta*=  1.397, z = -251.992
[jh_simplexP]     Iteracio   24 : iout =   6, q =   8, B(p) =   5, theta*=  0.702, z = -296.160
[jh_simplexP]     Iteracio   25 : iout =   7, q =  17, B(p) =   2, theta*= 273.616, z = -323.140
[jh_simplexP]     Iteracio   26 : iout =   1, q =  19, B(p) =   6, theta*= 147.251, z = -363.300
[jh_simplexP]     Iteracio   27 : iout =   5, q =   2, B(p) =   1, theta*=  0.761, z = -531.478
[jh_simplexP]     Iteracio   28 : iout =   7, q =   5, B(p) =  17, theta*=  0.832, z = -580.361
[jh_simplexP]     Iteracio   29 : iout =   3, q =   7, B(p) =  11, theta*=  1.360, z = -588.083
[jh_simplexP]     Iteracio   30 : iout =   7, q =  17, B(p) =   5, theta*= 356.161, z = -630.546
[jh_simplexP]     Iteracio   31 : iout =   5, q =  20, B(p) =   2, theta*= 104.747, z = -639.085
[jh_simplexP]     Iteracio   32 : iout =   2, q =   5, B(p) =  13, theta*=  2.203, z = -695.736
[jh_simplexP]     Iteracio   33 : iout =   7, q =   1, B(p) =  17, theta*=  0.004, z = -787.471
[jh_simplexP]     Iteracio   34 : iout =   -, q =   -, B(p) =   -, theta*=     -, z = -787.487
[jh_simplexP]     Solucio optima trobada, iteracio  34, z = -787.486630
[jh_simplexP] Fi simplex primal 

Solucio optima:

vb = 19  5  7 10 20  8  1 12 16  3
xb =  764.4    2.2    1.5    2.5  291.8    3.6    0.0    2.2  322.8    0.5
z = -787.4866
r =    0.0   31.4    0.0  110.0   -0.0  192.2   -0.0    0.0  176.1   -0.0    5.5    0.0  114.3  186.0    0.4    0.0    0.0    0.5    0.0    0.0
```
### Problema 2 conjunt 34

```
c=
 -82 -63  73  15 -51  69 100  94  61  81  18  41 -33 -46   0   0   0   0   0   0 

A=
  53 -61 -53  23   7  17 -23 -87  64  18  73  33 -20 -38   0   0   0   0   0   0 
 -45 -81  63 -85  38 -67  87  53 -61  41  11  66 -23  82   0   0   0   0   0   0 
 -74  91  16   4 -31 -42  51  38 -75  94  69 -51  81  76   0   0   0   0   0   0 
 -65 -75  85 -22 -76  95 -13   2  57  68  10  78  57 -54   0   0   0   0   0   0 
  69  57  94  77  98  63  74  51  51  74  58  53  78  68   1   0   0   0   0   0 
  65 -78  75 -14  38 -23 -59 -29 -67  61  74  -6  -6  75   0  -1   0   0   0   0 
 -38 -21 -10  10 -66  75 -11  34  37  17  74  10  75  99   0   0  -1   0   0   0 
  58  27 -10 -23  60   2  72 -90  47 -77 -67  24  64 -80   0   0   0   1   0   0 
  26 -33  66   4  29  11  62  31  47 -54  79  89 -98 -37   0   0   0   0   1   0 
  20  37  96  28 -95  32 -31  43  41  49 -22 -90  68 -11   0   0   0   0   0   1 

b=
   6  79 247 147 966 105 284   8 223 166 
```
#### Output

```
[jh_simplexP] Inici simplex primal amb regla de Bland
[jh_simplexP]   Fase I
[jh_simplexP]     Iteracio    1 : iout =   1, q =   1, B(p) =  21, theta*=  0.113, z = 2231.000
[jh_simplexP]     Iteracio    2 : iout =   8, q =   3, B(p) =  28, theta*=  0.030, z = 2223.189
[jh_simplexP]     Iteracio    3 : iout =   6, q =   4, B(p) =  26, theta*=  0.951, z = 2208.520
[jh_simplexP]     Iteracio    4 : iout =  10, q =   2, B(p) =  30, theta*=  0.145, z = 1766.561
[jh_simplexP]     Iteracio    5 : iout =   5, q =   5, B(p) =  25, theta*=  2.494, z = 1724.330
[jh_simplexP]     Iteracio    6 : iout =   1, q =   6, B(p) =   1, theta*=  0.132, z = 1441.831
[jh_simplexP]     Iteracio    7 : iout =   1, q =   7, B(p) =   6, theta*=  0.108, z = 1412.946
[jh_simplexP]     Iteracio    8 : iout =   9, q =   8, B(p) =  29, theta*=  0.008, z = 1393.373
[jh_simplexP]     Iteracio    9 : iout =   1, q =   6, B(p) =   7, theta*=  0.666, z = 1379.835
[jh_simplexP]     Iteracio   10 : iout =   4, q =  10, B(p) =  24, theta*=  0.394, z = 1314.303
[jh_simplexP]     Iteracio   11 : iout =   8, q =   1, B(p) =   3, theta*=  2.818, z =  894.204
[jh_simplexP]     Iteracio   12 : iout =  10, q =   7, B(p) =   2, theta*=  1.790, z =  888.772
[jh_simplexP]     Iteracio   13 : iout =   8, q =   3, B(p) =   1, theta*=  1.718, z =  411.338
[jh_simplexP]     Iteracio   14 : iout =   2, q =   9, B(p) =  22, theta*=  0.237, z =  392.790
[jh_simplexP]     Iteracio   15 : iout =   3, q =  11, B(p) =  23, theta*=  0.050, z =  313.545
[jh_simplexP]     Iteracio   16 : iout =   2, q =   1, B(p) =   9, theta*=  0.049, z =  282.348
[jh_simplexP]     Iteracio   17 : iout =   6, q =   2, B(p) =   4, theta*=  1.079, z =  278.716
[jh_simplexP]     Iteracio   18 : iout =   7, q =  13, B(p) =  27, theta*=  1.030, z =  136.476
[jh_simplexP]     Iteracio   19 : iout =   -, q =   -, B(p) =   -, theta*=     -, z =    0.000
[jh_simplexP]     Solucio basica factible trobada, iteracio 19
[jh_simplexP]   Fase II
[jh_simplexP]     Iteracio   20 : iout =  10, q =   9, B(p) =   7, theta*=  2.023, z =  417.139
[jh_simplexP]     Iteracio   21 : iout =   1, q =  12, B(p) =   6, theta*=  2.362, z =  267.479
[jh_simplexP]     Iteracio   22 : iout =   8, q =   7, B(p) =   3, theta*=  0.057, z =  135.921
[jh_simplexP]     Iteracio   23 : iout =   5, q =  14, B(p) =   5, theta*=  0.100, z =  119.775
[jh_simplexP]     Iteracio   24 : iout =   8, q =   3, B(p) =   7, theta*=  0.520, z =  108.443
[jh_simplexP]     Iteracio   25 : iout =   9, q =   4, B(p) =   8, theta*=  0.316, z =   82.939
[jh_simplexP]     Iteracio   26 : iout =   9, q =   6, B(p) =   4, theta*=  0.180, z =   56.912
[jh_simplexP]     Iteracio   27 : iout =   9, q =  15, B(p) =   6, theta*= 105.586, z =   54.896
[jh_simplexP]     Iteracio   28 : iout =   4, q =  16, B(p) =  10, theta*= 73.570, z =   53.658
[jh_simplexP]     Iteracio   29 : iout =  10, q =   4, B(p) =   9, theta*=  1.514, z =   18.470
[jh_simplexP]     Iteracio   30 : iout =  10, q =   6, B(p) =   4, theta*=  0.622, z =  -23.129
[jh_simplexP]     Iteracio   31 : iout =   9, q =  10, B(p) =  15, theta*=  0.926, z =  -32.778
[jh_simplexP]     Iteracio   32 : iout =   9, q =  17, B(p) =  10, theta*= 150.172, z =  -42.298
[jh_simplexP]     Iteracio   33 : iout =   8, q =  19, B(p) =   3, theta*= 444.748, z =  -93.057
[jh_simplexP]     Iteracio   34 : iout =   -, q =   -, B(p) =   -, theta*=     -, z = -250.134
[jh_simplexP]     Solucio optima trobada, iteracio  34, z = -250.133658
[jh_simplexP] Fi simplex primal 

Solucio optima:

vb = 12  1 11 16 14  2 13 19 17  6
xb =    1.4    2.4    0.3  294.0    4.3    0.1    3.0  444.7  520.0    2.7
z = -250.1337
r =    0.0   -0.0   90.3   16.0   32.7    0.0  189.2   91.8   48.4   57.6    0.0    0.0   -0.0   -0.0    0.2    0.0    0.0    0.5    0.0    0.4
```
### Problema 3 conjunt 34

```
c=
  29  80 -61 -85  11  78  -1  60  67  97  46 -62 -94  24   0   0   0   0   0   0 

A=
  28 -60   0  -7  -3 -15 -69  -3 -22  72  57  89  92  28   0   0   0   0   0   0 
 -94  92  67 -97 -94  10  83   2  46 -61 -89  56  73  43   0   0   0   0   0   0 
  30  83  18  99 -24   7  28 -99 -57  38  21  29  78 -11   0   0   0   0   0   0 
 -24 -99  -5  51  41 -61 -55  50 -91  84  23  43  67   8   0   0   0   0   0   0 
 -58  31 -52 -58  75 -76  69 -31 -10  31  81  38  40  36  -1   0   0   0   0   0 
 -37 -64  44 -47  36  46  83 -94  40  -3 -92  37  96 -21   0  -1   0   0   0   0 
  31 -20  80  14  93  99 -54   7  80 -69  78 -57  77  -9   0   0   1   0   0   0 
  80  64  96  74  34  37  34   6  89  96  96   7  60   6   0   0   0   1   0   0 
  55  34  36  78  94   2  60  99   9  76   7  23  84   6   0   0   0   0   1   0 
  25  33  19  91  15  47  99  15  94  99  78  18  54  26   0   0   0   0   0   1 

b=
 187  37 240  32 115  23 351  29  97  10 
```
#### Output

```
[jh_simplexP] Inici simplex primal amb regla de Bland
[jh_simplexP]   Fase I
[jh_simplexP]     Iteracio    1 : iout =   8, q =   1, B(p) =  28, theta*=  0.362, z = 1121.000
[jh_simplexP]     Iteracio    2 : iout =  10, q =   2, B(p) =  30, theta*=  0.072, z = 1107.950
[jh_simplexP]     Iteracio    3 : iout =   8, q =   3, B(p) =   1, theta*=  0.162, z = 1103.248
[jh_simplexP]     Iteracio    4 : iout =  10, q =   4, B(p) =   2, theta*=  0.056, z = 1052.098
[jh_simplexP]     Iteracio    5 : iout =   6, q =   5, B(p) =  26, theta*=  0.487, z = 1031.453
[jh_simplexP]     Iteracio    6 : iout =  10, q =   1, B(p) =   4, theta*=  0.008, z =  951.743
[jh_simplexP]     Iteracio    7 : iout =  10, q =   2, B(p) =   1, theta*=  0.006, z =  950.516
[jh_simplexP]     Iteracio    8 : iout =  10, q =   8, B(p) =   2, theta*=  0.006, z =  950.032
[jh_simplexP]     Iteracio    9 : iout =  10, q =  16, B(p) =   8, theta*=  0.809, z =  949.273
[jh_simplexP]     Iteracio   10 : iout =   7, q =  17, B(p) =  27, theta*= 293.523, z =  948.020
[jh_simplexP]     Iteracio   11 : iout =  10, q =  12, B(p) =  16, theta*=  0.112, z =  654.497
[jh_simplexP]     Iteracio   12 : iout =   8, q =  18, B(p) =   3, theta*= 11.785, z =  638.813
[jh_simplexP]     Iteracio   13 : iout =   9, q =  19, B(p) =  29, theta*= 48.817, z =  621.054
[jh_simplexP]     Iteracio   14 : iout =   6, q =  26, B(p) =   5, theta*=  2.444, z =  572.237
[jh_simplexP]     Iteracio   15 : iout =   -, q =   -, B(p) =   -, theta*=     -, z =  471.778
[jh_simplexP]     Problema infactible.
```

### Problema 4 conjunt 34

```
c=
 -89 -86  -1 -57 -99 -53 -51 -81  -8  -5  -5  -2 -60 -75   0   0   0   0   0   0   0   0   0   0 

A=
  34  39  23  84  91  84  95  67  52  57   7  86  86  16  -1   0   0   0   0   0   0   0   0   0 
  21  38  87   2  96  75  47  83  17  58   7  35  15 100   0  -1   0   0   0   0   0   0   0   0 
  25  61   5  61  80  51  47  30  97  63  11  63  79  15   0   0  -1   0   0   0   0   0   0   0 
  74  43  10   7  51  21  47  55  64  92  94  92  97  74   0   0   0  -1   0   0   0   0   0   0 
  21  72  63   4  65  83  16  83  79 100   6  16  43  62   0   0   0   0  -1   0   0   0   0   0 
  18  66  30  48  82  79  51  30  20  67  23  32  27  91   0   0   0   0   0  -1   0   0   0   0 
  68  68  60  60  31   2   4  71  47  66  11  95  58  31   0   0   0   0   0   0  -1   0   0   0 
  97  26  77  51  45  91  55  42  49  47  90  25   4  93   0   0   0   0   0   0   0  -1   0   0 
  17  41  69  46  72  94  45  53  17  19  10  58   6   9   0   0   0   0   0   0   0   0  -1   0 
  60  87  45  69  61  66  30  93  52  93  96  74  67  24   0   0   0   0   0   0   0   0   0  -1 

b=
 820 680 687 820 712 663 671 791 555 916 
```
#### Output

```
[jh_simplexP] Inici simplex primal amb regla de Bland
[jh_simplexP]   Fase I
[jh_simplexP]     Iteracio    1 : iout =   8, q =   1, B(p) =  32, theta*=  8.155, z = 7315.000
[jh_simplexP]     Iteracio    2 : iout =   7, q =   2, B(p) =  31, theta*=  2.340, z = 3767.732
[jh_simplexP]     Iteracio    3 : iout =   9, q =   3, B(p) =  33, theta*=  6.480, z = 2774.501
[jh_simplexP]     Iteracio    4 : iout =   2, q =   5, B(p) =  26, theta*=  0.965, z = 2305.679
[jh_simplexP]     Iteracio    5 : iout =   6, q =   4, B(p) =  30, theta*=  0.222, z = 1933.786
[jh_simplexP]     Iteracio    6 : iout =   9, q =   8, B(p) =   3, theta*=  0.099, z =  392.172
[jh_simplexP]     Iteracio    7 : iout =   6, q =   6, B(p) =   4, theta*=  0.225, z =  372.775
[jh_simplexP]     Iteracio    8 : iout =   6, q =   7, B(p) =   6, theta*=  1.002, z =  368.020
[jh_simplexP]     Iteracio    9 : iout =   3, q =   9, B(p) =  27, theta*=  0.214, z =  313.193
[jh_simplexP]     Iteracio   10 : iout =   1, q =  10, B(p) =  25, theta*=  0.312, z =  276.343
[jh_simplexP]     Iteracio   11 : iout =   5, q =   3, B(p) =  29, theta*=  1.031, z =  156.435
[jh_simplexP]     Iteracio   12 : iout =   5, q =   6, B(p) =   3, theta*=  0.486, z =  110.614
[jh_simplexP]     Iteracio   13 : iout =  10, q =  11, B(p) =  34, theta*=  0.251, z =   61.507
[jh_simplexP]     Iteracio   14 : iout =   5, q =   3, B(p) =   6, theta*=  0.593, z =   27.792
[jh_simplexP]     Iteracio   15 : iout =   4, q =  12, B(p) =  28, theta*=  0.423, z =   25.700
[jh_simplexP]     Iteracio   16 : iout =   -, q =   -, B(p) =   -, theta*=     -, z =    0.000
[jh_simplexP]     Solucio basica factible trobada, iteracio 16
[jh_simplexP]   Fase II
[jh_simplexP]     Iteracio   17 : iout =   5, q =   6, B(p) =   3, theta*=  0.687, z = -774.572
[jh_simplexP]     Iteracio   18 : iout =   3, q =  13, B(p) =   9, theta*=  0.924, z = -894.881
[jh_simplexP]     Iteracio   19 : iout =   1, q =  14, B(p) =  10, theta*=  0.423, z = -982.164
[jh_simplexP]     Iteracio   20 : iout =   4, q =  15, B(p) =  12, theta*=  2.057, z = -1072.642
[jh_simplexP]     Iteracio   21 : iout =   4, q =  16, B(p) =  15, theta*=  1.275, z = -1073.653
[jh_simplexP]     Iteracio   22 : iout =   6, q =  18, B(p) =   7, theta*= 18.095, z = -1074.433
[jh_simplexP]     Iteracio   23 : iout =   5, q =   9, B(p) =   6, theta*=  0.307, z = -1088.044
[jh_simplexP]     Iteracio   24 : iout =   7, q =  15, B(p) =   2, theta*= 45.483, z = -1103.173
[jh_simplexP]     Iteracio   25 : iout =   5, q =   6, B(p) =   9, theta*=  0.122, z = -1122.540
[jh_simplexP]     Iteracio   26 : iout =   5, q =  20, B(p) =   6, theta*=  5.499, z = -1123.208
[jh_simplexP]     Iteracio   27 : iout =   7, q =   2, B(p) =  15, theta*=  1.014, z = -1128.115
[jh_simplexP]     Iteracio   28 : iout =   9, q =  17, B(p) =   8, theta*= 91.735, z = -1132.582
[jh_simplexP]     Iteracio   29 : iout =   7, q =  15, B(p) =   2, theta*= 200.097, z = -1183.095
[jh_simplexP]     Iteracio   30 : iout =  10, q =   7, B(p) =  11, theta*=  3.975, z = -1229.467
[jh_simplexP]     Iteracio   31 : iout =   8, q =  19, B(p) =   1, theta*= 198.423, z = -1259.729
[jh_simplexP]     Iteracio   32 : iout =  10, q =  11, B(p) =   7, theta*=  0.361, z = -1335.323
[jh_simplexP]     Iteracio   33 : iout =  10, q =  21, B(p) =  11, theta*= 30.122, z = -1337.273
[jh_simplexP]     Iteracio   34 : iout =   3, q =  22, B(p) =  13, theta*= 1934.565, z = -1383.652
[jh_simplexP]     Iteracio   35 : iout =   2, q =   7, B(p) =   5, theta*=  6.267, z = -2468.023
[jh_simplexP]     Iteracio   36 : iout =   2, q =  24, B(p) =   7, theta*= 564.000, z = -2594.600
[jh_simplexP]     Iteracio   37 : iout =   3, q =  13, B(p) =  22, theta*= 85.241, z = -4625.000
[jh_simplexP]     Iteracio   38 : iout =   1, q =  23, B(p) =  14, theta*= 631.500, z = -5477.414
[jh_simplexP]     Iteracio   39 : problema no acotat detectat.
[jh_simplexP]     Problema no acotat.
```
### Problema 1 conjunt 37

```
c=
   53  -65   31   70   30   63   74  -29   -1  -24  -37   -6   94    4    0    0    0    0    0    0

A=
  -87  -65  -31   84   63   97   96  -79  -73  -93   35   93   72   83    0    0    0    0    0    0
   45   12   54   82   89   16 -100   83  -86  -38    2   88   24  -74    0    0    0    0    0    0
   84  -50  -67   62  -72   89   33  -74   55   31   33  -60   23  -42    0    0    0    0    0    0
   39   28   42   27   19  -36  -81   73  -13   79   99   76  -71   48    0    0    0    0    0    0
   91   52   51   86   62   69   92   78   70   67   76   75   72   50    1    0    0    0    0    0
  -92   87  -45   15  -24   62   13  -12   19  -75   90  -44  100   41    0   -1    0    0    0    0
   97  -48   62    0   73   57   23   59 -100   32  -42  -51  -63  -93    0    0    1    0    0    0
   64  -34   19   61  -44   25   -6  -10   33  -89   61   96    5  -94    0    0    0   -1    0    0
  -58   89    6  -72   96   30   40   91   -2  -97   38   91   22  -32    0    0    0    0   -1    0
   47  -22  -48   29   64   36   -7  -89  -91   80   28   73  -29  -20    0    0    0    0    0    1

b=
  195  197   45  329  992  134    7   86  241   52
```

#### Output

```
[jh_simplexP] Inici simplex primal amb regla de Bland
[jh_simplexP]   Fase I
[jh_simplexP]     Iteracio    1 : iout =   7, q =   1, B(p) =  27, theta*=  0.072, z = 2278.000
[jh_simplexP]     Iteracio    2 : iout =   6, q =   2, B(p) =  26, theta*=  3.391, z = 2261.402
[jh_simplexP]     Iteracio    3 : iout =   3, q =   4, B(p) =  23, theta*=  1.038, z = 1709.298
[jh_simplexP]     Iteracio    4 : iout =   2, q =   3, B(p) =  22, theta*=  0.038, z = 1382.148
[jh_simplexP]     Iteracio    5 : iout =   2, q =   5, B(p) =   3, theta*=  0.035, z = 1366.560
[jh_simplexP]     Iteracio    6 : iout =  10, q =   6, B(p) =  30, theta*=  0.166, z = 1361.862
[jh_simplexP]     Iteracio    7 : iout =   7, q =   3, B(p) =   1, theta*=  0.183, z = 1325.968
[jh_simplexP]     Iteracio    8 : iout =   1, q =   7, B(p) =  21, theta*=  0.663, z = 1291.896
[jh_simplexP]     Iteracio    9 : iout =  10, q =   1, B(p) =   6, theta*=  0.012, z = 1106.070
[jh_simplexP]     Iteracio   10 : iout =   9, q =   9, B(p) =  29, theta*=  1.593, z = 1099.446
[jh_simplexP]     Iteracio   11 : iout =   5, q =  10, B(p) =  25, theta*=  0.773, z =  584.210
[jh_simplexP]     Iteracio   12 : iout =   8, q =  11, B(p) =  28, theta*=  1.221, z =  277.690
[jh_simplexP]     Iteracio   13 : iout =  10, q =   6, B(p) =   1, theta*=  0.127, z =   65.586
[jh_simplexP]     Iteracio   14 : iout =  10, q =   8, B(p) =   6, theta*=  0.188, z =   46.863
[jh_simplexP]     Iteracio   15 : iout =   7, q =  12, B(p) =   3, theta*=  0.864, z =   32.430
[jh_simplexP]     Iteracio   16 : iout =   4, q =  13, B(p) =  24, theta*=  0.640, z =   25.819
[jh_simplexP]     Iteracio   17 : iout =   -, q =   -, B(p) =   -, theta*=     -, z =    0.000
[jh_simplexP]     Solucio basica factible trobada, iteracio 17
[jh_simplexP]   Fase II
[jh_simplexP]     Iteracio   18 : iout =   4, q =  14, B(p) =  13, theta*=  0.065, z =  183.243
[jh_simplexP]     Iteracio   19 : iout =   2, q =   6, B(p) =   5, theta*=  0.572, z =  163.725
[jh_simplexP]     Iteracio   20 : iout =   1, q =  15, B(p) =   7, theta*= 374.780, z =  145.879
[jh_simplexP]     Iteracio   21 : iout =   9, q =  16, B(p) =   9, theta*= 111.238, z =   51.560
[jh_simplexP]     Iteracio   22 : iout =   7, q =   5, B(p) =  12, theta*=  0.618, z =    7.871
[jh_simplexP]     Iteracio   23 : iout =   5, q =  17, B(p) =  10, theta*= 100.232, z =   -5.342
[jh_simplexP]     Iteracio   24 : iout =  10, q =   7, B(p) =   8, theta*=  0.628, z =  -30.623
[jh_simplexP]     Iteracio   25 : iout =   7, q =  12, B(p) =   5, theta*=  0.181, z =  -46.593
[jh_simplexP]     Iteracio   26 : iout =  10, q =   9, B(p) =   7, theta*=  0.214, z =  -57.585
[jh_simplexP]     Iteracio   27 : iout =   1, q =  19, B(p) =  15, theta*= 529.954, z =  -60.506
[jh_simplexP]     Iteracio   28 : iout =   3, q =   8, B(p) =   4, theta*=  0.465, z = -172.346
[jh_simplexP]     Iteracio   29 : iout =   3, q =  13, B(p) =   8, theta*=  1.365, z = -193.197
[jh_simplexP]     Iteracio   30 : iout =   -, q =   -, B(p) =   -, theta*=     -, z = -197.476
[jh_simplexP]     Solucio optima trobada, iteracio  30, z = -197.476081
[jh_simplexP] Fi simplex primal 

Solucio optima:

vb = 19  6 13 14 17  2 12 11 16  9
xb =  651.2    3.4    1.4    0.4  433.8    6.6    0.8    2.9 1038.4    0.2
z = -197.4761
r =  152.8    0.0   43.4   25.9   52.8    0.0    2.4    9.2   -0.0  166.9    0.0    0.0    0.0    0.0    0.2    0.0    0.0    0.2    0.0    0.7
```
### Problema 2 conjunt 37

```
c=
  -37   78   98    8   56   36   11  -99  -31   -8   -9  -71   33   74    0    0    0    0    0    0

A=
    6   -9    4   80  -24   49  -83  -54   90   -2  -99   26   19   38    0    0    0    0    0    0
   69   38   39  -76    8  -48  -27  -24  -90   25  -15   69   21   37    0    0    0    0    0    0
  -18   61  -48   98   -4   68   -4   57  -91   83  -34  -22    0  -56    0    0    0    0    0    0
   35   80  -92  -49   23  100   -6  -10  -74  -38   90   72   19   52    0    0    0    0    0    0
   92  100   53   54   86   57   50   93   51   91   71   82   98   93    1    0    0    0    0    0
  -67   98  -15   80   46   42   34  -77   13 -100   29   82  -83   32    0   -1    0    0    0    0
   52  -87  -28   20   65  -17  -36  -32   15   70   -4   29  -93   95    0    0    1    0    0    0
   98   36  -15   83  -69   35   35   63  -23   62  -57  -59   47   15    0    0    0   -1    0    0
   14   93   94  -29  -71   54   96   23  -43   76  -47   94   49  -74    0    0    0    0    1    0
   56   15   50   37   39  -78   72   -4   36   34   24   68   47   56    0    0    0    0    0   -1

b=
   41   26   90  202 1072  113   50  250  330  451
```

#### Output

```
[jh_simplexP] Inici simplex primal amb regla de Bland
[jh_simplexP]   Fase I
[jh_simplexP]     Iteracio    1 : iout =   2, q =   1, B(p) =  22, theta*=  0.377, z = 2625.000
[jh_simplexP]     Iteracio    2 : iout =   2, q =   2, B(p) =   1, theta*=  0.684, z = 2498.014
[jh_simplexP]     Iteracio    3 : iout =   6, q =   4, B(p) =  26, theta*=  0.166, z = 2334.211
[jh_simplexP]     Iteracio    4 : iout =   3, q =   1, B(p) =  23, theta*=  0.175, z = 2143.096
[jh_simplexP]     Iteracio    5 : iout =   1, q =   3, B(p) =  21, theta*=  0.384, z = 2040.842
[jh_simplexP]     Iteracio    6 : iout =   7, q =   5, B(p) =  27, theta*=  0.864, z = 1945.333
[jh_simplexP]     Iteracio    7 : iout =   4, q =   6, B(p) =  24, theta*=  1.436, z = 1801.254
[jh_simplexP]     Iteracio    8 : iout =   9, q =   7, B(p) =  29, theta*=  1.060, z = 1641.509
[jh_simplexP]     Iteracio    9 : iout =   8, q =   9, B(p) =  28, theta*=  0.162, z =  855.444
[jh_simplexP]     Iteracio   10 : iout =   5, q =   8, B(p) =  25, theta*=  1.387, z =  717.123
[jh_simplexP]     Iteracio   11 : iout =   6, q =  10, B(p) =   4, theta*=  0.473, z =  184.750
[jh_simplexP]     Iteracio   12 : iout =   5, q =  11, B(p) =   8, theta*=  1.371, z =  156.310
[jh_simplexP]     Iteracio   13 : iout =   1, q =  12, B(p) =   3, theta*=  0.069, z =   49.013
[jh_simplexP]     Iteracio   14 : iout =   4, q =   4, B(p) =   6, theta*=  0.256, z =   39.014
[jh_simplexP]     Iteracio   15 : iout =   9, q =   3, B(p) =   7, theta*=  0.230, z =   18.269
[jh_simplexP]     Iteracio   16 : iout =  10, q =   8, B(p) =  30, theta*=  1.620, z =   13.140
[jh_simplexP]     Iteracio   17 : iout =   -, q =   -, B(p) =   -, theta*=     -, z =    0.000
[jh_simplexP]     Solucio basica factible trobada, iteracio 17
[jh_simplexP]   Fase II
[jh_simplexP]     Iteracio   18 : iout =   6, q =   6, B(p) =  10, theta*=  0.206, z = -246.297
[jh_simplexP]     Iteracio   19 : iout =   9, q =   7, B(p) =   3, theta*=  0.339, z = -519.979
[jh_simplexP]     Iteracio   20 : iout =   6, q =  14, B(p) =   6, theta*=  0.357, z = -537.351
[jh_simplexP]     Iteracio   21 : iout =   7, q =   3, B(p) =   5, theta*=  0.181, z = -557.947
[jh_simplexP]     Iteracio   22 : iout =   6, q =  15, B(p) =  14, theta*= 47.146, z = -562.608
[jh_simplexP]     Iteracio   23 : iout =   6, q =  20, B(p) =  15, theta*= 14.253, z = -580.288
[jh_simplexP]     Iteracio   24 : iout =   -, q =   -, B(p) =   -, theta*=     -, z = -586.610
[jh_simplexP]     Solucio optima trobada, iteracio  24, z = -586.610343
[jh_simplexP] Fi simplex primal 

Solucio optima:

vb = 12  2  1  4 11 20  3  9  7  8
xb =    3.1    0.7    1.8    1.5    1.2   14.3    0.3    1.5    0.0    3.4
z = -586.6103
r =    0.0    0.0   -0.0   -0.0   40.5  304.0   -0.0    0.0   -0.0  248.2    0.0    0.0   83.2   44.7    0.1    1.5    1.3    0.8    1.3    0.0
```

### Problema 3 conjunt 37

```
c=
  74  20 -36 -92 -83  45 -54  68 -61 -28  19  -9 -28 -57   0   0   0   0   0   0 

A=
  75 -69  12 -31   5  56 -42 -59  36  91  24  62   1  12   0   0   0   0   0   0 
  20 -97  92  40  52  40 -65 -30 -82  70  35  64  29  -4   0   0   0   0   0   0 
 -50  42  90  59 -84  98 -22  80 -65  19 -40  97  71  77   0   0   0   0   0   0 
  54  83  68  58 -32  94  47  60  45  70 -15  96 -79 -85   0   0   0   0   0   0 
 -63  88 -74 -89  19 -96  95  37  70  40 -45 -67  14  89   1   0   0   0   0   0 
  32 -70  -8  -9 -35  15  92 -28 -11  31  97  89  19  93   0  -1   0   0   0   0 
  45 -25  23 -11  98 -96  46 -58  25 -28  80 -40  61  12   0   0   1   0   0   0 
  58  69  51  25   5  60  86  43  70  48  50  88  84  28   0   0   0   1   0   0 
  75  97  50  74  95  46  25  69  88  60  63  58  27  36   0   0   0   0   1   0 
  66  94  57  32  74  31  79  14  23  14   6 100  34   7   0   0   0   0   0   1 

b=
 173 164 372 464  19 306 133  34  66  70 
```

#### Output

```
[jh_simplexP] Inici simplex primal amb regla de Bland
[jh_simplexP]   Fase I
[jh_simplexP]     Iteracio    1 : iout =   8, q =   1, B(p) =  28, theta*=  0.586, z = 1801.000
[jh_simplexP]     Iteracio    2 : iout =   8, q =   3, B(p) =   1, theta*=  0.667, z = 1618.103
[jh_simplexP]     Iteracio    3 : iout =   9, q =   5, B(p) =  29, theta*=  0.363, z = 1560.333
[jh_simplexP]     Iteracio    4 : iout =   5, q =  10, B(p) =  25, theta*=  0.556, z = 1501.739
[jh_simplexP]     Iteracio    5 : iout =   5, q =  15, B(p) =  10, theta*= 58.814, z = 1472.841
[jh_simplexP]     Iteracio    6 : iout =   7, q =  17, B(p) =  27, theta*= 82.953, z = 1442.925
[jh_simplexP]     Iteracio    7 : iout =   9, q =   4, B(p) =   5, theta*=  0.660, z = 1359.973
[jh_simplexP]     Iteracio    8 : iout =   8, q =  10, B(p) =   3, theta*=  0.422, z = 1343.945
[jh_simplexP]     Iteracio    9 : iout =   5, q =  19, B(p) =  15, theta*= 19.867, z = 1342.596
[jh_simplexP]     Iteracio   10 : iout =   9, q =   3, B(p) =   4, theta*=  0.080, z = 1340.436
[jh_simplexP]     Iteracio   11 : iout =  10, q =  20, B(p) =  30, theta*= 56.709, z = 1340.237
[jh_simplexP]     Iteracio   12 : iout =   9, q =   4, B(p) =   3, theta*=  0.085, z = 1283.528
[jh_simplexP]     Iteracio   13 : iout =   -, q =   -, B(p) =   -, theta*=     -, z = 1282.453
[jh_simplexP]     Problema infactible.
```
### Problema 4 conjunt 37

```
c=
  -37   78   98    8   56   36   11  -99  -31   -8   -9  -71   33   74    0    0    0    0    0    0

A=
    6   -9    4   80  -24   49  -83  -54   90   -2  -99   26   19   38    0    0    0    0    0    0
   69   38   39  -76    8  -48  -27  -24  -90   25  -15   69   21   37    0    0    0    0    0    0
  -18   61  -48   98   -4   68   -4   57  -91   83  -34  -22    0  -56    0    0    0    0    0    0
   35   80  -92  -49   23  100   -6  -10  -74  -38   90   72   19   52    0    0    0    0    0    0
   92  100   53   54   86   57   50   93   51   91   71   82   98   93    1    0    0    0    0    0
  -67   98  -15   80   46   42   34  -77   13 -100   29   82  -83   32    0   -1    0    0    0    0
   52  -87  -28   20   65  -17  -36  -32   15   70   -4   29  -93   95    0    0    1    0    0    0
   98   36  -15   83  -69   35   35   63  -23   62  -57  -59   47   15    0    0    0   -1    0    0
   14   93   94  -29  -71   54   96   23  -43   76  -47   94   49  -74    0    0    0    0    1    0
   56   15   50   37   39  -78   72   -4   36   34   24   68   47   56    0    0    0    0    0   -1

b=
   41   26   90  202 1072  113   50  250  330  451
```

#### Output

```
[jh_simplexP] Inici simplex primal amb regla de Bland
[jh_simplexP]   Fase I
[jh_simplexP]     Iteracio    1 : iout =   8, q =   1, B(p) =  32, theta*=  8.866, z = 6828.000
[jh_simplexP]     Iteracio    2 : iout =   5, q =   2, B(p) =  29, theta*=  4.076, z = 2483.670
[jh_simplexP]     Iteracio    3 : iout =   5, q =   3, B(p) =   2, theta*=  1.467, z = 2397.856
[jh_simplexP]     Iteracio    4 : iout =   2, q =   4, B(p) =  26, theta*=  2.368, z = 2050.292
[jh_simplexP]     Iteracio    5 : iout =   4, q =   5, B(p) =  28, theta*=  3.489, z = 1347.373
[jh_simplexP]     Iteracio    6 : iout =   2, q =   2, B(p) =   4, theta*=  0.034, z =  875.187
[jh_simplexP]     Iteracio    7 : iout =   2, q =   7, B(p) =   2, theta*=  0.102, z =  813.890
[jh_simplexP]     Iteracio    8 : iout =   2, q =   8, B(p) =   7, theta*=  0.162, z =  813.641
[jh_simplexP]     Iteracio    9 : iout =   6, q =   9, B(p) =  30, theta*=  0.954, z =  783.651
[jh_simplexP]     Iteracio   10 : iout =  10, q =   4, B(p) =  34, theta*=  0.052, z =  506.720
[jh_simplexP]     Iteracio   11 : iout =  10, q =   6, B(p) =   4, theta*=  0.078, z =  489.532
[jh_simplexP]     Iteracio   12 : iout =  10, q =   7, B(p) =   6, theta*=  0.165, z =  485.329
[jh_simplexP]     Iteracio   13 : iout =   2, q =  10, B(p) =   8, theta*=  0.251, z =  483.995
[jh_simplexP]     Iteracio   14 : iout =   2, q =  11, B(p) =  10, theta*=  0.158, z =  481.791
[jh_simplexP]     Iteracio   15 : iout =  10, q =   4, B(p) =   7, theta*=  1.330, z =  451.599
[jh_simplexP]     Iteracio   16 : iout =   1, q =   6, B(p) =  25, theta*=  1.977, z =  412.496
[jh_simplexP]     Iteracio   17 : iout =   7, q =  12, B(p) =  31, theta*=  0.711, z =  248.833
[jh_simplexP]     Iteracio   18 : iout =   5, q =  10, B(p) =   3, theta*=  0.782, z =  204.689
[jh_simplexP]     Iteracio   19 : iout =   9, q =  13, B(p) =  33, theta*=  0.592, z =  187.369
[jh_simplexP]     Iteracio   20 : iout =   3, q =   2, B(p) =  27, theta*=  0.071, z =    7.419
[jh_simplexP]     Iteracio   21 : iout =   -, q =   -, B(p) =   -, theta*=     -, z =    0.000
[jh_simplexP]     Solucio basica factible trobada, iteracio 21
[jh_simplexP]   Fase II
[jh_simplexP]     Iteracio   22 : iout =   3, q =   7, B(p) =   2, theta*=  0.102, z = -788.970
[jh_simplexP]     Iteracio   23 : iout =   8, q =   8, B(p) =   1, theta*=  0.999, z = -791.546
[jh_simplexP]     Iteracio   24 : iout =   7, q =  14, B(p) =  12, theta*=  0.473, z = -876.151
[jh_simplexP]     Iteracio   25 : iout =  10, q =   2, B(p) =   4, theta*=  0.153, z = -908.193
[jh_simplexP]     Iteracio   26 : iout =  10, q =  17, B(p) =   2, theta*= 110.029, z = -915.049
[jh_simplexP]     Iteracio   27 : iout =   4, q =  12, B(p) =   5, theta*=  0.045, z = -934.298
[jh_simplexP]     Iteracio   28 : iout =   2, q =  18, B(p) =  11, theta*= 83.214, z = -941.734
[jh_simplexP]     Iteracio   29 : iout =   7, q =   4, B(p) =  14, theta*=  0.084, z = -975.353
[jh_simplexP]     Iteracio   30 : iout =   3, q =  15, B(p) =   7, theta*=  7.563, z = -979.763
[jh_simplexP]     Iteracio   31 : iout =   7, q =  19, B(p) =   4, theta*= 18.090, z = -981.786
[jh_simplexP]     Iteracio   32 : iout =   4, q =  11, B(p) =  12, theta*=  0.546, z = -993.564
[jh_simplexP]     Iteracio   33 : iout =   4, q =  16, B(p) =  11, theta*= 75.522, z = -1017.486
[jh_simplexP]     Iteracio   34 : iout =  10, q =   4, B(p) =  17, theta*=  0.182, z = -1025.912
[jh_simplexP]     Iteracio   35 : iout =   9, q =  20, B(p) =  13, theta*= 75.212, z = -1025.990
[jh_simplexP]     Iteracio   36 : iout =  10, q =   2, B(p) =   4, theta*=  1.770, z = -1099.114
[jh_simplexP]     Iteracio   37 : iout =   8, q =  17, B(p) =   8, theta*= 234.920, z = -1136.402
[jh_simplexP]     Iteracio   38 : iout =  10, q =  21, B(p) =   2, theta*= 46.640, z = -1290.780
[jh_simplexP]     Iteracio   39 : iout =   6, q =  11, B(p) =   9, theta*=  3.879, z = -1341.322
[jh_simplexP]     Iteracio   40 : iout =   6, q =  22, B(p) =  11, theta*= 1170.705, z = -1410.919
[jh_simplexP]     Iteracio   41 : iout =   5, q =   2, B(p) =  10, theta*= 16.187, z = -2293.982
[jh_simplexP]     Iteracio   42 : iout =   9, q =  13, B(p) =  20, theta*=  3.701, z = -2384.052
[jh_simplexP]     Iteracio   43 : iout =   5, q =  23, B(p) =   2, theta*= 167.687, z = -2725.429
[jh_simplexP]     Iteracio   44 : iout =   7, q =  20, B(p) =  19, theta*= 725.863, z = -2962.109
[jh_simplexP]     Iteracio   45 : iout =   1, q =   2, B(p) =   6, theta*=  1.983, z = -4267.537
[jh_simplexP]     Iteracio   46 : iout =   1, q =  24, B(p) =   2, theta*= 132.077, z = -4486.508
[jh_simplexP]     Iteracio   47 : iout =   9, q =   4, B(p) =  13, theta*= 799.000, z = -5777.385
[jh_simplexP]     Iteracio   48 : problema no acotat detectat.
[jh_simplexP]     Problema no acotat.
```