
## 1. Descripció resumida de la implementació

Hem implementat la solució en Python a través de la classe `SimplexSolver`, aplicant el **mètode del símplex primal estructurat en dues fases** i utilitzant la **regla de Bland** en totes les decisions per evitar la degeneració.

Respecte l'aplicació del model, s'han utilitzat dos conjunts de dades: el conjunt 34 i el conjunt 37, resolent un total de 8 problemes d'optimització.
La seqüència d'execució del nostre codi és la següent:

1. **Preprocessament:** Forcem els termes independents a ser no negatius ($b_i \geq 0$) canviant el signe de tota la fila corresponent de la matriu $A$ si és necessari. Ademés, comprovem si la matriu conté una base trivial (identitat), cosa que ens permet passar a la fase II directament.

2. **Fase I (Solució Bàsica Factible Inicial):** Per assegurar l'arrencada independentment de les restriccions del problema, estenem la matriu afegint variables artificials (una matriu identitat) que componen la base inicial. Executem el nucli del símplex amb una funció objectiu que minimitza la suma d'aquestes variables artificials. Si el valor òptim d'aquesta fase és major que $0$ (considerant una tolerància de $10^{-8}$), declarem el problema com a infactible. En cas de degeneració (variables artificials encara a la base amb valor zero), intentem pivotar-les fora abans de passar a la fase II.

3. **Fase II (Optimització real):** Un cop garantida la factibilitat, eliminem les variables artificials no bàsiques i restablim el vector de costos original $c$. A partir d’aquí, continuem amb el símplex fins trobar una solució òptima o detectar que el problema és no acotat.

4. **Nucli iteratiu:** En cada iteració, calculem la solució bàsica $x_B = B^{-1} b$, els multiplicadors duals $w = c_B^T B^{-1}$ i la direcció $d_B = -B^{-1} A_q$.

5. **Pivoteig (Regla de Bland):** Per entrar a la base calculem els costos reduïts ($r_j$) manualment; d'entre els que són negatius (més enllà del llindar de tolerància), seleccionem el candidat amb l'índex $j$ més petit. Per sortir de la base apliquem el test del quocient $\min \{x_B/d\}$ considerant $d > 0$, escollint novament l'índex més baix en cas d'empat per complir amb la Regla de Bland.

<div style="page-break-before: always;"></div>

## 2. Solució obtinguda

A continuació es detalla el registre de l'execució pas a pas de cada iteració i les solucions òptimes per als 8 problemes dels conjunts indicats.

### Conjunt 34

#### Problema 1

```
[jh_simplexP]   Fase I
[jh_simplexP]     Iteració  1 : iout = 0, q =  1, B(p) = 27, theta*=  0.659, z = 2815.000
[jh_simplexP]     Iteració  2 : iout = 0, q =  2, B(p) = 26, theta*=  0.845, z = 2707.000
[jh_simplexP]     Iteració  3 : iout = 0, q =  3, B(p) = 23, theta*=  1.218, z = 2443.242
[jh_simplexP]     Iteració  4 : iout = 0, q =  4, B(p) = 29, theta*=  1.161, z = 1991.135
[jh_simplexP]     Iteració  5 : iout = 0, q =  5, B(p) =  2, theta*=  1.134, z = 1624.557
[jh_simplexP]     Iteració  6 : iout = 0, q =  6, B(p) =  3, theta*=  1.349, z = 1468.862
[jh_simplexP]     Iteració  7 : iout = 0, q =  7, B(p) = 28, theta*=  0.383, z = 1415.671
[jh_simplexP]     Iteració  8 : iout = 0, q =  3, B(p) = 30, theta*=  2.740, z = 1222.967
[jh_simplexP]     Iteració  9 : iout = 0, q =  8, B(p) = 25, theta*=  1.514, z =  957.663
[jh_simplexP]     Iteració 10 : iout = 0, q =  9, B(p) =  6, theta*=  2.064, z =  735.863
[jh_simplexP]     Iteració 11 : iout = 0, q = 10, B(p) = 24, theta*=  1.571, z =  682.006
[jh_simplexP]     Iteració 12 : iout = 0, q =  2, B(p) = 21, theta*=  0.057, z =  474.694
[jh_simplexP]     Iteració 13 : iout = 0, q =  6, B(p) =  2, theta*=  0.150, z =  463.824
[jh_simplexP]     Iteració 14 : iout = 0, q = 12, B(p) =  7, theta*=  0.717, z =  463.043
[jh_simplexP]     Iteració 15 : iout = 0, q = 11, B(p) =  9, theta*=  0.654, z =  364.841
[jh_simplexP]     Iteració 16 : iout = 0, q =  2, B(p) =  1, theta*=  0.104, z =  187.990
[jh_simplexP]     Iteració 17 : iout = 0, q =  9, B(p) =  5, theta*=  1.050, z =  157.658
[jh_simplexP]     Iteració 18 : iout = 0, q = 13, B(p) = 22, theta*=  0.714, z =  156.481
[jh_simplexP]     Iteració 19 : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z =    0.000
[jh_simplexP]     Solució bàsica factible trobada, iteració  19
[jh_simplexP]   Fase II
[jh_simplexP]     Iteració 20 : iout = 0, q =  5, B(p) =  9, theta*=  1.706, z =  -77.097
[jh_simplexP]     Iteració 21 : iout = 0, q = 15, B(p) =  8, theta*= 179.166, z = -132.510
[jh_simplexP]     Iteració 22 : iout = 0, q = 16, B(p) =  4, theta*=  9.381, z = -245.214
[jh_simplexP]     Iteració 23 : iout = 0, q =  1, B(p) = 15, theta*=  1.397, z = -251.992
[jh_simplexP]     Iteració 24 : iout = 0, q =  8, B(p) =  5, theta*=  0.702, z = -296.160
[jh_simplexP]     Iteració 25 : iout = 0, q = 17, B(p) =  2, theta*= 273.616, z = -323.140
[jh_simplexP]     Iteració 26 : iout = 0, q = 19, B(p) =  6, theta*= 147.251, z = -363.300
[jh_simplexP]     Iteració 27 : iout = 0, q =  2, B(p) =  1, theta*=  0.761, z = -531.478
[jh_simplexP]     Iteració 28 : iout = 0, q =  5, B(p) = 17, theta*=  0.832, z = -580.361
[jh_simplexP]     Iteració 29 : iout = 0, q =  7, B(p) = 11, theta*=  1.360, z = -588.083
[jh_simplexP]     Iteració 30 : iout = 0, q = 17, B(p) =  5, theta*= 356.161, z = -630.546
[jh_simplexP]     Iteració 31 : iout = 0, q = 20, B(p) =  2, theta*= 104.747, z = -639.085
[jh_simplexP]     Iteració 32 : iout = 0, q =  5, B(p) = 13, theta*=  2.203, z = -695.736
[jh_simplexP]     Iteració 33 : iout = 0, q =  1, B(p) = 17, theta*=  0.004, z = -787.471
[jh_simplexP]     Iteració 34 : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z = -787.487
[jh_simplexP]     Solució òptima trobada, iteració  34, z = -787.486630
[jh_simplexP] Fi simplex primal

Solució òptima:

vb =  19   5   7  10  20   8   1  12  16   3
xb = 764.42  2.21  1.48  2.47 291.82  3.60  0.00  2.21 322.82  0.52
z  = -787.486630
r  = 31.36 109.99 192.25 176.13  5.52 114.31 186.01  0.37  0.01  0.53
```

#### Problema 2

```
[jh_simplexP]   Fase I
[jh_simplexP]     Iteració  1 : iout = 0, q =  1, B(p) = 21, theta*=  0.113, z = 2231.000
[jh_simplexP]     Iteració  2 : iout = 0, q =  3, B(p) = 28, theta*=  0.030, z = 2223.189
[jh_simplexP]     Iteració  3 : iout = 0, q =  4, B(p) = 26, theta*=  0.951, z = 2208.520
[jh_simplexP]     Iteració  4 : iout = 0, q =  2, B(p) = 30, theta*=  0.145, z = 1766.561
[jh_simplexP]     Iteració  5 : iout = 0, q =  5, B(p) = 25, theta*=  2.494, z = 1724.330
[jh_simplexP]     Iteració  6 : iout = 0, q =  6, B(p) =  1, theta*=  0.132, z = 1441.831
[jh_simplexP]     Iteració  7 : iout = 0, q =  7, B(p) =  6, theta*=  0.108, z = 1412.946
[jh_simplexP]     Iteració  8 : iout = 0, q =  8, B(p) = 29, theta*=  0.008, z = 1393.373
[jh_simplexP]     Iteració  9 : iout = 0, q =  6, B(p) =  7, theta*=  0.666, z = 1379.835
[jh_simplexP]     Iteració 10 : iout = 0, q = 10, B(p) = 24, theta*=  0.394, z = 1314.303
[jh_simplexP]     Iteració 11 : iout = 0, q =  1, B(p) =  3, theta*=  2.818, z =  894.204
[jh_simplexP]     Iteració 12 : iout = 0, q =  7, B(p) =  2, theta*=  1.790, z =  888.772
[jh_simplexP]     Iteració 13 : iout = 0, q =  3, B(p) =  1, theta*=  1.718, z =  411.338
[jh_simplexP]     Iteració 14 : iout = 0, q =  9, B(p) = 22, theta*=  0.237, z =  392.790
[jh_simplexP]     Iteració 15 : iout = 0, q = 11, B(p) = 23, theta*=  0.050, z =  313.545
[jh_simplexP]     Iteració 16 : iout = 0, q =  1, B(p) =  9, theta*=  0.049, z =  282.348
[jh_simplexP]     Iteració 17 : iout = 0, q =  2, B(p) =  4, theta*=  1.079, z =  278.716
[jh_simplexP]     Iteració 18 : iout = 0, q = 13, B(p) = 27, theta*=  1.030, z =  136.476
[jh_simplexP]     Iteració 19 : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z =    0.000
[jh_simplexP]     Solució bàsica factible trobada, iteració  19
[jh_simplexP]   Fase II
[jh_simplexP]     Iteració 20 : iout = 0, q =  9, B(p) =  7, theta*=  2.023, z =  417.139
[jh_simplexP]     Iteració 21 : iout = 0, q = 12, B(p) =  6, theta*=  2.362, z =  267.479
[jh_simplexP]     Iteració 22 : iout = 0, q =  7, B(p) =  3, theta*=  0.057, z =  135.921
[jh_simplexP]     Iteració 23 : iout = 0, q = 14, B(p) =  5, theta*=  0.100, z =  119.775
[jh_simplexP]     Iteració 24 : iout = 0, q =  3, B(p) =  7, theta*=  0.520, z =  108.443
[jh_simplexP]     Iteració 25 : iout = 0, q =  4, B(p) =  8, theta*=  0.316, z =   82.939
[jh_simplexP]     Iteració 26 : iout = 0, q =  6, B(p) =  4, theta*=  0.180, z =   56.912
[jh_simplexP]     Iteració 27 : iout = 0, q = 15, B(p) =  6, theta*= 105.586, z =   54.896
[jh_simplexP]     Iteració 28 : iout = 0, q = 16, B(p) = 10, theta*= 73.570, z =   53.658
[jh_simplexP]     Iteració 29 : iout = 0, q =  4, B(p) =  9, theta*=  1.514, z =   18.470
[jh_simplexP]     Iteració 30 : iout = 0, q =  6, B(p) =  4, theta*=  0.622, z =  -23.129
[jh_simplexP]     Iteració 31 : iout = 0, q = 10, B(p) = 15, theta*=  0.926, z =  -32.778
[jh_simplexP]     Iteració 32 : iout = 0, q = 17, B(p) = 10, theta*= 150.172, z =  -42.298
[jh_simplexP]     Iteració 33 : iout = 0, q = 19, B(p) =  3, theta*= 444.748, z =  -93.057
[jh_simplexP]     Iteració 34 : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z = -250.134
[jh_simplexP]     Solució òptima trobada, iteració  34, z = -250.133658
[jh_simplexP] Fi simplex primal

Solució òptima:

vb =  12   1  11  16  14   2  13  19  17   6
xb =  1.40  2.36  0.30 293.96  4.32  0.15  3.03 444.75 520.01  2.74
z  = -250.133658
r  = 90.27 16.02 32.73 189.23 91.76 48.41 57.56  0.24  0.51  0.41
```
#### Problema 3

```
[jh_simplexP]   Fase I
[jh_simplexP]     Iteració  1 : iout = 0, q =  1, B(p) = 28, theta*=  0.362, z = 1121.000
[jh_simplexP]     Iteració  2 : iout = 0, q =  2, B(p) = 30, theta*=  0.072, z = 1107.950
[jh_simplexP]     Iteració  3 : iout = 0, q =  3, B(p) =  1, theta*=  0.162, z = 1103.248
[jh_simplexP]     Iteració  4 : iout = 0, q =  4, B(p) =  2, theta*=  0.056, z = 1052.098
[jh_simplexP]     Iteració  5 : iout = 0, q =  5, B(p) = 26, theta*=  0.487, z = 1031.453
[jh_simplexP]     Iteració  6 : iout = 0, q =  1, B(p) =  4, theta*=  0.008, z =  951.743
[jh_simplexP]     Iteració  7 : iout = 0, q =  2, B(p) =  1, theta*=  0.006, z =  950.516
[jh_simplexP]     Iteració  8 : iout = 0, q =  8, B(p) =  2, theta*=  0.006, z =  950.032
[jh_simplexP]     Iteració  9 : iout = 0, q = 16, B(p) =  8, theta*=  0.809, z =  949.273
[jh_simplexP]     Iteració 10 : iout = 0, q = 17, B(p) = 27, theta*= 293.523, z =  948.020
[jh_simplexP]     Iteració 11 : iout = 0, q = 12, B(p) = 16, theta*=  0.112, z =  654.497
[jh_simplexP]     Iteració 12 : iout = 0, q = 18, B(p) =  3, theta*= 11.785, z =  638.813
[jh_simplexP]     Iteració 13 : iout = 0, q = 19, B(p) = 29, theta*= 48.817, z =  621.054
[jh_simplexP]     Iteració 14 : iout = 0, q = 26, B(p) =  5, theta*=  2.444, z =  572.237
[jh_simplexP]     Iteració 15 : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z =  471.778
[jh_simplexP]     Problema infactible.
[jh_simplexP] Fi simplex primal
```
#### Problema 4

```
[jh_simplexP]   Fase I
[jh_simplexP]     Iteració  1 : iout = 0, q =  1, B(p) = 32, theta*=  8.155, z = 7315.000
[jh_simplexP]     Iteració  2 : iout = 0, q =  2, B(p) = 31, theta*=  2.340, z = 3767.732
[jh_simplexP]     Iteració  3 : iout = 0, q =  3, B(p) = 33, theta*=  6.480, z = 2774.501
[jh_simplexP]     Iteració  4 : iout = 0, q =  5, B(p) = 26, theta*=  0.965, z = 2305.679
[jh_simplexP]     Iteració  5 : iout = 0, q =  4, B(p) = 30, theta*=  0.222, z = 1933.786
[jh_simplexP]     Iteració  6 : iout = 0, q =  8, B(p) =  3, theta*=  0.099, z =  392.172
[jh_simplexP]     Iteració  7 : iout = 0, q =  6, B(p) =  4, theta*=  0.225, z =  372.775
[jh_simplexP]     Iteració  8 : iout = 0, q =  7, B(p) =  6, theta*=  1.002, z =  368.020
[jh_simplexP]     Iteració  9 : iout = 0, q =  9, B(p) = 27, theta*=  0.214, z =  313.193
[jh_simplexP]     Iteració 10 : iout = 0, q = 10, B(p) = 25, theta*=  0.312, z =  276.343
[jh_simplexP]     Iteració 11 : iout = 0, q =  3, B(p) = 29, theta*=  1.031, z =  156.435
[jh_simplexP]     Iteració 12 : iout = 0, q =  6, B(p) =  3, theta*=  0.486, z =  110.614
[jh_simplexP]     Iteració 13 : iout = 0, q = 11, B(p) = 34, theta*=  0.251, z =   61.507
[jh_simplexP]     Iteració 14 : iout = 0, q =  3, B(p) =  6, theta*=  0.593, z =   27.792
[jh_simplexP]     Iteració 15 : iout = 0, q = 12, B(p) = 28, theta*=  0.423, z =   25.700
[jh_simplexP]     Iteració 16 : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z =    0.000
[jh_simplexP]     Solució bàsica factible trobada, iteració  16
[jh_simplexP]   Fase II
[jh_simplexP]     Iteració 17 : iout = 0, q =  6, B(p) =  3, theta*=  0.687, z = -774.572
[jh_simplexP]     Iteració 18 : iout = 0, q = 13, B(p) =  9, theta*=  0.924, z = -894.881
[jh_simplexP]     Iteració 19 : iout = 0, q = 14, B(p) = 10, theta*=  0.423, z = -982.164
[jh_simplexP]     Iteració 20 : iout = 0, q = 15, B(p) = 12, theta*=  2.057, z = -1072.642
[jh_simplexP]     Iteració 21 : iout = 0, q = 16, B(p) = 15, theta*=  1.275, z = -1073.653
[jh_simplexP]     Iteració 22 : iout = 0, q = 18, B(p) =  7, theta*= 18.095, z = -1074.433
[jh_simplexP]     Iteració 23 : iout = 0, q =  9, B(p) =  6, theta*=  0.307, z = -1088.044
[jh_simplexP]     Iteració 24 : iout = 0, q = 15, B(p) =  2, theta*= 45.483, z = -1103.173
[jh_simplexP]     Iteració 25 : iout = 0, q =  6, B(p) =  9, theta*=  0.122, z = -1122.540
[jh_simplexP]     Iteració 26 : iout = 0, q = 20, B(p) =  6, theta*=  5.499, z = -1123.208
[jh_simplexP]     Iteració 27 : iout = 0, q =  2, B(p) = 15, theta*=  1.014, z = -1128.115
[jh_simplexP]     Iteració 28 : iout = 0, q = 17, B(p) =  8, theta*= 91.735, z = -1132.582
[jh_simplexP]     Iteració 29 : iout = 0, q = 15, B(p) =  2, theta*= 200.097, z = -1183.095
[jh_simplexP]     Iteració 30 : iout = 0, q =  7, B(p) = 11, theta*=  3.975, z = -1229.467
[jh_simplexP]     Iteració 31 : iout = 0, q = 19, B(p) =  1, theta*= 198.423, z = -1259.729
[jh_simplexP]     Iteració 32 : iout = 0, q = 11, B(p) =  7, theta*=  0.361, z = -1335.323
[jh_simplexP]     Iteració 33 : iout = 0, q = 21, B(p) = 11, theta*= 30.122, z = -1337.273
[jh_simplexP]     Iteració 34 : iout = 0, q = 22, B(p) = 13, theta*= 1934.565, z = -1383.652
[jh_simplexP]     Iteració 35 : iout = 0, q =  7, B(p) =  5, theta*=  6.267, z = -2468.023
[jh_simplexP]     Iteració 36 : iout = 0, q = 24, B(p) =  7, theta*= 564.000, z = -2594.600
[jh_simplexP]     Iteració 37 : iout = 0, q = 13, B(p) = 22, theta*= 85.241, z = -4625.000
[jh_simplexP]     Iteració 38 : iout = 0, q = 23, B(p) = 14, theta*= 631.500, z = -5477.414
[jh_simplexP]     Problema no acotat.
[jh_simplexP] Fi simplex primal
```
### Conjunt 37

#### Problema 1

```
[jh_simplexP]   Fase I
[jh_simplexP]     Iteració  1 : iout = 0, q =  1, B(p) = 27, theta*=  0.072, z = 2278.000
[jh_simplexP]     Iteració  2 : iout = 0, q =  2, B(p) = 26, theta*=  3.391, z = 2261.402
[jh_simplexP]     Iteració  3 : iout = 0, q =  4, B(p) = 23, theta*=  1.038, z = 1709.298
[jh_simplexP]     Iteració  4 : iout = 0, q =  3, B(p) = 22, theta*=  0.038, z = 1382.148
[jh_simplexP]     Iteració  5 : iout = 0, q =  5, B(p) =  3, theta*=  0.035, z = 1366.560
[jh_simplexP]     Iteració  6 : iout = 0, q =  6, B(p) = 30, theta*=  0.166, z = 1361.862
[jh_simplexP]     Iteració  7 : iout = 0, q =  3, B(p) =  1, theta*=  0.183, z = 1325.968
[jh_simplexP]     Iteració  8 : iout = 0, q =  7, B(p) = 21, theta*=  0.663, z = 1291.896
[jh_simplexP]     Iteració  9 : iout = 0, q =  1, B(p) =  6, theta*=  0.012, z = 1106.070
[jh_simplexP]     Iteració 10 : iout = 0, q =  9, B(p) = 29, theta*=  1.593, z = 1099.446
[jh_simplexP]     Iteració 11 : iout = 0, q = 10, B(p) = 25, theta*=  0.773, z =  584.210
[jh_simplexP]     Iteració 12 : iout = 0, q = 11, B(p) = 28, theta*=  1.221, z =  277.690
[jh_simplexP]     Iteració 13 : iout = 0, q =  6, B(p) =  1, theta*=  0.127, z =   65.586
[jh_simplexP]     Iteració 14 : iout = 0, q =  8, B(p) =  6, theta*=  0.188, z =   46.863
[jh_simplexP]     Iteració 15 : iout = 0, q = 12, B(p) =  3, theta*=  0.864, z =   32.430
[jh_simplexP]     Iteració 16 : iout = 0, q = 13, B(p) = 24, theta*=  0.640, z =   25.819
[jh_simplexP]     Iteració 17 : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z =    0.000
[jh_simplexP]     Solució bàsica factible trobada, iteració  17
[jh_simplexP]   Fase II
[jh_simplexP]     Iteració 18 : iout = 0, q = 14, B(p) = 13, theta*=  0.065, z =  183.243
[jh_simplexP]     Iteració 19 : iout = 0, q =  6, B(p) =  5, theta*=  0.572, z =  163.725
[jh_simplexP]     Iteració 20 : iout = 0, q = 15, B(p) =  7, theta*= 374.780, z =  145.879
[jh_simplexP]     Iteració 21 : iout = 0, q = 16, B(p) =  9, theta*= 111.238, z =   51.560
[jh_simplexP]     Iteració 22 : iout = 0, q =  5, B(p) = 12, theta*=  0.618, z =    7.871
[jh_simplexP]     Iteració 23 : iout = 0, q = 17, B(p) = 10, theta*= 100.232, z =   -5.342
[jh_simplexP]     Iteració 24 : iout = 0, q =  7, B(p) =  8, theta*=  0.628, z =  -30.623
[jh_simplexP]     Iteració 25 : iout = 0, q = 12, B(p) =  5, theta*=  0.181, z =  -46.593
[jh_simplexP]     Iteració 26 : iout = 0, q =  9, B(p) =  7, theta*=  0.214, z =  -57.585
[jh_simplexP]     Iteració 27 : iout = 0, q = 19, B(p) = 15, theta*= 529.954, z =  -60.506
[jh_simplexP]     Iteració 28 : iout = 0, q =  8, B(p) =  4, theta*=  0.465, z = -172.346
[jh_simplexP]     Iteració 29 : iout = 0, q = 13, B(p) =  8, theta*=  1.365, z = -193.197
[jh_simplexP]     Iteració 30 : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z = -197.476
[jh_simplexP]     Solució òptima trobada, iteració  30, z = -197.476081
[jh_simplexP] Fi simplex primal

Solució òptima:

vb =  19   6  13  14  17   2  12  11  16   9
xb = 651.17  3.43  1.36  0.41 433.77  6.64  0.79  2.89 1038.43  0.18
z  = -197.476081
r  = 152.79 43.35 25.94 52.76  2.43  9.20 166.89  0.22  0.20  0.70
```
#### Problema 2

```
[jh_simplexP]   Fase I
[jh_simplexP]     Iteració  1 : iout = 0, q =  1, B(p) = 22, theta*=  0.377, z = 2625.000
[jh_simplexP]     Iteració  2 : iout = 0, q =  2, B(p) =  1, theta*=  0.684, z = 2498.014
[jh_simplexP]     Iteració  3 : iout = 0, q =  4, B(p) = 26, theta*=  0.166, z = 2334.211
[jh_simplexP]     Iteració  4 : iout = 0, q =  1, B(p) = 23, theta*=  0.175, z = 2143.096
[jh_simplexP]     Iteració  5 : iout = 0, q =  3, B(p) = 21, theta*=  0.384, z = 2040.842
[jh_simplexP]     Iteració  6 : iout = 0, q =  5, B(p) = 27, theta*=  0.864, z = 1945.333
[jh_simplexP]     Iteració  7 : iout = 0, q =  6, B(p) = 24, theta*=  1.436, z = 1801.254
[jh_simplexP]     Iteració  8 : iout = 0, q =  7, B(p) = 29, theta*=  1.060, z = 1641.509
[jh_simplexP]     Iteració  9 : iout = 0, q =  9, B(p) = 28, theta*=  0.162, z =  855.444
[jh_simplexP]     Iteració 10 : iout = 0, q =  8, B(p) = 25, theta*=  1.387, z =  717.123
[jh_simplexP]     Iteració 11 : iout = 0, q = 10, B(p) =  4, theta*=  0.473, z =  184.750
[jh_simplexP]     Iteració 12 : iout = 0, q = 11, B(p) =  8, theta*=  1.371, z =  156.310
[jh_simplexP]     Iteració 13 : iout = 0, q = 12, B(p) =  3, theta*=  0.069, z =   49.013
[jh_simplexP]     Iteració 14 : iout = 0, q =  4, B(p) =  6, theta*=  0.256, z =   39.014
[jh_simplexP]     Iteració 15 : iout = 0, q =  3, B(p) =  7, theta*=  0.230, z =   18.269
[jh_simplexP]     Iteració 16 : iout = 0, q =  8, B(p) = 30, theta*=  1.620, z =   13.140
[jh_simplexP]     Iteració 17 : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z =    0.000
[jh_simplexP]     Solució bàsica factible trobada, iteració  17
[jh_simplexP]   Fase II
[jh_simplexP]     Iteració 18 : iout = 0, q =  6, B(p) = 10, theta*=  0.206, z = -246.297
[jh_simplexP]     Iteració 19 : iout = 0, q =  7, B(p) =  3, theta*=  0.339, z = -519.979
[jh_simplexP]     Iteració 20 : iout = 0, q = 14, B(p) =  6, theta*=  0.357, z = -537.351
[jh_simplexP]     Iteració 21 : iout = 0, q =  3, B(p) =  5, theta*=  0.181, z = -557.947
[jh_simplexP]     Iteració 22 : iout = 0, q = 15, B(p) = 14, theta*= 47.146, z = -562.608
[jh_simplexP]     Iteració 23 : iout = 0, q = 20, B(p) = 15, theta*= 14.253, z = -580.288
[jh_simplexP]     Iteració 24 : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z = -586.610
[jh_simplexP]     Solució òptima trobada, iteració  24, z = -586.610343
[jh_simplexP] Fi simplex primal

Solució òptima:

vb =  12   2   1   4  11  20   3   9   7   8
xb =  3.12  0.72  1.79  1.50  1.19 14.25  0.31  1.54  0.01  3.43
z  = -586.610343
r  = 40.50 303.97 248.20 83.22 44.69  0.13  1.53  1.26  0.80  1.31
```

#### Problema 3

```
[jh_simplexP]   Fase I
[jh_simplexP]     Iteració  1 : iout = 0, q =  1, B(p) = 28, theta*=  0.586, z = 1801.000
[jh_simplexP]     Iteració  2 : iout = 0, q =  3, B(p) =  1, theta*=  0.667, z = 1618.103
[jh_simplexP]     Iteració  3 : iout = 0, q =  5, B(p) = 29, theta*=  0.363, z = 1560.333
[jh_simplexP]     Iteració  4 : iout = 0, q = 10, B(p) = 25, theta*=  0.556, z = 1501.739
[jh_simplexP]     Iteració  5 : iout = 0, q = 15, B(p) = 10, theta*= 58.814, z = 1472.841
[jh_simplexP]     Iteració  6 : iout = 0, q = 17, B(p) = 27, theta*= 82.953, z = 1442.925
[jh_simplexP]     Iteració  7 : iout = 0, q =  4, B(p) =  5, theta*=  0.660, z = 1359.973
[jh_simplexP]     Iteració  8 : iout = 0, q = 10, B(p) =  3, theta*=  0.422, z = 1343.945
[jh_simplexP]     Iteració  9 : iout = 0, q = 19, B(p) = 15, theta*= 19.867, z = 1342.596
[jh_simplexP]     Iteració 10 : iout = 0, q =  3, B(p) =  4, theta*=  0.080, z = 1340.436
[jh_simplexP]     Iteració 11 : iout = 0, q = 20, B(p) = 30, theta*= 56.709, z = 1340.237
[jh_simplexP]     Iteració 12 : iout = 0, q =  4, B(p) =  3, theta*=  0.085, z = 1283.528
[jh_simplexP]     Iteració 13 : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z = 1282.453
[jh_simplexP]     Problema infactible.
[jh_simplexP] Fi simplex primal
```
#### Problema 4

```
[jh_simplexP]   Fase I
[jh_simplexP]     Iteració  1 : iout = 0, q =  1, B(p) = 32, theta*=  8.866, z = 6828.000
[jh_simplexP]     Iteració  2 : iout = 0, q =  2, B(p) = 29, theta*=  4.076, z = 2483.670
[jh_simplexP]     Iteració  3 : iout = 0, q =  3, B(p) =  2, theta*=  1.467, z = 2397.856
[jh_simplexP]     Iteració  4 : iout = 0, q =  4, B(p) = 26, theta*=  2.368, z = 2050.292
[jh_simplexP]     Iteració  5 : iout = 0, q =  5, B(p) = 28, theta*=  3.489, z = 1347.373
[jh_simplexP]     Iteració  6 : iout = 0, q =  2, B(p) =  4, theta*=  0.034, z =  875.187
[jh_simplexP]     Iteració  7 : iout = 0, q =  7, B(p) =  2, theta*=  0.102, z =  813.890
[jh_simplexP]     Iteració  8 : iout = 0, q =  8, B(p) =  7, theta*=  0.162, z =  813.641
[jh_simplexP]     Iteració  9 : iout = 0, q =  9, B(p) = 30, theta*=  0.954, z =  783.651
[jh_simplexP]     Iteració 10 : iout = 0, q =  4, B(p) = 34, theta*=  0.052, z =  506.720
[jh_simplexP]     Iteració 11 : iout = 0, q =  6, B(p) =  4, theta*=  0.078, z =  489.532
[jh_simplexP]     Iteració 12 : iout = 0, q =  7, B(p) =  6, theta*=  0.165, z =  485.329
[jh_simplexP]     Iteració 13 : iout = 0, q = 10, B(p) =  8, theta*=  0.251, z =  483.995
[jh_simplexP]     Iteració 14 : iout = 0, q = 11, B(p) = 10, theta*=  0.158, z =  481.791
[jh_simplexP]     Iteració 15 : iout = 0, q =  4, B(p) =  7, theta*=  1.330, z =  451.599
[jh_simplexP]     Iteració 16 : iout = 0, q =  6, B(p) = 25, theta*=  1.977, z =  412.496
[jh_simplexP]     Iteració 17 : iout = 0, q = 12, B(p) = 31, theta*=  0.711, z =  248.833
[jh_simplexP]     Iteració 18 : iout = 0, q = 10, B(p) =  3, theta*=  0.782, z =  204.689
[jh_simplexP]     Iteració 19 : iout = 0, q = 13, B(p) = 33, theta*=  0.592, z =  187.369
[jh_simplexP]     Iteració 20 : iout = 0, q =  2, B(p) = 27, theta*=  0.071, z =    7.419
[jh_simplexP]     Iteració 21 : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z =    0.000
[jh_simplexP]     Solució bàsica factible trobada, iteració  21
[jh_simplexP]   Fase II
[jh_simplexP]     Iteració 22 : iout = 0, q =  7, B(p) =  2, theta*=  0.102, z = -788.970
[jh_simplexP]     Iteració 23 : iout = 0, q =  8, B(p) =  1, theta*=  0.999, z = -791.546
[jh_simplexP]     Iteració 24 : iout = 0, q = 14, B(p) = 12, theta*=  0.473, z = -876.151
[jh_simplexP]     Iteració 25 : iout = 0, q =  2, B(p) =  4, theta*=  0.153, z = -908.193
[jh_simplexP]     Iteració 26 : iout = 0, q = 17, B(p) =  2, theta*= 110.029, z = -915.049
[jh_simplexP]     Iteració 27 : iout = 0, q = 12, B(p) =  5, theta*=  0.045, z = -934.298
[jh_simplexP]     Iteració 28 : iout = 0, q = 18, B(p) = 11, theta*= 83.214, z = -941.734
[jh_simplexP]     Iteració 29 : iout = 0, q =  4, B(p) = 14, theta*=  0.084, z = -975.353
[jh_simplexP]     Iteració 30 : iout = 0, q = 15, B(p) =  7, theta*=  7.563, z = -979.763
[jh_simplexP]     Iteració 31 : iout = 0, q = 19, B(p) =  4, theta*= 18.090, z = -981.786
[jh_simplexP]     Iteració 32 : iout = 0, q = 11, B(p) = 12, theta*=  0.546, z = -993.564
[jh_simplexP]     Iteració 33 : iout = 0, q = 16, B(p) = 11, theta*= 75.522, z = -1017.486
[jh_simplexP]     Iteració 34 : iout = 0, q =  4, B(p) = 17, theta*=  0.182, z = -1025.912
[jh_simplexP]     Iteració 35 : iout = 0, q = 20, B(p) = 13, theta*= 75.212, z = -1025.990
[jh_simplexP]     Iteració 36 : iout = 0, q =  2, B(p) =  4, theta*=  1.770, z = -1099.114
[jh_simplexP]     Iteració 37 : iout = 0, q = 17, B(p) =  8, theta*= 234.920, z = -1136.402
[jh_simplexP]     Iteració 38 : iout = 0, q = 21, B(p) =  2, theta*= 46.640, z = -1290.780
[jh_simplexP]     Iteració 39 : iout = 0, q = 11, B(p) =  9, theta*=  3.879, z = -1341.322
[jh_simplexP]     Iteració 40 : iout = 0, q = 22, B(p) = 11, theta*= 1170.705, z = -1410.919
[jh_simplexP]     Iteració 41 : iout = 0, q =  2, B(p) = 10, theta*= 16.187, z = -2293.982
[jh_simplexP]     Iteració 42 : iout = 0, q = 13, B(p) = 20, theta*=  3.701, z = -2384.052
[jh_simplexP]     Iteració 43 : iout = 0, q = 23, B(p) =  2, theta*= 167.687, z = -2725.429
[jh_simplexP]     Iteració 44 : iout = 0, q = 20, B(p) = 19, theta*= 725.863, z = -2962.109
[jh_simplexP]     Iteració 45 : iout = 0, q =  2, B(p) =  6, theta*=  1.983, z = -4267.537
[jh_simplexP]     Iteració 46 : iout = 0, q = 24, B(p) =  2, theta*= 132.077, z = -4486.508
[jh_simplexP]     Iteració 47 : iout = 0, q =  4, B(p) = 13, theta*= 799.000, z = -5777.385
[jh_simplexP]     Problema no acotat.
[jh_simplexP] Fi simplex primal
```