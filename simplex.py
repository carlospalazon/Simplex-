class ResolvedorSimplex:
    def __init__(self, c, A, b):
        # Copia floats de c, A, b
        self.c = [float(val) for val in c]
        self.A = [[float(val) for val in fila] for fila in A]
        self.b = [float(val) for val in b]

        # Definim el nombre de restriccions
        self.m = len(self.A)
        # Definim el nombre de variables
        self.n = len(self.c)
        # Definim la tolerància, ja que els floats mai arriben a ser 0
        self.tol = 1e-8
        # Definim el número d'iteracions
        self.comptador_iter = 0

        # Garantim que b[i] sigui positiu. Si no ho és, multipliquem tota la fila per -1
        for i in range(self.m):
            if self.b[i] < 0:
                self.b[i] *= -1.0
                for j in range(self.n):
                    self.A[i][j] *= -1.0

        self.base_trivial = self._comprova_base_trivial()

    # Recorrem les últimes m columnes i comprovem que cada columna té exactament un 1 i m-1 zeros
    def _comprova_base_trivial(self):
        for j in range(self.n - self.m, self.n):
            col = [self.A[i][j] for i in range(self.m)]
            if sum(1 for x in col if abs(x - 1) < self.tol) != 1:
                return False
            if sum(1 for x in col if abs(x) < self.tol) != self.m - 1:
                return False
        return True

    # Multipliquem una matriu M per un vector v
    def _mat_vec(self, M, v):
        return [sum(M[i][k] * v[k] for k in range(self.m)) for i in range(self.m)]

    # Obtenim les variables bàsiques
    def _obte_x_B(self):
        return self._mat_vec(self.B_inv, self.b)

    # Calculem el vector w = c_B^T B^-1
    def _obte_w(self, c_B):
        return [sum(c_B[i] * self.B_inv[i][j] for i in range(self.m))
                for j in range(self.m)]

    # Calculem la direcció
    def _obte_d(self, col_A_q):
        Binv_Aq = self._mat_vec(self.B_inv, col_A_q)
        return [-val for val in Binv_Aq]

    # Actualitzem la inversa a través de la matriu eta
    def _actualitza_B_inv(self, fila_p, d_B):
        d_pq = d_B[fila_p]
        nova_B_inv = []
        for i in range(self.m):
            if i == fila_p:
                eta_p = -1.0 / d_pq
                nova_fila = [eta_p * self.B_inv[fila_p][j] for j in range(self.m)]
            else:
                eta_i = -d_B[i] / d_pq
                nova_fila = [self.B_inv[i][j] + eta_i * self.B_inv[fila_p][j]
                             for j in range(self.m)]
            nova_B_inv.append(nova_fila)
        self.B_inv = nova_B_inv


    def resol(self):

        # Si no tenim una matriu identitat, pasem a la fase I
        if not self.base_trivial:
            print("[jh_simplexP]   Fase I")
            estat = self._fase_u()
            # Si el problema és infactible, aturem
            if estat == "Infactible":
                print("[jh_simplexP]     Problema infactible.")
                print("[jh_simplexP] Fi simplex primal")
                return
            print(f"[jh_simplexP]     Solució bàsica factible trobada, iteració {self.comptador_iter:3d}")
        else:
            # Si hi ha una identitat, pasem a la fase II
            print("[jh_simplexP] Base trivial detectada, saltant Fase I")
            self.B_inv = [[1.0 if i == j else 0.0 for j in range(self.m)]
                          for i in range(self.m)]
            self.idx_B = list(range(self.n - self.m, self.n))
            self.idx_N = list(range(self.n - self.m))
            self.A_ext = [fila[:] for fila in self.A]

        print("[jh_simplexP]   Fase II")
        estat = self._fase_dos()
        if estat == "No acotat":
            print("[jh_simplexP]     Problema no acotat.")
            print("[jh_simplexP] Fi simplex primal")
        elif estat == "Òptim":
            self._imprimeix_solucio_final()


    def _fase_u(self):
        # Ampliem A amb variables artificials
        self.A_ext = []
        for i in range(self.m):
            fila_ext = self.A[i][:] + [1.0 if i == j else 0.0 for j in range(self.m)]
            self.A_ext.append(fila_ext)

        # Definim els nous costos amb les variables artificials
        self.c_ext = [0.0] * self.n + [1.0] * self.m
        # Definim els índexs de la nova base
        self.idx_B = list(range(self.n, self.n + self.m))
        # Definim els índexs de les variables no bàsiques
        self.idx_N = list(range(self.n))
        # B^-1 inicial = I
        self.B_inv = [[1.0 if i == j else 0.0 for j in range(self.m)]
                      for i in range(self.m)]

        self._nucli_simplex(self.A_ext, self.c_ext, fase=1)

        c_B = [self.c_ext[j] for j in self.idx_B]
        x_B = self._obte_x_B()
        z = sum(c_B[i] * x_B[i] for i in range(self.m))

        # z*_I > 0 → infactible (Imatge 4)
        if z > self.tol:
            return "Infactible"

        # z*_I = 0 però hi ha variables artificials a la base (cas degenerat)
        # → intentar pivotar-les fora (Imatge 4)
        for fila_p, idx_art in [(i, idx) for i, idx in enumerate(self.idx_B) if idx >= self.n]:
            for j in sorted(self.idx_N):
                if j >= self.n:
                    continue
                col_A = [self.A_ext[i][j] for i in range(self.m)]
                d_B = self._obte_d(col_A)
                if abs(d_B[fila_p]) > self.tol:
                    self._actualitza_B_inv(fila_p, d_B)
                    idx_q_en_N = self.idx_N.index(j)
                    self.idx_B[fila_p] = j
                    self.idx_N[idx_q_en_N] = idx_art
                    break

        return "Factible"


    def _fase_dos(self):
        # Eliminem variables artificials no bàsiques
        self.idx_N = [idx for idx in self.idx_N if idx < self.n]
        # Cost 0 per a les artificials que poguessin quedar a la base (degenerat)
        self.c_fase2 = self.c[:] + [0.0] * self.m
        return self._nucli_simplex(self.A_ext, self.c_fase2, fase=2)


    def _nucli_simplex(self, mat_A, vec_c, fase=1):
        while True:
            self.comptador_iter += 1
            c_B = [vec_c[j] for j in self.idx_B]
            x_B = self._obte_x_B()
            z_actual = sum(c_B[i] * x_B[i] for i in range(self.m))

            w = self._obte_w(c_B)
            r_N = []
            for j in self.idx_N:
                col_A = [mat_A[i][j] for i in range(self.m)]
                producte = sum(w[k] * col_A[k] for k in range(self.m))
                r_N.append(vec_c[j] - producte)

            candidats_entrada = [self.idx_N[k] for k, r in enumerate(r_N) if r < -self.tol]
            if not candidats_entrada:
                print(f"[jh_simplexP]     Iteració {self.comptador_iter:2d} : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z = {z_actual:8.3f}")
                return "Òptim"

            # Regla de Bland
            q = min(candidats_entrada)
            idx_q_en_N = self.idx_N.index(q)

            col_A_q = [mat_A[i][q] for i in range(self.m)]
            d_B = self._obte_d(col_A_q)

            if all(val >= -self.tol for val in d_B):
                return "No acotat"

            theta_star = float('inf')
            candidats_sortida = []
            for i in range(self.m):
                if d_B[i] < -self.tol:
                    theta = -x_B[i] / d_B[i]
                    if theta < theta_star - self.tol:
                        theta_star = theta
                        candidats_sortida = [self.idx_B[i]]
                    elif abs(theta - theta_star) <= self.tol:
                        candidats_sortida.append(self.idx_B[i])

            # Regla de Bland 
            p = min(candidats_sortida)
            idx_p_en_B = self.idx_B.index(p)

            print(f"[jh_simplexP]     Iteració {self.comptador_iter:2d} : iout = 0, q = {q+1:2d}, B(p) = {p+1:2d}, theta*= {theta_star:6.3f}, z = {z_actual:8.3f}")

            # Actualització de la inversa amb la matriu eta
            self._actualitza_B_inv(idx_p_en_B, d_B)

            # Actualització B i N
            self.idx_B[idx_p_en_B] = q
            self.idx_N[idx_q_en_N] = p


    def _imprimeix_solucio_final(self):
        c_B = [self.c_fase2[idx] for idx in self.idx_B]
        x_B = self._obte_x_B()
        z = float(sum(c_B[i] * x_B[i] for i in range(self.m)))

        print(f"[jh_simplexP]     Solució òptima trobada, iteració {self.comptador_iter:3d}, z = {z:.6f}")
        print(f"[jh_simplexP] Fi simplex primal")
        print()
        print("Solució òptima:")
        print()

        # Variables bàsiques originals i els seus valors
        vb_originals = [idx for idx in self.idx_B if idx < self.n]
        xb_originals = [x_B[i] for i, idx in enumerate(self.idx_B) if idx < self.n]

        # Costos reduïts de les variables no bàsiques originals
        w = self._obte_w(c_B)
        N_orig_ordenat = sorted([idx for idx in self.idx_N if idx < self.n])
        r_no_basiques = []
        for j in N_orig_ordenat:
            col_A = [self.A_ext[i][j] for i in range(self.m)]
            producte = sum(w[k] * col_A[k] for k in range(self.m))
            r_no_basiques.append(self.c_fase2[j] - producte)

        str_vb = " ".join([f"{idx+1:3d}" for idx in vb_originals])
        str_xb = " ".join([f"{val:5.2f}" for val in xb_originals])   
        str_r  = " ".join([f"{val:5.2f}" for val in r_no_basiques])  

        print(f"vb = {str_vb}")
        print(f"xb = {str_xb}")
        print(f"z  = {z:.6f}")  
        print(f"r  = {str_r}")



if __name__ == "__main__":

    c = [-35, -72, -21, -24, -70, -82, -68, -66, -26, -71, -50, -43, -94, -40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    A = [
        [33, 41, 76, 20, 26, 81, 58, 37, 57, 16, 67, 14, 64, 88, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [55, 100, 86, 78, 80, 83, 70, 32, 44, 94, 2, 83, 16, 19, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
        [58, 24, 13, 33, 32, 23, 63, 11, 77, 54, 49, 4, 98, 81, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
        [64, 14, 80, 37, 64, 66, 31, 29, 64, 86, 33, 69, 67, 93, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
        [81, 100, 81, 1, 44, 46, 43, 69, 71, 92, 100, 21, 13, 38, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
        [3, 36, 2, 12, 70, 5, 12, 0, 22, 57, 9, 48, 26, 47, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],
        [50, 55, 64, 20, 75, 77, 84, 96, 25, 30, 74, 95, 74, 51, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
        [97, 96, 31, 64, 81, 92, 97, 88, 49, 6, 38, 11, 89, 22, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
        [7, 22, 5, 55, 72, 8, 27, 55, 12, 41, 14, 26, 52, 80, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0],
        [42, 18, 14, 19, 50, 17, 47, 91, 78, 24, 54, 37, 11, 43, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]
    ]

    b = [677, 841, 619, 796, 799, 348, 869, 860, 475, 544]

    resoledor = ResolvedorSimplex(c, A, b)
    resoledor.resol()