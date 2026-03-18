class SimplexSolver:
    def __init__(self, c, A, b):
        self.c = [float(val) for val in c]
        self.A = [[float(val) for val in row] for row in A]
        self.b = [float(val) for val in b]
        self.m = len(self.A)
        self.n = len(self.c)
        self.tol = 1e-8
        self.iter_count = 0

        for i in range(self.m):
            if self.b[i] < 0:
                self.b[i] *= -1.0
                for j in range(self.n):
                    self.A[i][j] *= -1.0

        self.base_trivial = self._check_trivial_base()

    # ------------------------------------------------------------------
    # Funció per detectar identitat en A (base trivial)
    # ------------------------------------------------------------------
    def _check_trivial_base(self):
        for j in range(self.n - self.m, self.n):
            col = [self.A[i][j] for i in range(self.m)]
            if sum(1 for x in col if abs(x - 1) < self.tol) != 1:
                return False
            if sum(1 for x in col if abs(x) < self.tol) != self.m - 1:
                return False
        return True

    # ------------------------------------------------------------------
    # Funcions bàsiques amb B⁻¹ explícita
    # ------------------------------------------------------------------
    def _mat_vec(self, M, v):
        """Producte matriu M (m×m) per vector v (m)"""
        return [sum(M[i][k] * v[k] for k in range(self.m)) for i in range(self.m)]

    def _get_x_B(self):
        # x_B = B⁻¹ b
        return self._mat_vec(self.B_inv, self.b)

    def _get_w(self, c_B):
        # w' = C'_B · B⁻¹  →  w[j] = Σ_i c_B[i] * B_inv[i][j]
        return [sum(c_B[i] * self.B_inv[i][j] for i in range(self.m))
                for j in range(self.m)]

    def _get_d(self, A_col_q):
        # d_B = -B⁻¹ A_q  (Imatge 1, pas 3.1)
        Binv_Aq = self._mat_vec(self.B_inv, A_col_q)
        return [-val for val in Binv_Aq]

    # ------------------------------------------------------------------
    # Actualització de B⁻¹ via matriu eta (Imatge 2)
    # B⁻¹_actual = E · B⁻¹_prèvia
    # η_p(i) = -d_iq / d_pq  per i ≠ p
    # η_p(p) = -1 / d_pq
    # ------------------------------------------------------------------
    def _update_B_inv(self, p_row, d_B):
        d_pq = d_B[p_row]  # d_B[p] < 0 (condició de pivot)
        new_B_inv = []
        for i in range(self.m):
            if i == p_row:
                # η_p(p) = -1 / d_pq
                eta_p = -1.0 / d_pq
                new_row = [eta_p * self.B_inv[p_row][j] for j in range(self.m)]
            else:
                # η_p(i) = -d_iq / d_pq
                eta_i = -d_B[i] / d_pq
                new_row = [self.B_inv[i][j] + eta_i * self.B_inv[p_row][j]
                           for j in range(self.m)]
            new_B_inv.append(new_row)
        self.B_inv = new_B_inv

    # ------------------------------------------------------------------
    # Funció solve
    # ------------------------------------------------------------------
    def solve(self):
        print("[jh_simplexP] Inici simplex primal amb regla de Bland")

        if not self.base_trivial:
            print("[jh_simplexP]   Fase I")
            status = self._phase_one()
            if status == "Infactible":
                print("[jh_simplexP]     Problema infactible.")
                print("[jh_simplexP] Fi simplex primal")
                return
            print(f"[jh_simplexP]     Solució bàsica factible trobada, iteració {self.iter_count:3d}")
        else:
            print("[jh_simplexP] Base trivial detectada, saltant Fase I")
            self.B_inv = [[1.0 if i == j else 0.0 for j in range(self.m)]
                          for i in range(self.m)]
            self.B_idx = list(range(self.n - self.m, self.n))
            self.N_idx = list(range(self.n - self.m))
            self.A_ext = [row[:] for row in self.A]

        print("[jh_simplexP]   Fase II")
        status = self._phase_two()
        if status == "No acotat":
            print("[jh_simplexP]     Problema no acotat.")
            print("[jh_simplexP] Fi simplex primal")
        elif status == "Òptim":
            self._print_final_solution()

    # ------------------------------------------------------------------
    # Fase I
    # ------------------------------------------------------------------
    def _phase_one(self):
        # Ampliem A amb variables artificials → B inicial = I
        self.A_ext = []
        for i in range(self.m):
            row_ext = self.A[i][:] + [1.0 if i == j else 0.0 for j in range(self.m)]
            self.A_ext.append(row_ext)

        self.c_ext = [0.0] * self.n + [1.0] * self.m
        self.B_idx = list(range(self.n, self.n + self.m))
        self.N_idx = list(range(self.n))
        # B⁻¹ inicial = I
        self.B_inv = [[1.0 if i == j else 0.0 for j in range(self.m)]
                      for i in range(self.m)]

        self._simplex_core(self.A_ext, self.c_ext, phase=1)

        c_B = [self.c_ext[j] for j in self.B_idx]
        x_B = self._get_x_B()
        z = sum(c_B[i] * x_B[i] for i in range(self.m))

        # z*_I > 0 → infactible (Imatge 4)
        if z > self.tol:
            return "Infactible"

        # z*_I = 0 però hi ha variables artificials a la base (cas degenerat)
        # → intentar pivotar-les fora (Imatge 4)
        for p_row, art_idx in [(i, idx) for i, idx in enumerate(self.B_idx) if idx >= self.n]:
            for j in sorted(self.N_idx):
                if j >= self.n:
                    continue
                A_col = [self.A_ext[i][j] for i in range(self.m)]
                d_B = self._get_d(A_col)
                if abs(d_B[p_row]) > self.tol:
                    self._update_B_inv(p_row, d_B)
                    q_idx_in_N = self.N_idx.index(j)
                    self.B_idx[p_row] = j
                    self.N_idx[q_idx_in_N] = art_idx
                    break

        return "Factible"

    # ------------------------------------------------------------------
    # Fase II
    # ------------------------------------------------------------------
    def _phase_two(self):
        # Eliminem variables artificials no bàsiques
        self.N_idx = [idx for idx in self.N_idx if idx < self.n]
        # Cost 0 per a les artificials que poguessin quedar a la base (degenerat)
        self.c_phase2 = self.c[:] + [0.0] * self.m
        return self._simplex_core(self.A_ext, self.c_phase2, phase=2)

    # ------------------------------------------------------------------
    # Nucli del símplex
    # ------------------------------------------------------------------
    def _simplex_core(self, A_mat, c_vec, phase=1):
        while True:
            self.iter_count += 1
            c_B = [c_vec[j] for j in self.B_idx]
            x_B = self._get_x_B()
            z_current = sum(c_B[i] * x_B[i] for i in range(self.m))

            # Pas 2.1: costos reduïts r' = C'_N - C'_B B⁻¹ A_N
            w = self._get_w(c_B)
            r_N = []
            for j in self.N_idx:
                A_col = [A_mat[i][j] for i in range(self.m)]
                dot = sum(w[k] * A_col[k] for k in range(self.m))
                r_N.append(c_vec[j] - dot)

            # Pas 2.2: si r' ≥ [0] → SBF òptima (iout=2)
            candidates_in = [self.N_idx[k] for k, r in enumerate(r_N) if r < -self.tol]
            if not candidates_in:
                print(f"[jh_simplexP]     Iteració {self.iter_count:2d} : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z = {z_current:8.3f}")
                return "Òptim"

            # Regla de Bland: variable d'entrada amb índex mínim (Imatge 3)
            q = min(candidates_in)
            q_idx_in_N = self.N_idx.index(q)

            # Pas 3.1: d_B = -B⁻¹ A_q
            A_col_q = [A_mat[i][q] for i in range(self.m)]
            d_B = self._get_d(A_col_q)

            # Pas 3.2: si d_B ≥ [0] → no acotat
            if all(val >= -self.tol for val in d_B):
                return "No acotat"

            # Pas 4.1: θ* = min{i|d_B(i)<0} { -x_B(i)/d_B(i) }
            theta_star = float('inf')
            candidates_out = []
            for i in range(self.m):
                if d_B[i] < -self.tol:
                    theta = -x_B[i] / d_B[i]
                    if theta < theta_star - self.tol:
                        theta_star = theta
                        candidates_out = [self.B_idx[i]]
                    elif abs(theta - theta_star) <= self.tol:
                        candidates_out.append(self.B_idx[i])

            # Pas 4.2: Bland en empat (Imatge 3)
            p = min(candidates_out)
            p_idx_in_B = self.B_idx.index(p)

            print(f"[jh_simplexP]     Iteració {self.iter_count:2d} : iout = 0, q = {q+1:2d}, B(p) = {p+1:2d}, theta*= {theta_star:6.3f}, z = {z_current:8.3f}")

            # Pas 5.1: actualitzar B⁻¹ via matriu eta (Imatge 2)
            self._update_B_inv(p_idx_in_B, d_B)

            # Pas 5.2: actualitzar conjunts B i N
            self.B_idx[p_idx_in_B] = q
            self.N_idx[q_idx_in_N] = p

    # ------------------------------------------------------------------
    # Impressió de la solució final
    # ------------------------------------------------------------------
    def _print_final_solution(self):
        c_B = [self.c_phase2[idx] for idx in self.B_idx]
        x_B = self._get_x_B()
        z = float(sum(c_B[i] * x_B[i] for i in range(self.m)))

        print(f"[jh_simplexP]     Solució òptima trobada, iteració {self.iter_count:3d}, z = {z:.6f}")
        print(f"[jh_simplexP] Fi simplex primal")
        print()
        print("Solució òptima:")
        print()

        vb_originals = [idx for idx in self.B_idx if idx < self.n]
        xb_originals = [x_B[i] for i, idx in enumerate(self.B_idx) if idx < self.n]

        # Costos reduïts de les variables no bàsiques originals
        w = self._get_w(c_B)
        N_orig_sorted = sorted([idx for idx in self.N_idx if idx < self.n])
        r_nonbasic = []
        for j in N_orig_sorted:
            A_col = [self.A_ext[i][j] for i in range(self.m)]
            dot = sum(w[k] * A_col[k] for k in range(self.m))
            r_nonbasic.append(self.c_phase2[j] - dot)

        vb_str = " ".join([f"{idx+1:3d}" for idx in vb_originals])
        xb_str = " ".join([f"{val:5.1f}" for val in xb_originals])
        r_str  = " ".join([f"{val:5.1f}" for val in r_nonbasic])

        print(f"vb = {vb_str}")
        print(f"xb = {xb_str}")
        print(f"z = {z:.1f}")
        print(f"r = {r_str}")



# =======================================================
# EXECUCIÓ - Dades de prova de l'enunciat (Problema 1)
# =======================================================
if __name__ == "__main__":

    ## Problema inicial 

    c = [-85, -65, -34, 92, -83, 93, 22, 76, 48, -71, 25, -50, 79, -99, 0, 0, 0, 0, 0, 0]
    A = [
        [-13, -29, 23, -27, 57, -52, -70, 81, -2, 64, -54, 9, 68, 10, 0, 0, 0, 0, 0, 0],
        [-6, -53, -37, 36, 85, -22, 35, 87, 32, 34, 65, -64, 65, -25, 0, 0, 0, 0, 0, 0],
        [-13, 93, 31, 16, 45, 23, -11, 83, 75, -65, 67, 10, -92, 88, 0, 0, 0, 0, 0, 0],
        [67, -8, -27, -22, 47, 44, -3, -23, -16, 81, -24, -38, -67, 24, 0, 0, 0, 0, 0, 0],
        [80, 57, 77, 72, 52, 58, 53, 84, 77, 92, 57, 75, 68, 88, 1, 0, 0, 0, 0, 0],
        [-81, 21, -14, 73, 11, -10, -7, 45, 95, -92, 97, -79, 89, -92, 0, -1, 0, 0, 0, 0],
        [96, -7, 85, -63, -24, -95, 42, -71, 100, 42, 43, -23, 83, 60, 0, 0, -1, 0, 0, 0],
        [-27, 92, -84, 77, 68, 73, 28, 53, 34, -71, 59, -17, -36, -1, 0, 0, 0, 1, 0, 0],
        [-1, -28, 96, 74, 66, -82, -24, 91, 63, 65, 64, -29, 73, 93, 0, 0, 0, 0, -1, 0],
        [-60, -34, 21, -65, 84, 4, 12, 98, -38, -1, -99, 25, 21, 82, 0, 0, 0, 0, 0, -1]
    ]
    b = [65, 232, 350, 35, 991, 55, 267, 249, 520, 49]

        ## Problema 1 conjunto 34
    c = [8, -5, 6, 77, -77, 96, -74, -37, 91, -100, -10, -59, 28, 85, 0, 0, 0, 0, 0, 0]
    A = [
        [-67, 13, -26, 24, -32, 36, -36, 46, 11, 77, -64, -86, 90, 42, 0, 0, 0, 0, 0, 0],
        [-42, -70, 30, -74, 59, 35, -13, 6, 86, -54, 12, 91, 68, 82, 0, 0, 0, 0, 0, 0],
        [97, 80, -70, 45, 45, -34, 22, 55, -39, -20, 70, 31, 58, -27, 0, 0, 0, 0, 0, 0],
        [-45, -33, -35, 44, 34, 6, 100, -65, 30, 56, 38, -25, -5, -46, 0, 0, 0, 0, 0, 0],
        [59, 60, 92, 58, 83, 75, 74, 87, 70, 59, 63, 97, 52, 85, 1, 0, 0, 0, 0, 0],
        [-40, 89, -9, 79, 32, 3, 44, -36, 25, -86, 30, -33, 24, -84, 0, 1, 0, 0, 0, 0],
        [82, -24, -100, 30, -74, 1, 95, -57, -49, 55, -89, 90, 45, 50, 0, 0, -1, 0, 0, 0],
        [92, 97, -44, 35, 35, 48, 77, -17, -9, 74, 16, 70, -59, 29, 0, 0, 0, 1, 0, 0],
        [70, 22, 2, 69, 4, 92, -74, -83, 56, -69, 88, 69, 32, 69, 0, 0, 0, 0, 1, 0],
        [-42, 30, 98, 94, 8, -46, 70, -50, -14, -38, 30, 51, 86, 25, 0, 0, 0, 0, 0, 1]
    ]
    b = [28, 216, 313, 54, 1015, 39, 54, 445, 348, 303]

    ## Problema 2 conjunt 34

    c = [-82, -63, 73, 15, -51, 69, 100, 94, 61, 81, 18, 41, -33, -46, 0, 0, 0, 0, 0, 0]
    A = [
        [53, -61, -53, 23, 7, 17, -23, -87, 64, 18, 73, 33, -20, -38, 0, 0, 0, 0, 0, 0],
        [-45, -81, 63, -85, 38, -67, 87, 53, -61, 41, 11, 66, -23, 82, 0, 0, 0, 0, 0, 0],
        [-74, 91, 16, 4, -31, -42, 51, 38, -75, 94, 69, -51, 81, 76, 0, 0, 0, 0, 0, 0],
        [-65, -75, 85, -22, -76, 95, -13, 2, 57, 68, 10, 78, 57, -54, 0, 0, 0, 0, 0, 0],
        [69, 57, 94, 77, 98, 63, 74, 51, 51, 74, 58, 53, 78, 68, 1, 0, 0, 0, 0, 0],
        [65, -78, 75, -14, 38, -23, -59, -29, -67, 61, 74, -6, -6, 75, 0, -1, 0, 0, 0, 0],
        [-38, -21, -10, 10, -66, 75, -11, 34, 37, 17, 74, 10, 75, 99, 0, 0, -1, 0, 0, 0],
        [58, 27, -10, -23, 60, 2, 72, -90, 47, -77, -67, 24, 64, -80, 0, 0, 0, 1, 0, 0],
        [26, -33, 66, 4, 29, 11, 62, 31, 47, -54, 79, 89, -98, -37, 0, 0, 0, 0, 1, 0],
        [20, 37, 96, 28, -95, 32, -31, 43, 41, 49, -22, -90, 68, -11, 0, 0, 0, 0, 0, 1]
    ]
    b = [6, 79, 247, 147, 966, 105, 284, 8, 223, 166]

    ## Problema 3 conjunt 34

    c = [29, 80, -61, -85, 11, 78, -1, 60, 67, 97, 46, -62, -94, 24, 0, 0, 0, 0, 0, 0]
    A = [
        [28, -60, 0, -7, -3, -15, -69, -3, -22, 72, 57, 89, 92, 28, 0, 0, 0, 0, 0, 0],
        [-94, 92, 67, -97, -94, 10, 83, 2, 46, -61, -89, 56, 73, 43, 0, 0, 0, 0, 0, 0],
        [30, 83, 18, 99, -24, 7, 28, -99, -57, 38, 21, 29, 78, -11, 0, 0, 0, 0, 0, 0],
        [-24, -99, -5, 51, 41, -61, -55, 50, -91, 84, 23, 43, 67, 8, 0, 0, 0, 0, 0, 0],
        [-58, 31, -52, -58, 75, -76, 69, -31, -10, 31, 81, 38, 40, 36, -1, 0, 0, 0, 0, 0],
        [-37, -64, 44, -47, 36, 46, 83, -94, 40, -3, -92, 37, 96, -21, 0, -1, 0, 0, 0, 0],
        [31, -20, 80, 14, 93, 99, -54, 7, 80, -69, 78, -57, 77, -9, 0, 0, 1, 0, 0, 0],
        [80, 64, 96, 74, 34, 37, 34, 6, 89, 96, 96, 7, 60, 6, 0, 0, 0, 1, 0, 0],
        [55, 34, 36, 78, 94, 2, 60, 99, 9, 76, 7, 23, 84, 6, 0, 0, 0, 0, 1, 0],
        [25, 33, 19, 91, 15, 47, 99, 15, 94, 99, 78, 18, 54, 26, 0, 0, 0, 0, 0, 1]
    ]
    b = [187, 37, 240, 32, 115, 23, 351, 29, 97, 10]

    ## Problema 4 conjunt 34

    c = [-89, -86, -1, -57, -99, -53, -51, -81, -8, -5, -5, -2, -60, -75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    A = [
        [34, 39, 23, 84, 91, 84, 95, 67, 52, 57, 7, 86, 86, 16, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [21, 38, 87, 2, 96, 75, 47, 83, 17, 58, 7, 35, 15, 100, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
        [25, 61, 5, 61, 80, 51, 47, 30, 97, 63, 11, 63, 79, 15, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
        [74, 43, 10, 7, 51, 21, 47, 55, 64, 92, 94, 92, 97, 74, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
        [21, 72, 63, 4, 65, 83, 16, 83, 79, 100, 6, 16, 43, 62, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
        [18, 66, 30, 48, 82, 79, 51, 30, 20, 67, 23, 32, 27, 91, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],
        [68, 68, 60, 60, 31, 2, 4, 71, 47, 66, 11, 95, 58, 31, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
        [97, 26, 77, 51, 45, 91, 55, 42, 49, 47, 90, 25, 4, 93, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0],
        [17, 41, 69, 46, 72, 94, 45, 53, 17, 19, 10, 58, 6, 9, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0],
        [60, 87, 45, 69, 61, 66, 30, 93, 52, 93, 96, 74, 67, 24, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1]
    ]
    b = [820, 680, 687, 820, 712, 663, 671, 791, 555, 916]

    ## Problema 1 conjunto 37

    c = [53, -65, 31, 70, 30, 63, 74, -29, -1, -24, -37, -6, 94, 4, 0, 0, 0, 0, 0, 0]
    A = [
        [-87, -65, -31, 84, 63, 97, 96, -79, -73, -93, 35, 93, 72, 83, 0, 0, 0, 0, 0, 0],
        [45, 12, 54, 82, 89, 16, -100, 83, -86, -38, 2, 88, 24, -74, 0, 0, 0, 0, 0, 0],
        [84, -50, -67, 62, -72, 89, 33, -74, 55, 31, 33, -60, 23, -42, 0, 0, 0, 0, 0, 0],
        [39, 28, 42, 27, 19, -36, -81, 73, -13, 79, 99, 76, -71, 48, 0, 0, 0, 0, 0, 0],
        [91, 52, 51, 86, 62, 69, 92, 78, 70, 67, 76, 75, 72, 50, 1, 0, 0, 0, 0, 0],
        [-92, 87, -45, 15, -24, 62, 13, -12, 19, -75, 90, -44, 100, 41, 0, -1, 0, 0, 0, 0],
        [97, -48, 62, 0, 73, 57, 23, 59, -100, 32, -42, -51, -63, -93, 0, 0, 1, 0, 0, 0],
        [64, -34, 19, 61, -44, 25, -6, -10, 33, -89, 61, 96, 5, -94, 0, 0, 0, -1, 0, 0],
        [-58, 89, 6, -72, 96, 30, 40, 91, -2, -97, 38, 91, 22, -32, 0, 0, 0, 0, -1, 0],
        [47, -22, -48, 29, 64, 36, -7, -89, -91, 80, 28, 73, -29, -20, 0, 0, 0, 0, 0, 1]
    ]
    b = [195, 197, 45, 329, 992, 134, 7, 86, 241, 52]

    ## Problema 2 conjunto 37
    c = [-37, 78, 98, 8, 56, 36, 11, -99, -31, -8, -9, -71, 33, 74, 0, 0, 0, 0, 0, 0]
    A = [
        [6, -9, 4, 80, -24, 49, -83, -54, 90, -2, -99, 26, 19, 38, 0, 0, 0, 0, 0, 0],
        [69, 38, 39, -76, 8, -48, -27, -24, -90, 25, -15, 69, 21, 37, 0, 0, 0, 0, 0, 0],
        [-18, 61, -48, 98, -4, 68, -4, 57, -91, 83, -34, -22, 0, -56, 0, 0, 0, 0, 0, 0],
        [35, 80, -92, -49, 23, 100, -6, -10, -74, -38, 90, 72, 19, 52, 0, 0, 0, 0, 0, 0],
        [92, 100, 53, 54, 86, 57, 50, 93, 51, 91, 71, 82, 98, 93, 1, 0, 0, 0, 0, 0],
        [-67, 98, -15, 80, 46, 42, 34, -77, 13, -100, 29, 82, -83, 32, 0, -1, 0, 0, 0, 0],
        [52, -87, -28, 20, 65, -17, -36, -32, 15, 70, -4, 29, -93, 95, 0, 0, 1, 0, 0, 0],
        [98, 36, -15, 83, -69, 35, 35, 63, -23, 62, -57, -59, 47, 15, 0, 0, 0, -1, 0, 0],
        [14, 93, 94, -29, -71, 54, 96, 23, -43, 76, -47, 94, 49, -74, 0, 0, 0, 0, 1, 0],
        [56, 15, 50, 37, 39, -78, 72, -4, 36, 34, 24, 68, 47, 56, 0, 0, 0, 0, 0, -1]
    ]
    b = [41, 26, 90, 202, 1072, 113, 50, 250, 330, 451]

    ## Problema 3 conjunto 37

    c = [74, 20, -36, -92, -83, 45, -54, 68, -61, -28, 19, -9, -28, -57, 0, 0, 0, 0, 0, 0]
    A = [
        [75, -69, 12, -31, 5, 56, -42, -59, 36, 91, 24, 62, 1, 12, 0, 0, 0, 0, 0, 0],
        [20, -97, 92, 40, 52, 40, -65, -30, -82, 70, 35, 64, 29, -4, 0, 0, 0, 0, 0, 0],
        [-50, 42, 90, 59, -84, 98, -22, 80, -65, 19, -40, 97, 71, 77, 0, 0, 0, 0, 0, 0],
        [54, 83, 68, 58, -32, 94, 47, 60, 45, 70, -15, 96, -79, -85, 0, 0, 0, 0, 0, 0],
        [-63, 88, -74, -89, 19, -96, 95, 37, 70, 40, -45, -67, 14, 89, 1, 0, 0, 0, 0, 0],
        [32, -70, -8, -9, -35, 15, 92, -28, -11, 31, 97, 89, 19, 93, 0, -1, 0, 0, 0, 0],
        [45, -25, 23, -11, 98, -96, 46, -58, 25, -28, 80, -40, 61, 12, 0, 0, 1, 0, 0, 0],
        [58, 69, 51, 25, 5, 60, 86, 43, 70, 48, 50, 88, 84, 28, 0, 0, 0, 1, 0, 0],
        [75, 97, 50, 74, 95, 46, 25, 69, 88, 60, 63, 58, 27, 36, 0, 0, 0, 0, 1, 0],
        [66, 94, 57, 32, 74, 31, 79, 14, 23, 14, 6, 100, 34, 7, 0, 0, 0, 0, 0, 1]
    ]
    b = [173, 164, 372, 464, 19, 306, 133, 34, 66, 70]

    ## Problema 4 conjunto 37
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


    solver = SimplexSolver(c, A, b)
    solver.solve()
