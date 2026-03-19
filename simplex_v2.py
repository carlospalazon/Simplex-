class SimplexSolver:
    def __init__(self, c, A, b):
        # Copia floats de c, A, b
        self.c = [float(val) for val in c]
        self.A = [[float(val) for val in row] for row in A]
        self.b = [float(val) for val in b]

        # Definim el nombre de restriccions
        self.m = len(self.A)
        # Definim el nombre de variables
        self.n = len(self.c)
        # Definim la tolerància, ja que els floats mai arriben a ser 0
        self.tol = 1e-8
        # Definim el número d'iteracions
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


    def _mat_vec(self, M, v):
        return [sum(M[i][k] * v[k] for k in range(self.m)) for i in range(self.m)]

    def _get_x_B(self):
        return self._mat_vec(self.B_inv, self.b)

    def _get_w(self, c_B):
        return [sum(c_B[i] * self.B_inv[i][j] for i in range(self.m))
                for j in range(self.m)]

    def _get_d(self, A_col_q):
        Binv_Aq = self._mat_vec(self.B_inv, A_col_q)
        return [-val for val in Binv_Aq]


    def _update_B_inv(self, p_row, d_B):
        # d_B[p] < 0 (condició de pivot)
        d_pq = d_B[p_row] 
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


    def solve(self):
        print("[jh_simplexP] Inici simplex primal amb regla de Bland")

        # Si no tenim una matriu identitat, pasem a la fase I
        if not self.base_trivial:
            print("[jh_simplexP]   Fase I")
            status = self._phase_one()
            # Si el problema és infactible, aturem
            if status == "Infactible":
                print("[jh_simplexP]     Problema infactible.")
                print("[jh_simplexP] Fi simplex primal")
                return
            print(f"[jh_simplexP]     Solució bàsica factible trobada, iteració {self.iter_count:3d}")
        else:
            # Si hi ha una identitat, pasem a la fase II
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


    def _phase_one(self):
        # Ampliem A amb variables artificials
        self.A_ext = []
        for i in range(self.m):
            row_ext = self.A[i][:] + [1.0 if i == j else 0.0 for j in range(self.m)]
            self.A_ext.append(row_ext)

        # Definim els nous costos amb les variables artificials
        self.c_ext = [0.0] * self.n + [1.0] * self.m
        # Definim els índexs de la nova base
        self.B_idx = list(range(self.n, self.n + self.m))
        # Definim els índexs de les variables no bàsiques
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


    def _phase_two(self):
        # Eliminem variables artificials no bàsiques
        self.N_idx = [idx for idx in self.N_idx if idx < self.n]
        # Cost 0 per a les artificials que poguessin quedar a la base (degenerat)
        self.c_phase2 = self.c[:] + [0.0] * self.m
        return self._simplex_core(self.A_ext, self.c_phase2, phase=2)


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


    def _print_final_solution(self):
        c_B = [self.c_phase2[idx] for idx in self.B_idx]
        x_B = self._get_x_B()
        z = float(sum(c_B[i] * x_B[i] for i in range(self.m)))

        print(f"[jh_simplexP]     Solució òptima trobada, iteració {self.iter_count:3d}, z = {z:.6f}")
        print(f"[jh_simplexP] Fi simplex primal")
        print()
        print("Solució òptima:")
        print()

        # Variables bàsiques originals i els seus valors
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
        xb_str = " ".join([f"{val:5.2f}" for val in xb_originals])   # Arrodonim variables bàsiques
        r_str  = " ".join([f"{val:5.2f}" for val in r_nonbasic])      # Arrodonim costos reduïts

        print(f"vb = {vb_str}")
        print(f"xb = {xb_str}")
        print(f"z  = {z:.6f}")  # Mostrem z complet
        print(f"r  = {r_str}")



if __name__ == "__main__":

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

    solver = SimplexSolver(c, A, b)
    solver.solve()