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

        # Garantim que b[i] sigui positiu. Si no ho és, multipliquem tota la fila per -1
        for i in range(self.m):
            if self.b[i] < 0:
                self.b[i] *= -1.0
                for j in range(self.n):
                    self.A[i][j] *= -1.0


        self.base_trivial = self._check_trivial_base()

    # Recorre les últimes m columnes i comprovem que cada columna té exactament un 1 i m-1 zeros
    def _check_trivial_base(self):
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
    def _get_x_B(self):
        return self._mat_vec(self.B_inv, self.b)

    # Calculem el vector w = c_B^T B^-1
    def _get_w(self, c_B):
        return [sum(c_B[i] * self.B_inv[i][j] for i in range(self.m))
                for j in range(self.m)]

    # Calculem la direcció
    def _get_d(self, A_col_q):
        Binv_Aq = self._mat_vec(self.B_inv, A_col_q)
        return [-val for val in Binv_Aq]

    # Actualitzem la inversa a través de la matriu eta
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
        # B^-1 inicial = I
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

            w = self._get_w(c_B)
            r_N = []
            for j in self.N_idx:
                A_col = [A_mat[i][j] for i in range(self.m)]
                dot = sum(w[k] * A_col[k] for k in range(self.m))
                r_N.append(c_vec[j] - dot)

            candidates_in = [self.N_idx[k] for k, r in enumerate(r_N) if r < -self.tol]
            if not candidates_in:
                print(f"[jh_simplexP]     Iteració {self.iter_count:2d} : iout = 2, q =  0, B(p) =  0, theta*=  0.000, z = {z_current:8.3f}")
                return "Òptim"

            # Regla de Bland
            q = min(candidates_in)
            q_idx_in_N = self.N_idx.index(q)

            A_col_q = [A_mat[i][q] for i in range(self.m)]
            d_B = self._get_d(A_col_q)

            if all(val >= -self.tol for val in d_B):
                return "No acotat"

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