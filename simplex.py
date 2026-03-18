class SimplexSolver:
    def __init__(self, c, A, b):

        # Fem una copia profunda
        self.c = [float(val) for val in c]
        self.A = [[float(val) for val in row] for row in A]
        self.b = [float(val) for val in b]

        # Nombre de restriccions
        self.m = len(self.A)
        # Nombre de variables
        self.n = len(self.c)
        # Afegim una tolerància, ja que els floats no arriben a ser 0
        self.tol = 1e-8
        # Nombre d'iteracions
        self.iter_count = 0

        # Garantim termes independents positius
        for i in range(self.m):
            if self.b[i] < 0:
                self.b[i] *= -1.0
                for j in range(self.n):
                    self.A[i][j] *= -1.0

        # Detectar si ja hi ha base trivial (identitat)
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
    # Funcions per fer les operacions bàsiques
    # ------------------------------------------------------------------
    def _get_x_B(self):
        return self.b[:]

    def _get_w(self, c_B):
        return c_B[:]

    def _get_d(self, A_mat, q):
        return [A_mat[i][q] for i in range(self.m)]

    # ------------------------------------------------------------------
    # Funció solve 
    # ------------------------------------------------------------------
    def solve(self):
        print("[jh_simplexP] Inici simplex primal amb regla de Bland")

        # Si no tenim matriu identitat, anem a la fase I
        if not self.base_trivial:
            print("[jh_simplexP]   Fase I")
            status = self._phase_one()
            # Si la solució és infactible, aturem
            if status == "Infactible":
                print("[jh_simplexP]     Problema infactible.")
                return
            print(f"[jh_simplexP]     Solució bàsica factible trobada, iteració {self.iter_count}")
        # Si tenim matriu identitat, passem a la fase II directament
        else:
            print("[jh_simplexP] Base trivial detectada, saltant Fase I")

        print("[jh_simplexP]   Fase II")
        status = self._phase_two()
        if status == "No acotat":
            print("[jh_simplexP]     Problema no acotat.")
        # Si hem arribat a l'òptim, escrivim la solució
        elif status == "Òptim":
            self._print_final_solution()

    # ------------------------------------------------------------------
    # Fase I 
    # ------------------------------------------------------------------

    def _phase_one(self):

        # Definim una nova matriu A, que tindrà les variables artificials
        self.A_ext = []
        # Recorrem totes les restriccions
        for i in range(self.m):
            # Afegim les variables artificials a la fila
            row_ext = self.A[i][:] + [1.0 if i == j else 0.0 for j in range(self.m)]
            # Afegim la nova fila
            self.A_ext.append(row_ext)

        # Calculem els nous costos
        self.c_ext = [0.0] * self.n + [1.0] * self.m
        # Calculem els nous índexs
        self.B_idx = list(range(self.n, self.n + self.m))
        self.N_idx = list(range(self.n))
        # Canviem l'estat
        status = self._simplex_core(self.A_ext, self.c_ext, phase=1)

        # Definim els costos de les variables bàsiques
        c_B = [self.c_ext[j] for j in self.B_idx]
        # Calculem les variables bàsiques
        x_B = self._get_x_B()
        # Calculem z i comparem amb la tolerància (0)
        z = sum(c_B[i] * x_B[i] for i in range(self.m))
        # Si z > 0 el problema és infactible, ja que les variables artificials han de ser 0 per a que el problema sigui factible
        if z > self.tol:
            return "Infactible"
        
        # Si no es infactible, pues serà factible
        return "Factible"

    # ------------------------------------------------------------------
    # Fase II (igual que abans)
    # ------------------------------------------------------------------

    def _phase_two(self):
        # Eliminem les variables artificials que no són bàsiques
        self.N_idx = [idx for idx in self.N_idx if idx < self.n]
        # Calculem el cost. Per a les variables artificials que queden a la base assignem 0
        c_phase2 = self.c[:] + [0.0] * self.m
        # Cridem al nucli del símplex
        return self._simplex_core(self.A_ext, c_phase2, phase=2)

    # ------------------------------------------------------------------
    # Nucli del simplex (igual que abans)
    # ------------------------------------------------------------------

    def _simplex_core(self, A_mat, c_vec, phase=1):
        while True:
            # +1 iteració
            self.iter_count += 1
            # Obtenim els costos de les variables bàsiques
            c_B = [c_vec[j] for j in self.B_idx]
            # Calculem les solucions bàsiques
            x_B = self._get_x_B()
            # Calculem els multiplicadors
            w = self._get_w(c_B)

            # Costos reduïts
            r_N = []
            for j in self.N_idx:
                dot = sum(w[i] * A_mat[i][j] for i in range(self.m))
                r_N.append(c_vec[j] - dot)

            candidates_in = [self.N_idx[k] for k, r in enumerate(r_N) if r < -self.tol]
            if not candidates_in:
                z_final = 0.0 if phase == 1 else sum(c_B[i] * x_B[i] for i in range(self.m))
                print(f"[jh_simplexP] Iteració {self.iter_count:4d} : z = {z_final:8.3f}")
                return "Òptim"

            q = min(candidates_in)
            q_idx_in_N = self.N_idx.index(q)
            d = self._get_d(A_mat, q)
            if all(val <= self.tol for val in d):
                return "No acotat"

            theta_star = float('inf')
            candidates_out = []
            for i in range(self.m):
                if d[i] > self.tol:
                    theta = x_B[i] / d[i]
                    if theta < theta_star - self.tol:
                        theta_star = theta
                        candidates_out = [self.B_idx[i]]
                    elif abs(theta - theta_star) <= self.tol:
                        candidates_out.append(self.B_idx[i])

            p = min(candidates_out)
            p_idx_in_B = self.B_idx.index(p)

            self._pivot_update(A_mat, p_idx_in_B, d)

            self.B_idx[p_idx_in_B] = q
            self.N_idx[q_idx_in_N] = p

    def _pivot_update(self, A_mat, p_row, d):
        pivot = d[p_row]
        for j in range(len(A_mat[p_row])):
            A_mat[p_row][j] /= pivot
        self.b[p_row] /= pivot

        for i in range(self.m):
            if i == p_row:
                continue
            factor = d[i]
            if abs(factor) < self.tol:
                continue
            for j in range(len(A_mat[i])):
                A_mat[i][j] -= factor * A_mat[p_row][j]
            self.b[i] -= factor * self.b[p_row]

    def _print_final_solution(self):
        c_B = [self.c[idx] if idx < self.n else 0.0 for idx in self.B_idx]
        x_B = self._get_x_B()
        z = float(sum(c_B[i] * x_B[i] for i in range(self.m)))
        print(f"[jh_simplexP] Solució òptima trobada, iteració {self.iter_count}, z = {z:.6f}\n")
        vb_originals = [idx for idx in self.B_idx if idx < self.n]
        xb_originals = [x_B[i] for i, idx in enumerate(self.B_idx) if idx < self.n]
        vb_str = " ".join([f"{idx+1:2d}" for idx in vb_originals])
        xb_str = " ".join([f"{val:6.1f}" for val in xb_originals])
        print(f"vb = {vb_str}")
        print(f"xb = {xb_str}")
        print(f"z = {z:7.4f}")


# =======================================================
# EXECUCIÓ
# =======================================================

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