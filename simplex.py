class SimplexSolver:
    def __init__(self, c, A, b):
        self.c = [float(val) for val in c]
        self.A = [[float(val) for val in row] for row in A]
        self.b = [float(val) for val in b]
        
        self.m = len(self.A)
        self.n = len(self.c)
        self.tol = 1e-8
        self.iter_count = 0
        
        # Garantim termes independents positius (b >= 0)
        for i in range(self.m):
            if self.b[i] < 0:
                self.b[i] *= -1.0
                for j in range(self.n):
                    self.A[i][j] *= -1.0

    def solve(self):
        print("[jh_simplexP] Inici simplex primal amb regla de Bland")
        print("[jh_simplexP]   Fase I")
        status = self._phase_one()
        
        if status == "Infactible":
            print("[jh_simplexP]     Problema infactible.")
            return
            
        print(f"[jh_simplexP]     Solució bàsica factible trobada, iteració {self.iter_count}")
        print("[jh_simplexP]   Fase II")
        status = self._phase_two()
        
        if status == "No acotat":
            print("[jh_simplexP]     Problema no acotat.")
        elif status == "Òptim":
            self._print_final_solution()

    # ------------------------------------------------------------------
    # Operacions bàsiques assumint B = I
    # x_B = b  (directament)
    # w    = c_B  (directament)
    # d    = A_q  (directament, columna q de A)
    # ------------------------------------------------------------------

    def _get_x_B(self):
        """x_B = B^{-1} b = b  (B = I)"""
        return self.b[:]

    def _get_w(self, c_B):
        """w = B^{-T} c_B = c_B  (B = I)"""
        return c_B[:]

    def _get_d(self, A_mat, q):
        """d = B^{-1} A_q = A_q  (B = I), columna q de A_mat"""
        return [A_mat[i][q] for i in range(self.m)]

    # ------------------------------------------------------------------

    def _phase_one(self):
        # Crear el problema estès afegint variables artificials (matriu identitat)
        self.A_ext = []
        for i in range(self.m):
            row_ext = self.A[i][:] + [1.0 if i == j else 0.0 for j in range(self.m)]
            self.A_ext.append(row_ext)
            
        # Vector de costos Fase I (0 per originals, 1 per artificials)
        self.c_ext = [0.0] * self.n + [1.0] * self.m
        
        # Índexs inicials: base = variables artificials (columnes identitat)
        self.B_idx = list(range(self.n, self.n + self.m))
        self.N_idx = list(range(self.n))
        
        status = self._simplex_core(self.A_ext, self.c_ext, phase=1)
        
        # Avaluar factibilitat: z = c_B · x_B
        c_B = [self.c_ext[j] for j in self.B_idx]
        x_B = self._get_x_B()
        z = sum(c_B[i] * x_B[i] for i in range(self.m))
        
        if z > self.tol:
            return "Infactible"
        return "Factible"

    def _phase_two(self):
        # Eliminar variables artificials de no bàsiques
        self.N_idx = [idx for idx in self.N_idx if idx < self.n]
        
        # Cost original (zeros per artificials degenerades que puguin quedar a la base)
        c_phase2 = self.c[:] + [0.0] * self.m
        
        return self._simplex_core(self.A_ext, c_phase2, phase=2)

    def _simplex_core(self, A_mat, c_vec, phase=1):
        while True:
            self.iter_count += 1
            
            c_B = [c_vec[j] for j in self.B_idx]

            # 1. Solució bàsica i multiplicadors (trivials si B = I)
            x_B = self._get_x_B()
            w    = self._get_w(c_B)
            
            # 2. Costos reduïts per a les variables no bàsiques
            r_N = []
            for j in self.N_idx:
                dot = sum(w[i] * A_mat[i][j] for i in range(self.m))
                r_N.append(c_vec[j] - dot)
            
            candidates_in = [self.N_idx[k] for k, r in enumerate(r_N) if r < -self.tol]
            
            if not candidates_in:
                z_final = 0.0 if phase == 1 else sum(c_B[i] * x_B[i] for i in range(self.m))
                print(f"[jh_simplexP]     Iteració {self.iter_count:4d} : iout =  2, q =   0, "
                      f"B(p) =   0, theta*=  0.000, z = {z_final:8.3f}")
                return "Òptim"
            
            # 3. Regla de Bland: variable d'entrada
            q = min(candidates_in)
            q_idx_in_N = self.N_idx.index(q)
            
            # 4. Direcció de descens (trivial si B = I)
            d = self._get_d(A_mat, q)
            
            if all(val <= self.tol for val in d):
                return "No acotat"
            
            # 5. Test del quocient + Regla de Bland: variable de sortida
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
            
            # 6. Log de la iteració
            z = sum(c_B[i] * x_B[i] for i in range(self.m))
            print(f"[jh_simplexP]     Iteració {self.iter_count:4d} : iout =  0, q = {q+1:3d}, "
                  f"B(p) = {p+1:3d}, theta*= {theta_star:6.3f}, z = {z:8.3f}")
            
            # 7. Actualitzar base i actualitzar b = B^{-1} b per al nou B
            #    Quan B deixa de ser la identitat cal fer el pas de pivot
            #    sobre b i sobre les columnes de A_mat per mantenir B^{-1}·A i B^{-1}·b
            #    actualitzats (tableau explícit).
            self._pivot_update(A_mat, p_idx_in_B, d)
            
            self.B_idx[p_idx_in_B] = q
            self.N_idx[q_idx_in_N] = p

    def _pivot_update(self, A_mat, p_row, d):
        """
        Actualitza A_mat i self.b amb el pas de pivot de manera que
        la nova B segueixi sent la identitat (tableau explícit).
        
        Fila pivot: p_row
        Columna pivot: d[p_row]  (element de la direcció de descens)
        """
        pivot = d[p_row]
        
        # Normalitzar fila pivot
        for j in range(len(A_mat[p_row])):
            A_mat[p_row][j] /= pivot
        self.b[p_row] /= pivot
        
        # Eliminar la columna de la variable entrant en totes les altres files
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
        
        print(f"[jh_simplexP]     Solució òptima trobada, iteració  {self.iter_count}, z = {z:.6f}")
        print("[jh_simplexP] Fi simplex primal \n")
        print("Solució òptima:\n")
        
        # Multiplicadors i costos reduïts
        w = self._get_w(c_B)
        r = []
        for j in range(self.n):
            dot = sum(w[i] * self.A_ext[i][j] for i in range(self.m))
            r.append(self.c[j] - dot)
        
        vb_originals = [idx for idx in self.B_idx if idx < self.n]
        xb_originals = [x_B[i] for i, idx in enumerate(self.B_idx) if idx < self.n]
        
        vb_str = " ".join([f"{idx+1:2d}" for idx in vb_originals])
        xb_str = " ".join([f"{val:6.1f}" for val in xb_originals])
        r_str  = " ".join([f"{val:6.1f}" for val in r])
        
        print(f"vb = {vb_str}")
        print(f"xb = {xb_str}")
        print(f"z = {z:7.4f}")
        print(f"r = {r_str}")


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