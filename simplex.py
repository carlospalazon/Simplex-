class SimplexSolver:
    def __init__(self, c, A, b):
        # Convertim les dades a llistes de floats (còpies profundes per no modificar l'original)
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

    def _solve_linear_system(self, A_mat, b_vec):
        """Resol Ax = b utilitzant el mètode de Gauss-Jordan amb pivoteig parcial."""
        n = len(A_mat)
        # Matriu ampliada [A | b]
        M = [row[:] + [b_vec[i]] for i, row in enumerate(A_mat)]
        
        for i in range(n):
            # Pivoteig parcial
            max_row = max(range(i, n), key=lambda r: abs(M[r][i]))
            M[i], M[max_row] = M[max_row], M[i]
            
            pivot = M[i][i]
            if abs(pivot) < 1e-10:
                raise ValueError("Matriu singular o gairebé singular detectada.")
                
            # Fer 1 a la diagonal
            M[i] = [val / pivot for val in M[i]]
            
            # Fer zeros a la resta de la columna i
            for j in range(n):
                if i != j:
                    factor = M[j][i]
                    M[j] = [M[j][k] - factor * M[i][k] for k in range(n + 1)]
                    
        # Retornem l'última columna (la solució)
        return [M[i][n] for i in range(n)]

    def _phase_one(self):
        # 1. Crear el problema estès afegint variables artificials (Matriu Identitat)
        self.A_ext = []
        for i in range(self.m):
            row_ext = self.A[i][:] + [1.0 if i == j else 0.0 for j in range(self.m)]
            self.A_ext.append(row_ext)
            
        # 2. Vector de costos Fase I (0 per originals, 1 per artificials)
        self.c_ext = [0.0] * self.n + [1.0] * self.m
        
        # 3. Índexs inicials
        self.B_idx = list(range(self.n, self.n + self.m))
        self.N_idx = list(range(self.n))
        
        # 4. Executar Fase I
        status = self._simplex_core(self.A_ext, self.c_ext, phase=1)
        
        # 5. Avaluar factibilitat
        B_mat = [[self.A_ext[i][j] for j in self.B_idx] for i in range(self.m)]
        x_B = self._solve_linear_system(B_mat, self.b)
        c_B = [self.c_ext[j] for j in self.B_idx]
        
        z = sum(c_B[i] * x_B[i] for i in range(self.m))
        
        if z > self.tol:
            return "Infactible"
            
        return "Factible"

    def _phase_two(self):
        # 1. Eliminar variables artificials de la llista de no bàsiques
        self.N_idx = [idx for idx in self.N_idx if idx < self.n]
        
        # 2. Estendre el vector de costos original amb zeros 
        # (per tolerar artificials degenerades que hagin quedat a la base)
        c_phase2 = self.c[:] + [0.0] * self.m
        
        # 3. Executar Fase II usant l'A estesa però amb el cost original
        return self._simplex_core(self.A_ext, c_phase2, phase=2)

    def _simplex_core(self, A_mat, c_vec, phase=1):
        while True:
            self.iter_count += 1
            
            # 1. Construir matrius
            B_mat = [[A_mat[i][j] for j in self.B_idx] for i in range(self.m)]
            c_B = [c_vec[j] for j in self.B_idx]
            
            # 2. Resoldre sistemes
            x_B = self._solve_linear_system(B_mat, self.b)
            
            B_mat_T = [[B_mat[j][i] for j in range(self.m)] for i in range(self.m)]
            w = self._solve_linear_system(B_mat_T, c_B)
            
            # 3. Càlcul de costos reduïts (r_N) manualment
            r_N = []
            for j in self.N_idx:
                dot_product = sum(w[i] * A_mat[i][j] for i in range(self.m))
                r_N.append(c_vec[j] - dot_product)
            
            # Cerca de candidats d'entrada
            candidates_in = [self.N_idx[idx] for idx, r in enumerate(r_N) if r < -self.tol]
            
            if not candidates_in:
                # Impressió final iteració "buida" quan ja és òptim
                if phase == 1:
                    z_final = -0.000
                else:
                    z_final = sum(c_B[i] * x_B[i] for i in range(self.m))
                    
                print(f"[jh_simplexP]     Iteració {self.iter_count:4d} : iout =  2, q =   0, "
                      f"B(p) =   0, theta*=  0.000, z = {z_final:8.3f}")
                return "Òptim"
            
            # 4. Regla de Bland: Variable d'Entrada
            q = min(candidates_in)
            q_idx_in_N = self.N_idx.index(q)
            
            # 5. Direcció de descens
            A_q = [A_mat[i][q] for i in range(self.m)]
            d = self._solve_linear_system(B_mat, A_q)
            
            # Test no acotat
            if all(val <= self.tol for val in d):
                return "No acotat"
            
            # 6. Test del quocient i Regla de Bland: Variable de Sortida
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
            
            # 7. Imprimir registre de la iteració
            z = sum(c_B[i] * x_B[i] for i in range(self.m))
            print(f"[jh_simplexP]     Iteració {self.iter_count:4d} : iout =  0, q = {q+1:3d}, "
                  f"B(p) = {p+1:3d}, theta*= {theta_star:6.3f}, z = {z:8.3f}")
            
            # 8. Actualitzar base
            self.B_idx[p_idx_in_B] = q
            self.N_idx[q_idx_in_N] = p

    def _print_final_solution(self):
        # Avaluem sobre la matriu original estesa (fase 2)
        B_mat = [[self.A_ext[i][j] for j in self.B_idx] for i in range(self.m)]
        x_B = self._solve_linear_system(B_mat, self.b)
        
        # Filtrem els costos (només ens interessen les originals)
        c_B = [self.c[idx] if idx < self.n else 0.0 for idx in self.B_idx]
        z = float(sum(c_B[i] * x_B[i] for i in range(self.m)))
        
        print(f"[jh_simplexP]     Solució òptima trobada, iteració  {self.iter_count}, z = {z:.6f}")
        print("[jh_simplexP] Fi simplex primal \n")
        print("Solució òptima:\n")
        
        # Multiplicadors del símplex i costos reduïts globals
        B_mat_T = [[B_mat[j][i] for j in range(self.m)] for i in range(self.m)]
        w = self._solve_linear_system(B_mat_T, c_B)
        
        r = []
        for j in range(self.n):
            dot_product = sum(w[i] * self.A[i][j] for i in range(self.m))
            r.append(self.c[j] - dot_product)
        
        # Filtrem variables originals a la base
        vb_originals = [idx for idx in self.B_idx if idx < self.n]
        xb_originals = [x_B[i] for i, idx in enumerate(self.B_idx) if idx < self.n]
        
        # Formateig final
        vb_str = " ".join([f"{idx+1:2d}" for idx in vb_originals])
        xb_str = " ".join([f"{val:6.1f}" for val in xb_originals])
        r_str = " ".join([f"{val:6.1f}" for val in r])
        
        print(f"vb = {vb_str}")
        print(f"xb = {xb_str}")
        print(f"z = {z:7.4f}")
        print(f"r = {r_str}")

# =======================================================
# EXECUCIÓ - Dades de prova de l'enunciat (Problema 1)
# =======================================================
if __name__ == "__main__":
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
