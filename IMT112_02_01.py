import numpy as np
from scipy.optimize import fsolve

# ==========================================
# SCRIPT 0: ANÁLISIS DE POSICIÓN (4 Barras)
# ==========================================

# 1. PARÁMETROS GEOMÉTRICOS (cm y radianes)
# L1 es el suelo (distancia entre O2 y O4)
L1, L2, L3, L4 = 20.0, 5.0, 15.0, 10.0 

# Ángulo de entrada del motor (El único dato que conocemos)
th2 = np.radians(60)

# 2. DEFINIR LAS ECUACIONES DEL LAZO VECTORIAL
# Queremos que estas dos ecuaciones sean iguales a CERO
def ecuaciones_lazo(incognitas):
    th3, th4 = incognitas
    
    # Ecuación en X: L2*cos(th2) + L3*cos(th3) - L1 - L4*cos(th4) = 0
    eq_x = L2 * np.cos(th2) + L3 * np.cos(th3) - L1 - L4 * np.cos(th4)
    
    # Ecuación en Y: L2*sin(th2) + L3*sin(th3) - 0 - L4*sin(th4) = 0
    eq_y = L2 * np.sin(th2) + L3 * np.sin(th3) - L4 * np.sin(th4)
    
    return [eq_x, eq_y]

# 3. LA "ADIVINANZA" INICIAL (Seed / Initial Guess)
# Los métodos numéricos necesitan una pista de por dónde empezar a buscar
# Asumimos visualmente que th3 está cerca de 20° y th4 cerca de 110°
semilla_inicial = [np.radians(25), np.radians(120)]

# 4. QUE LA COMPUTADORA HAGA LA MAGIA
solucion_pos = fsolve(ecuaciones_lazo, semilla_inicial)

# Convertimos de radianes a grados para que los humanos lo entiendan
th3_final_grados = np.degrees(solucion_pos[0])
th4_final_grados = np.degrees(solucion_pos[1])

print("--- 0. RESULTADOS DE POSICIÓN ---")
print(f"Ángulo del Acoplador (Theta 3): {th3_final_grados:.2f} grados")
print(f"Ángulo del Balancín (Theta 4):  {th4_final_grados:.2f} grados")