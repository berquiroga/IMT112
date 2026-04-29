import numpy as np
from numpy import sin, cos

# ==========================================
# SCRIPT 1: ANÁLISIS DE VELOCIDAD (4 Barras)
# ==========================================
# 1. PARÁMETROS GEOMÉTRICOS (cm y radianes)
L2, L3, L4 = 5.0, 15.0, 10.0 
th2 = np.radians(60)
th3 = np.radians(20)   
th4 = np.radians(110)  

# Velocidad del Motor de entrada (rad/s)
w2 = 10.0  

# 2. SISTEMA DE ECUACIONES (A * X = B)
# Incógnitas X = [w3, w4]

# Matriz A: Coeficientes de la geometría (¡Atentos a esta matriz!)
A = np.array([
    [-L3 * sin(th3),  L4 * sin(th4)],  # Eje X
    [ L3 * cos(th3), -L4 * cos(th4)]   # Eje Y
])

# Vector B: Lo que conocemos (El empuje del motor)
B = np.array([
     L2 * w2 * sin(th2),
    -L2 * w2 * cos(th2)
])

# 3. RESOLUCIÓN
solucion_vel = np.linalg.solve(A, B)
w3, w4 = solucion_vel

print("--- 1. RESULTADOS DE VELOCIDAD ---")
print(f"Velocidad Angular Acoplador (w3): {w3:.2f} rad/s")
print(f"Velocidad Angular Balancín (w4):  {w4:.2f} rad/s")