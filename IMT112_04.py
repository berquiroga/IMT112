import numpy as np

# ==========================================
# 1. DATOS DE ENTRADA (Calculados previamente)
# ==========================================
# Longitudes (cm)
L2, L3, L4 = 5.0, 12.0, 15.0 

# Ángulos en radianes (Analizados en Posición)
theta2 = np.radians(60)
theta3 = np.radians(20)   # Hipotético para este instante
theta4 = np.radians(110)  # Hipotético para este instante

# Velocidades angulares en rad/s (Analizadas en Velocidad)
w2 = 10.0  # Motor (Constante)
w3 = -2.5  # Acoplador
w4 = 3.2   # Balancín

# Aceleración de entrada del motor (rad/s^2)
# ¡OJO AQUÍ! Es cero porque el motor ya alcanzó su velocidad máxima y es constante
alpha2 = 0.0 

# ==========================================
# 2. CONSTRUCCIÓN DE MATRICES (A * X = B)
# ==========================================
# Incógnitas X = [alpha3, alpha4]

# Matriz A: Coeficientes de las incógnitas (Igual que en velocidad)
# Despejamos alpha3 y alpha4 al lado izquierdo de la ecuación
A = np.array([
    [-L3 * np.sin(theta3),  L4 * np.sin(theta4)],  # Fila del Eje X
    [ L3 * np.cos(theta3), -L4 * np.cos(theta4)]   # Fila del Eje Y
])

# Vector B: Todos los términos conocidos (Las aceleraciones centrípetas/normales)
# Mandamos todo lo que tiene w^2 y alpha2 al lado derecho
B_x = ( L2 * alpha2 * np.sin(theta2) + L2 * w2**2 * np.cos(theta2) + 
        L3 * w3**2 * np.cos(theta3)  - L4 * w4**2 * np.cos(theta4) )

B_y = (-L2 * alpha2 * np.cos(theta2) + L2 * w2**2 * np.sin(theta2) + 
        L3 * w3**2 * np.sin(theta3)  - L4 * w4**2 * np.sin(theta4) )

B = np.array([B_x, B_y])

# ==========================================
# 3. RESOLUCIÓN MATRICIAL
# ==========================================
# La computadora hace la magia algebraica pesada
solucion = np.linalg.solve(A, B)
alpha3, alpha4 = solucion

print("--- RESULTADOS DE ACELERACIÓN ---")
print(f"Aceleración Angular Biela (alpha3):    {alpha3:.2f} rad/s^2")
print(f"Aceleración Angular Balancín (alpha4): {alpha4:.2f} rad/s^2")