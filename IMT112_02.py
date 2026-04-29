import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. PARÁMETROS DEL MECANISMO (Biela-Manivela)
# ==========================================
L2 = 5.0   # Longitud de la Manivela (cm)
L3 = 12.0  # Longitud de la Biela (cm)
offset = 0 # Desfase vertical de la corredera (cm)

# Vector de ángulos de la manivela (0 a 360 grados)
theta2 = np.linspace(0, 2*np.pi, 100)

# ==========================================
# 2. CÁLCULO DE POSICIÓN (Ecuaciones de Lazo)
# ==========================================
# Coordenadas del Punto A (Unión Manivela-Biela)
x_A = L2 * np.cos(theta2)
y_A = L2 * np.sin(theta2)

# Coordenadas del Punto B (El Pistón / Corredera)
# Usando el Teorema de Pitágoras en el triángulo formado por la biela
x_B = x_A + np.sqrt(L3**2 - (y_A - offset)**2)
y_B = np.full_like(x_B, offset) # La corredera siempre está en Y = offset

# ==========================================
# 3. GRAFICACIÓN DEL MOVIMIENTO
# ==========================================
plt.figure(figsize=(10, 5))

# Graficar la trayectoria del Pasador A (Un círculo)
plt.plot(x_A, y_A, 'b--', label='Trayectoria Manivela (A)', alpha=0.5)

# Graficar la trayectoria del Pistón B (Una línea recta)
plt.plot(x_B, y_B, 'r-', linewidth=3, label='Trayectoria Pistón (B)')

# DIBUJAR EL MECANISMO EN UNA POSICIÓN ESPECÍFICA (ej. índice 25 de 100)
pos = 25 
plt.plot([0, x_A[pos]], [0, y_A[pos]], 'g-o', linewidth=4, markersize=8, label='Manivela (L2)')
plt.plot([x_A[pos], x_B[pos]], [y_A[pos], y_B[pos]], 'm-o', linewidth=4, markersize=8, label='Biela (L3)')
plt.plot(x_B[pos], y_B[pos], 'ks', markersize=15, label='Pistón')

# Formato del gráfico
plt.title('Análisis de Posición: Biela-Manivela-Corredera', fontsize=14)
plt.xlabel('Posición X (cm)')
plt.ylabel('Posición Y (cm)')
plt.plot(0, 0, 'kx', markersize=10, label='Bancada (Origen)') # Origen
plt.axis('equal') # Para que el círculo no parezca un óvalo
plt.grid(True)
plt.legend(loc='upper right')

plt.show()