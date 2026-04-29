import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. PARÁMETROS DEL MECANISMO (Biela-Manivela)
# ==========================================
# Definimos las longitudes físicas de los eslabones en centímetros
L2 = 5.0   # Longitud de la Manivela (Eslabón de entrada)
L3 = 12.0  # Longitud de la Biela (Eslabón acoplador)
offset = 0.0 # Desfase vertical de la línea de deslizamiento de la corredera

# Creamos un arreglo de 100 posiciones angulares para la manivela (de 0 a 360 grados)
theta2_vals = np.linspace(0, 2 * np.pi, 100)

# ==========================================
# 2. ANÁLISIS DE POSICIÓN (Vectores)
# ==========================================
# Coordenadas del Punto A (Pasador entre Manivela y Biela)
# Esto es simplemente descomponer el vector L2 en X e Y
X_A = L2 * np.cos(theta2_vals)
Y_A = L2 * np.sin(theta2_vals)

# Coordenadas del Punto B (Pasador de la Corredera / Pistón)
# Sabemos que la distancia entre A y B siempre es L3.
# Aplicamos el Teorema de Pitágoras: (X_B - X_A)^2 + (Y_B - Y_A)^2 = L3^2
# Como Y_B es fijo (igual al offset), despejamos X_B:
X_B = X_A + np.sqrt(L3**2 - (Y_A - offset)**2)
Y_B = np.full_like(X_B, offset) # La corredera siempre se mueve sobre la línea Y = offset

# ==========================================
# 3. GRAFICACIÓN DEL MECANISMO
# ==========================================
# Configuramos el tamaño del lienzo
plt.figure(figsize=(10, 5))

# 3.1 Graficar las trayectorias completas (El "camino" de los pasadores)
plt.plot(X_A, Y_A, 'b--', label='Trayectoria Punto A (Círculo)', alpha=0.5)
plt.plot(X_B, Y_B, 'r-', linewidth=3, label='Trayectoria Punto B (Línea Recta)')

# 3.2 Dibujar los eslabones en una posición específica (ej. frame 20)
frame = 80
# Dibujamos la Manivela (Línea desde el Origen hasta el Punto A)
plt.plot([0, X_A[frame]], [0, Y_A[frame]], 'g-o', linewidth=4, markersize=8, label='Manivela (L2)')
# Dibujamos la Biela (Línea desde el Punto A hasta el Punto B)
plt.plot([X_A[frame], X_B[frame]], [Y_A[frame], Y_B[frame]], 'm-o', linewidth=4, markersize=8, label='Biela (L3)')
# Dibujamos la Corredera (Un cuadrado negro en el Punto B)
plt.plot(X_B[frame], Y_B[frame], 'ks', markersize=15, label='Corredera')

# 3.3 Detalles estéticos del gráfico
plt.plot(0, 0, 'kx', markersize=10, label='Bancada (O12)') # Marca el origen
plt.title('Análisis de Posición: Biela-Manivela-Corredera', fontsize=14, fontweight='bold')
plt.xlabel('Posición X (cm)')
plt.ylabel('Posición Y (cm)')
plt.axis('equal') # IMPORTANTE: Evita que el círculo de la manivela parezca un óvalo
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right')

# Mostrar el gráfico en pantalla
plt.show()