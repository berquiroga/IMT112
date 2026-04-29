import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. PARÁMETROS DEL SISTEMA
# ==========================================
# Geometría (metros)
L2, L3 = 0.05, 0.12  
# Masas (kg) e Inercias (kg*m^2)
m2, m3, m4 = 1.0, 2.0, 1.5 
I2, I3 = 0.01, 0.05
g = 9.81
F_gas = 1500.0 # Fuerza de compresión constante (N)

# Condiciones de Operación
RPM = 300
w2 = RPM * (2 * np.pi / 60) # rad/s
alpha2 = 0.0 # Motor a velocidad constante

# Arrays para guardar resultados y graficar
angulos_motor = np.arange(0, 360, 1) # De 0 a 359 grados
torque_motor = []

# ==========================================
# 2. EL BUCLE DE 360 GRADOS
# ==========================================
for angulo in angulos_motor:
    th2 = np.radians(angulo)
    
    # --- A. CINEMÁTICA (La caja negra resuelta) ---
    # Ángulos y posiciones
    th3 = np.arcsin((-L2 * np.sin(th2)) / L3)
    # Velocidades
    w3 = (-L2 * w2 * np.cos(th2)) / (L3 * np.cos(th3))
    # Aceleraciones
    alpha3 = (L2 * w2**2 * np.sin(th2) + L3 * w3**2 * np.sin(th3)) / (L3 * np.cos(th3))
    apx = -L2 * w2**2 * np.cos(th2) - L3 * alpha3 * np.sin(th3) - L3 * w3**2 * np.cos(th3)
    
    # Aceleraciones de los Centros de Gravedad (Asumiendo CG a la mitad de las barras)
    aG2x = -(L2/2) * w2**2 * np.cos(th2)
    aG2y = -(L2/2) * w2**2 * np.sin(th2)
    
    aG3x = -L2 * w2**2 * np.cos(th2) - (L3/2) * alpha3 * np.sin(th3) - (L3/2) * w3**2 * np.cos(th3)
    aG3y = -L2 * w2**2 * np.sin(th2) + (L3/2) * alpha3 * np.cos(th3) - (L3/2) * w3**2 * np.sin(th3)
    
    # Distancias a los pasadores (R) desde los CG
    R12x, R12y = -(L2/2)*np.cos(th2), -(L2/2)*np.sin(th2)
    R32x, R32y =  (L2/2)*np.cos(th2),  (L2/2)*np.sin(th2)
    R23x, R23y = -(L3/2)*np.cos(th3), -(L3/2)*np.sin(th3)
    R43x, R43y =  (L3/2)*np.cos(th3),  (L3/2)*np.sin(th3)

    # --- B. LA SÚPER MATRIZ DINÁMICA DE D'ALEMBERT ---
    A = np.zeros((8, 8))
    B = np.zeros(8)
    # Orden de Incógnitas: [F12x, F12y, F23x, F23y, F34x, F34y, F14N, Tm]

    # Eslabón 4 (Pistón)
    A[0, 4] = 1; B[0] = F_gas + m4 * apx
    A[1, 5] = 1; A[1, 6] = 1; B[1] = m4 * g

    # Eslabón 3 (Biela)
    A[2, 2] = 1; A[2, 4] = -1; B[2] = m3 * aG3x
    A[3, 3] = 1; A[3, 5] = -1; B[3] = m3 * aG3y + m3 * g
    A[4, 2] = -R23y; A[4, 3] = R23x; A[4, 4] = R43y; A[4, 5] = -R43x; B[4] = I3 * alpha3

    # Eslabón 2 (Manivela)
    A[5, 0] = 1; A[5, 2] = -1; B[5] = m2 * aG2x
    A[6, 1] = 1; A[6, 3] = -1; B[6] = m2 * aG2y + m2 * g
    A[7, 0] = -R12y; A[7, 1] = R12x; A[7, 2] = R32y; A[7, 3] = -R32x; A[7, 7] = 1; B[7] = I2 * alpha2

    # --- C. RESOLVER ---
    X = np.linalg.solve(A, B)
    torque_motor.append(X[7]) # Guardamos el Tm (posición 7)

# ==========================================
# 3. RESULTADOS Y GRÁFICA
# ==========================================
torque_pico = max(np.abs(torque_motor))

print(f"TORQUE PICO REQUERIDO: {torque_pico:.2f} N*m")

plt.plot(angulos_motor, torque_motor, color='red', linewidth=2)
plt.title("Torque Dinámico Requerido vs Ángulo del Motor")
plt.xlabel("Ángulo del Motor (Grados)")
plt.ylabel("Torque (N*m)")
plt.grid(True)
plt.axhline(0, color='black', linewidth=1)
plt.show()