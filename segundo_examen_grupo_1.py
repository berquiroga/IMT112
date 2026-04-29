import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# 1. Parámetros del sistema
# Longitudes de los eslabones (m)
L1, L2, L3, L4 = 0.80, 0.20, 0.90, 0.60

# Valores de masas e inercias
m2 = 5 #kg
m3 = 35 #kg
m4 = 10 #kg
I2 = 0.02 #kg * m**2
I3 = 2.80 #kg * m**2
I4 = 0.30 #kg * m**2
# Gravedad
g = 9.81 # m/s2

# Coordenadas locales (respecto a A)
xG3p, yG3p = 0.45, 0.133
xCp, yCp = 0.45, 0.40

# Para detectar descenso
YC_anterior = None

# Condiciones de trabajo
alpha2 = 0 # Aceleracion 0
RPM = 120 # Velocidad del motor constante
omega2 = RPM * (2 * np.pi / 60)

# Un arreglo para pasar por todos los grados y verificar los torques por grado
torques_motor = []
fuerzas_A = []
omegas_4 = []
fuerzas_externas = []
angulos_motor = np.arange(0, 360, 1) # Estos son los ángulos theta2

# Usar una aproximación inicial
aproximacion_inicial = [np.radians(10), np.radians(120)]

# Ecuaciones de lazo vectorial del mecanismo
# r1 + r4 = r2 + r3
def obtener_ecuaciones_sistema(incognitas):
    theta3, theta4 = incognitas
    # Ecuacción en X: L2*cos(theta2) + L3*cos(theta3) = L1 + L4*cos(theta4)
    eq_x = L1 - L2 * np.cos(theta2) - L3 * np.cos(theta3) + L4 * np.cos(theta4)  
    # Ecuacción en X: L2*sin(theta2) + L3*sin(theta3) = L4*sin(theta4)
    eq_y = - L2 * np.sin(theta2) - L3 * np.sin(theta3) + L4 * np.sin(theta4)
    return [eq_x, eq_y]


for angulo in angulos_motor:
    # Actualizar el ángulo theta 2
    theta2 = np.radians(angulo) # En radianes para las operaciones

    # Resolver numéricamente el sistema
    solucion_pos = fsolve(obtener_ecuaciones_sistema, aproximacion_inicial)

    # Convertimos de radianes a grados para que los humanos lo entiendan
    theta3 = solucion_pos[0] # en radianes
    theta4 = solucion_pos[1] # en radianes

    # Actualizar la aproximación inicial
    aproximacion_inicial = [theta3, theta4]

    # Resolver para omega3 y omega4, un sistema 2x2 en forma matricial A*X=B
    # X = [omega3, omega4]
    # Matriz A: Coeficientes de la geometría
    A = np.array([
        [-L3 * np.sin(theta3),  L4 * np.sin(theta4)],  # Eje X
        [ L3 * np.cos(theta3), -L4 * np.cos(theta4)]   # Eje Y
    ])

    # Vector B: Los datos conocidos
    B = np.array([
        L2 * omega2 * np.sin(theta2),
        -L2 * omega2 * np.cos(theta2)
    ])

    # Resolución del sistema para las velocidades
    solucion_vel = np.linalg.solve(A, B)
    omega3, omega4 = solucion_vel
    omegas_4.append(omega4)

    # Resolver para alpha3 y alpha4, un sistema 2x2 en forma matricial A*X=B
    # Matriz A: Coeficientes de las incógnitas (Igual que en velocidad)
    # Despejamos alpha3 y alpha4 al lado izquierdo de la ecuación
    A = np.array([
        [-L3 * np.sin(theta3),  L4 * np.sin(theta4)],  # Fila del Eje X
        [ L3 * np.cos(theta3), -L4 * np.cos(theta4)]   # Fila del Eje Y
    ])

    # Vector B: Todos los términos conocidos (Las aceleraciones centrípetas/normales)
    # Mandamos todo lo que tiene w^2 y alpha2 al lado derecho
    B_x = ( L2 * alpha2 * np.sin(theta2) + L2 * omega2**2 * np.cos(theta2) + 
            L3 * omega3**2 * np.cos(theta3)  - L4 * omega4**2 * np.cos(theta4) )

    B_y = (-L2 * alpha2 * np.cos(theta2) + L2 * omega2**2 * np.sin(theta2) + 
            L3 * omega3**2 * np.sin(theta3)  - L4 * omega4**2 * np.sin(theta4) )

    B = np.array([B_x, B_y])

    # Resolución del sistema para las aceleraciones
    solucion_acel = np.linalg.solve(A, B)
    alpha3, alpha4 = solucion_acel

    # Cálculo de aceleraciones en centros de gravedad para cada eslabón
    # Eslabón 2
    a2x = - L2/2 * np.cos(theta2) * omega2**2 - L2/2 * np.sin(theta2) * alpha2
    a2y = - L2/2 * np.sin(theta2) * omega2**2 + L2/2 * np.cos(theta2) * alpha2

    # Eslabón 3
    a3x = (-L2*np.cos(theta2)*omega2**2 - xG3p*np.cos(theta3)*omega3**2 - yG3p*np.sin(theta3)*omega3**2 - xG3p*np.sin(theta3)*alpha3 + yG3p*np.cos(theta3)*alpha3)
    a3y = (-L2*np.sin(theta2)*omega2**2 - xG3p*np.sin(theta3)*omega3**2 + yG3p*np.cos(theta3)*omega3**2 + xG3p*np.cos(theta3)*alpha3 + yG3p*np.sin(theta3)*alpha3)

    # Eslabón 4
    a4x = - L4/2 * np.cos(theta4) * omega4**2 - L4/2 * np.sin(theta4) * alpha4
    a4y = - L4/2 * np.sin(theta4) * omega4**2 + L4/2 * np.cos(theta4) * alpha4

    # --- Matriz dinámica de D'Alembert ---
    A = np.zeros((9, 9))
    B = np.zeros(9)
    # Orden de Incógnitas: [0:F12x, 1:F12y, 2:F23x, 3:F23y, 4:F34x, 5:F34y, 6:F14x, 7:F14y, 8:Tm]

    # Determinar posiciones de centros de gravedad, pasadores y uniones a la bancada
    # Origen O2 en (0,0) y pivote O4 en (L1, 0)

    # Pasadores
    XO2, YO2 = 0.0, 0.0
    XA = L2 * np.cos(theta2)
    YA = L2 * np.sin(theta2)
    XB = L1 + L4 * np.cos(theta4)
    YB = L4 * np.sin(theta4)
    XO4, YO4 = L1, 0.0

    # Centros de Gravedad (Asumiendo a la mitad de las barras)
    XG2, YG2 = XA / 2, YA / 2
    XG3 = XA + xG3p*np.cos(theta3) - yG3p*np.sin(theta3)
    YG3 = YA + xG3p*np.sin(theta3) + yG3p*np.cos(theta3)
    XC = XA + xCp*np.cos(theta3) - yCp*np.sin(theta3)
    YC = YA + xCp*np.sin(theta3) + yCp*np.cos(theta3)
    XG4 = XO4 + (L4/2) * np.cos(theta4)
    YG4 = YO4 + (L4/2) * np.sin(theta4)

    # Analisis de Fuerza descendente de la cuchilla
    if YC_anterior is not None:
        descendiendo = YC < YC_anterior
    else:
        descendiendo = False

    # Aplicar condición
    if YC <= 0.20 and descendiendo:
        F_corte = 15000
    else:
        F_corte = 0

    # Recién aquí actualizas
    YC_anterior = YC

    
    fuerzas_externas.append(F_corte)

    # Cálculo de vectores, brazos de palanca
    # Eslabón 2
    R12x, R12y = XO2 - XG2, YO2 - YG2
    R32x, R32y = XA - XG2, YA - YG2

    # Eslabón 3
    R23x, R23y = XA - XG3, YA - YG3
    R43x, R43y = XB - XG3, YB - YG3

    RCx = XC - XG3
    RCy = YC - YG3

    # Eslabón 4
    R34x, R34y = XB - XG4, YB - YG4
    R14x, R14y = XO4 - XG4, YO4 - YG4

    # Eslabón 2
    A[0, 0] = 1
    A[0, 2] = -1
    B[0] = m2 * a2x
    A[1, 1] = 1
    A[1, 3] = -1
    B[1] = m2 * g + m2 * a2y
    A[2, 0] = -R12y
    A[2, 1] = R12x
    A[2, 2] = R32y
    A[2, 3] = -R32x
    A[2, 8] = 1
    B[2] = I2 * alpha2

    # Eslabón 3
    A[3, 2] = 1
    A[3, 4] = -1
    B[3] = m3 * a3x
    A[4, 3] = 1
    A[4, 5] = -1
    B[4] = m3 * g + m3 * a3y - F_corte
    A[5, 2] = -R23y
    A[5, 3] = R23x
    A[5, 4] = R43y
    A[5, 5] = -R43x
    B[5] = I3 * alpha3 - RCx * F_corte

    # Eslabón 4
    A[6, 4] = 1
    A[6, 6] = 1
    B[6] = m4 * a4x
    A[7, 5] = 1
    A[7, 7] = 1
    B[7] = m4 * g + m4 * a4y # La carga solamente actúa en un eslabón (no es acción-reacción)
    A[8, 4] = -R34y
    A[8, 5] = R34x
    A[8, 6] = -R14y
    A[8, 7] = R14x
    B[8] = I4 * alpha4 # La carga solamente actúa en un eslabón (no es acción-reacción)

    # Resolver el sistema
    X = np.linalg.solve(A, B)
    torques_motor.append(X[8])
    fuerza_A = np.sqrt(X[2]**2 + X[3]**2)
    fuerzas_A.append(fuerza_A)

# Mostrar resultados y graficar
torque_pico = max(np.abs(torques_motor))
fuerza_A_pico = max(np.abs(fuerzas_A))

# Cálculo de factor de seguridad
factor_seguridad = 1.5
ssy_SAE1045 = 310 * 10**6 # Pa
diametro = 0.018 # m
area = np.pi * diametro**2 / 4
tau = fuerza_A_pico / area
fs_SAE1045 = ssy_SAE1045 / tau

if fs_SAE1045 >= 1.5:
    print("DISEÑO SEGURO")
else:
    print("EL PASADOR FALLARA")


print(f"TORQUE PICO REQUERIDO: {torque_pico:.2f} N*m")
print(f"FUERZA PICO EJERCIDA EN A: {fuerza_A_pico:.2f} N")
print(f"FACTOR DE SEGURIDAD PARA LA ALEACIÓN ACERO ESTRUCTURAL SAE1045: {fs_SAE1045:.2f}")

# Gráficas
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
fig.suptitle("Auditoría Dinámica - Cizalla Automatizada - Análisis Dinámico", fontsize=16, fontweight='bold')

# Gráfica para el torque del motor en cada ángulo
axes[0, 0].plot(angulos_motor, torques_motor, color='red', linewidth=2)
axes[0, 0].set_title("Torque Dinámico Requerido vs Ángulo del Motor")
axes[0, 0].set_xlabel("Ángulo del Motor (Grados)")
axes[0, 0].set_ylabel("Torque (N*m)")
axes[0, 0].grid(True)
axes[0, 0].axhline(0, color='black', linewidth=1)

# Gráfica para la fuerza total en el pasador A
axes[0, 1].plot(angulos_motor, fuerzas_A, color='blue', linewidth=2)
axes[0, 1].set_title("Fuerza en Pasador A (F23)")
axes[0, 1].set_xlabel("Ángulo de Manivela (Grados)")
axes[0, 1].set_ylabel("Fuerza Absoluta (N)")
axes[0, 1].grid(True)

# Gráfica para la Velocidad Angular del Balancín (omega4)
axes[1, 0].plot(angulos_motor, omegas_4, color='green', linewidth=2)
axes[1, 0].set_title("Velocidad Angular del Balancín (\u03C94)") # \u03C9 es el símbolo omega
axes[1, 0].set_xlabel("Ángulo de Manivela (Grados)")
axes[1, 0].set_ylabel("Velocidad (rad/s)")
axes[1, 0].grid(True)
axes[1, 0].axhline(0, color='black', linewidth=1)
# Sombreado pedagógico: muestra la zona de avance
axes[1, 0].fill_between(angulos_motor, 0, omegas_4, where=(np.array(omegas_4) > 0), 
                        color='green', alpha=0.2, label='Carrera de Avance')
axes[1, 0].legend()

# Gráfica para las fuerzas externas a lo largo del movimiento
axes[1, 1].plot(angulos_motor, fuerzas_externas, color='purple', linewidth=2)
axes[1, 1].set_title("Perfil de Resistencia del Lodo")
axes[1, 1].set_xlabel("Ángulo de Manivela (Grados)")
axes[1, 1].set_ylabel("Fuerza Externa (N)")
axes[1, 1].grid(True)
axes[1, 1].fill_between(angulos_motor, 0, fuerzas_externas, color='purple', alpha=0.2)

plt.tight_layout()
fig.subplots_adjust(top=0.92) # Da espacio al título principal
plt.show()