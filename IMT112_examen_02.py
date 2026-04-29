import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# 1. Parámetros del sistema
# Longitudes de los eslabones (m)
L1, L2, L4 = 0.80, 0.20, 0.60
# Datos de masas e inercias
m2 = 5 # kg
m4 = 10 # kg
I2 = 0.02 # kg m2
I4 = 0.30 # kg m2

# Eslabón 3 - Trianguar
m3 = 35 # kg
I3 = 2.80 # kg m2
L3 = 0.90 # m
x3_cg = 0.45 # m
y3_cg = 0.133 # m
x3_c = 0.45 # m
y3_c = 0.40 # m

# Gravedad
g = 9.81 # m/s2

# Condiciones de trabajo
alpha2 = 0 # Aceleracion 0
RPM = 120 # Velocidad del motor constante
omega2 = RPM * (2 * np.pi / 60)

# Ecuaciones de lazo vectorial del mecanismo
# r1 + r4 = r2 + r3
def obtener_ecuaciones_sistema(incognitas, theta2_actual):
    theta3, theta4 = incognitas
    eq_x = L1 - L2 * np.cos(theta2_actual) - L3 * np.cos(theta3) + L4 * np.cos(theta4)  
    eq_y = - L2 * np.sin(theta2_actual) - L3 * np.sin(theta3) + L4 * np.sin(theta4)
    return [eq_x, eq_y]

# Ciclo principal
def simular_condicion(x3_cg, y3_cg):
    # Datos del eslabón 3
    H = np.sqrt(x3_cg**2 + y3_cg**2) # m
    beta = np.arctan2(y3_cg, x3_cg) # rad
    H_c = np.sqrt(x3_c**2 + y3_c**2) # m
    beta_c = np.arctan2(y3_c, x3_c) # rad

    # Usar aproximación inicial para ver si vértice C está bajando en ángulo theta2=359°
    angulo_anterior = 359
    theta2 = np.radians(angulo_anterior) # En radianes para las operaciones
    # Usar una aproximación inicial
    aproximacion_inicial = [np.radians(10), np.radians(120)]

    # Resolver numéricamente el sistema
    solucion_pos = fsolve(obtener_ecuaciones_sistema, aproximacion_inicial, args=(theta2,))
    # Convertimos de radianes a grados para que los humanos lo entiendan
    theta3 = solucion_pos[0] # en radianes
    theta4 = solucion_pos[1] # en radianes
    # Calculamos C en y3 en este ángulo
    c3y_anterior = L2 * np.sin(theta2) + H_c * np.sin(theta3 + beta_c)

    # Un arreglo para pasar por todos los grados y verificar los torques por grado
    torques_motor = []
    fuerzas_A = []
    fuerzas_externas = []
    alturas_c = []
    angulos_motor = np.arange(0, 360, 1) # Estos son los ángulos theta2
        
    for angulo in angulos_motor:
        # Actualizar el ángulo theta 2
        theta2 = np.radians(angulo) # En radianes para las operaciones

        # Usar una aproximación inicial
        aproximacion_inicial = [np.radians(10), np.radians(120)]

        # Resolver numéricamente el sistema
        solucion_pos = fsolve(obtener_ecuaciones_sistema, aproximacion_inicial, args=(theta2,))

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
        a3x = - L2 * np.cos(theta2) * omega2**2 - L2 * np.sin(theta2) * alpha2 - H * np.cos(theta3 + beta) * omega3**2 - H * np.sin(theta3 + beta) * alpha3
        a3y = - L2 * np.sin(theta2) * omega2**2 + L2 * np.cos(theta2) * alpha2 - H * np.sin(theta3 + beta) * omega3**2 + H * np.cos(theta3 + beta) * alpha3

        # Eslabón 4
        a4x = - L4/2 * np.cos(theta4) * omega4**2 - L4/2 * np.sin(theta4) * alpha4
        a4y = - L4/2 * np.sin(theta4) * omega4**2 + L4/2 * np.cos(theta4) * alpha4

        # Cálculo de coordenada C del eslabón 3
        c3x = L2 * np.cos(theta2) + H_c * np.cos(theta3 + beta_c)
        c3y = L2 * np.sin(theta2) + H_c * np.sin(theta3 + beta_c)
        alturas_c.append(c3y)

        # Si el punto C está bajando y se encuentra en un altura inferior a 20cm, hay fuerza en ese punto
        if c3y < c3y_anterior and c3y < 0.20:
            F_ext = 15000
        else:
            F_ext = 0
        fuerzas_externas.append(F_ext)
        c3y_anterior = c3y

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
        XC = XA + H_c * np.cos(theta3 + beta_c)
        YC = YA + H_c * np.sin(theta3 + beta_c)
        XO4, YO4 = L1, 0.0

        # Centros de Gravedad
        XG2, YG2 = XA / 2, YA / 2
        XG3 = XA + H * np.cos(theta3 + beta)
        YG3 = YA + H * np.sin(theta3 + beta)
        XG4 = XO4 + (L4/2) * np.cos(theta4)
        YG4 = YO4 + (L4/2) * np.sin(theta4)

        # Cálculo de vectores, brazos de palanca
        # Eslabón 2
        R12x, R12y = XO2 - XG2, YO2 - YG2
        R32x, R32y = XA - XG2, YA - YG2

        # Eslabón 3
        R23x, R23y = XA - XG3, YA - YG3
        R43x, R43y = XB - XG3, YB - YG3
        Rcx, Rcy = XC - XG3, YC - YG3

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
        B[4] = m3 * g + m3 * a3y + F_ext # La carga solamente actúa en un eslabón (no es acción-reacción)
        A[5, 2] = -R23y
        A[5, 3] = R23x
        A[5, 4] = R43y
        A[5, 5] = -R43x
        B[5] = I3 * alpha3 + Rcx * F_ext # La carga solamente actúa en un eslabón (no es acción-reacción)

        # Eslabón 4
        A[6, 4] = 1
        A[6, 6] = 1
        B[6] = m4 * a4x
        A[7, 5] = 1
        A[7, 7] = 1
        B[7] = m4 * g + m4 * a4y
        A[8, 4] = -R34y
        A[8, 5] = R34x
        A[8, 6] = -R14y
        A[8, 7] = R14x
        B[8] = I4 * alpha4

        # Resolver el sistema
        X = np.linalg.solve(A, B)
        torques_motor.append(X[8])
        fuerza_A = np.sqrt(X[2]**2 + X[3]**2)
        fuerzas_A.append(fuerza_A)
    return torques_motor, fuerzas_A, angulos_motor, alturas_c, fuerzas_externas

# Caso inicial
torques_base, f_A_base, angulos_motor, alt_c_base, f_ext_base = simular_condicion(x3_cg, y3_cg)
# Mostrar resultados finales
torque_pico_base = max(np.abs(torques_base))
fuerza_A_pico_base = max(np.abs(f_A_base))
print(f"TORQUE PICO REQUERIDO EN CASO BASE: {torque_pico_base:.2f} N*m")
print(f"FUERZA PICO EJERCIDA EN CASO BASE: {fuerza_A_pico_base:.2f} N")
print(f"CG base en: X'={x3_cg:.3f}m, Y'={y3_cg:.3f}m")

# ==========================================
# RUTINA DE OPTIMIZACIÓN DEL CG
# ==========================================
# Rango de -0.15m a +0.15m con pasos de 0.01m (1 cm)
x = np.arange(-0.15, 0.16, 0.01)
y = np.arange(-0.15, 0.16, 0.01)

# Queremos encontrar el MÍNIMO, así que inicializamos con un valor altísimo
torque_pico_minimo = float('inf') 
mejor_cg_x, mejor_cg_y = 0, 0

print("Iniciando optimización, por favor espere (evaluando 961 combinaciones)...")

for i in x:
    for j in y:
        if x3_cg + i > 0 and y3_cg + j > 0:
            torques, f_A, angulos_motor, alt_c, f_ext = simular_condicion(x3_cg + i, y3_cg + j)
            torque_pico_iteracion = max(np.abs(torques))
            
            # Si este diseño exige MENOS torque, lo guardamos como el nuevo mejor
            if torque_pico_iteracion < torque_pico_minimo:
                torque_pico_minimo = torque_pico_iteracion
                mejor_torques_motor = torques
                mejor_fuerzas_A = f_A
                mejor_alturas_c = alt_c
                mejor_fuerzas_externas = f_ext
                mejor_cg_x = x3_cg + i
                mejor_cg_y = y3_cg + j
            else:
                continue

print(f"Optimización completada. CG movido a: X'={mejor_cg_x:.3f}m, Y'={mejor_cg_y:.3f}m")

# Mostrar resultados finales
torque_pico = max(np.abs(mejor_torques_motor))
fuerza_A_pico = max(np.abs(mejor_fuerzas_A))

# Cálculo de factor de seguridad
factor_seguridad = 1.5
ssy_1045 = 310 * 10**6 # Pa
diametro = 0.018 # m
area = np.pi * diametro**2 / 4
tau = fuerza_A_pico / area
fs_1045 = ssy_1045 / tau

print(f"TORQUE PICO REQUERIDO EN MEJOR ITERACIÓN: {torque_pico:.2f} N*m")
print(f"FUERZA PICO EJERCIDA EN A EN MEJOR ITERACIÓN: {fuerza_A_pico:.2f} N")
print(f"FACTOR DE SEGURIDAD PARA LA ALEACIÓN ACERO SAE 1045 EN MEJOR ITERACIÓN: {fs_1045:.2f}")

# Gráficas
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
fig.suptitle("Auditoría Dinámica - Cizalla Industrial", fontsize=16, fontweight='bold')

# Gráfica para el torque del motor en cada ángulo
axes[0, 0].plot(angulos_motor, mejor_torques_motor, color='green', label="Diseño Optimizado", linewidth=2)
axes[0, 0].plot(angulos_motor, torques_base, color='red', label="Caso Base", linewidth=2, linestyle="--")
axes[0, 0].set_title("Torque Dinámico Requerido vs Ángulo del Motor")
axes[0, 0].set_xlabel("Ángulo del Motor (Grados)")
axes[0, 0].set_ylabel("Torque (N*m)")
axes[0, 0].grid(True)
axes[0, 0].legend()
axes[0, 0].axhline(0, color='black', linewidth=1)

# Gráfica para altura del pasador C
axes[0, 1].plot(angulos_motor, mejor_alturas_c, color='blue', label="Diseño Optimizado", linewidth=2)
axes[0, 1].plot(angulos_motor, alt_c_base, color='red', label="Caso Base", linewidth=2, linestyle="--")
axes[0, 1].set_title("Altura del punto C en el eslabón 3")
axes[0, 1].set_xlabel("Ángulo de Manivela (Grados)")
axes[0, 1].set_ylabel("Altura de C (m)")
axes[0, 1].grid(True)
axes[0, 1].legend()

# Gráfica para la fuerza total en el pasador A
axes[1, 0].plot(angulos_motor, mejor_fuerzas_A, color='blue', label="Diseño Optimizado", linewidth=2)
axes[1, 0].plot(angulos_motor, f_A_base, color='red', label="Caso Base", linewidth=2, linestyle="--")
axes[1, 0].set_title("Fuerza Cortante en Pasador A (F23)")
axes[1, 0].set_xlabel("Ángulo de Manivela (Grados)")
axes[1, 0].set_ylabel("Fuerza Absoluta (N)")
axes[1, 0].grid(True)
axes[1, 0].legend()

# Gráfica para las fuerzas externas a lo largo del movimiento
axes[1, 1].plot(angulos_motor, mejor_fuerzas_externas, color='purple', label="Diseño Optimizado", linewidth=2)
axes[1, 1].plot(angulos_motor, f_ext_base, color='red', label="Caso Base", linewidth=2, linestyle="--")
axes[1, 1].set_title("Perfil de Resistencia Generada por el Aluminio")
axes[1, 1].set_xlabel("Ángulo de Manivela (Grados)")
axes[1, 1].set_ylabel("Fuerza Externa (N)")
axes[1, 1].grid(True)
axes[1, 1].legend()
axes[1, 1].fill_between(angulos_motor, 0, mejor_fuerzas_externas, color='purple', alpha=0.2)

plt.tight_layout()
fig.subplots_adjust(top=0.92) # Da espacio al título principal
plt.show()