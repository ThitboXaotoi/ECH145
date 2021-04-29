import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# material properties
k = 167 #W/mK
density = 2770 #kg/m**3
Cp = 896 #J/KgK

alfa = k / (density * Cp)


#geometry
radius = 0.025 / 2 #m
height = 0.0618

A_s = 2 * np.pi * radius * height #m**2

volume = np.pi * radius**2 * height #m**3

mass = density * volume #kg


file_name = "/Users/Nguyen/Documents/ECH_145/Lab-2/data analysis/Lab2Data.xlsx"
raw_dat = pd.read_excel(file_name, sheet_name=1)
raw_dat = np.array(raw_dat)

exp_time = raw_dat[1:, 0]
cylinder_T = raw_dat[1:, 1]
ambient_T = raw_dat[1:, 2]


time = 100 #s

# initial condition
ic = cylinder_T.max() #deg C

T_ambient = np.mean(ambient_T)

dr = 0.00125
dt = 0.01

x_steps = np.arange(0, radius, step = dr)
t_steps = np.arange(0, time, step = dt)

n_x = len(x_steps)
n_t = len(t_steps)

s = alfa * dt / dr**2

if (s>0.5):
    print("bad bad change your dr dt")

h = 11 #W /(m**2*K)

BC2 = 0 #deg C





T = np.zeros((n_t, n_x))

T[0] = ic

for m in range(1, n_t):
    for i in range(1, n_x-1):
        a_i = s - s/(2*i)
        b_i = 1 - 2*s
        c_i = s + s/(2*i)
        T[m, i] = a_i * T[m-1, i-1] + b_i * T[m-1, i] + c_i * T[m-1, i+1]

        T[m, 0] = T[m, 1] #symmetry boundary condition

        #BC2
        # H = -h * A_s * dt /(mass * Cp)
        # T[m, -1] = H * (T[m-1, -1] - T_ambient) * T[m-1, -1] 
        beta = mass * Cp/dt
        T[m, -1] = ((beta + T[m-1, 1]) + (h * A_s * T_ambient)) / (beta + (h * A_s))   


# plt.plot(x_steps, T[0])
# plt.plot(x_steps, T[50])

# plt.plot(t_steps, T[:, -1])
# plt.plot(t_steps, T[:, 0])


plt.figure(0)
plt.plot(x_steps, T[0])
plt.plot(x_steps, T[50])

plt.figure(1)
# plt.plot(exp_time, cylinder_T)
plt.plot(t_steps, T[:, 0])


