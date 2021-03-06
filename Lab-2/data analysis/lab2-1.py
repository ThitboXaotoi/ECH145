"""
Credit to Toby Mea for a lot of stuffs in this code especially for the minimization of Numerical
Solution.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as so
import scipy.stats as ss
import scipy.integrate as si


# file_name = "/Users/Nguyen/Documents/ECH_145/Lab-2/data analysis/Lab2Data.xlsx"
file_name = "/Users/binhco/Documents/GitHub/ECH145/Lab-2/data analysis/Lab2Data.xlsx"

### Natural Convection
rawfile = pd.read_excel(file_name, "Natural Convection")
raw_time = np.array(rawfile)[2:, 1]
raw_time = np.array(raw_time, dtype=float)
raw_temps = np.array(rawfile)[2:, 3:5]
raw_temps = np.array(raw_temps, dtype=float)

raw_ambient_temp = raw_temps[:, 1]
raw_center_temp = raw_temps[:, 0]

i_start = np.argmax(raw_center_temp)
nat_time = np.array(raw_time[i_start:]) - np.array(raw_time[i_start])
nat_ambient_temp = np.array(raw_ambient_temp[i_start:])
nat_center_temp = np.array(raw_center_temp[i_start:])

nat_T_inf = np.mean(nat_ambient_temp)
nat_T_init = nat_center_temp[0]


### Forced Convection 1
rawfile = pd.read_excel(file_name, "Forced Convection 1")
raw_time = np.array(rawfile)[2:, 0]
raw_time = np.array(raw_time, dtype=float)
raw_temps = np.array(rawfile)[2:, 1:3]
raw_temps = np.array(raw_temps, dtype=float)

raw_ambient_temp = raw_temps[:, 1]
raw_center_temp = raw_temps[:, 0]

i_start = np.argmax(raw_center_temp)
for1_time = np.array(raw_time[i_start:]) - np.array(raw_time[i_start])
for1_ambient_temp = np.array(raw_ambient_temp[i_start:])
for1_center_temp = np.array(raw_center_temp[i_start:])

for1_T_inf = np.mean(for1_ambient_temp)
for1_T_init = for1_center_temp[0]


### Parameters
rho = 2700 #kg/m**3
diameter = 0.02536 #m
length = 0.06018

# diameter = 0.02537 #m
# length = 0.06019
volume = np.pi * diameter**2 * length / 4
sur_area = np.pi * diameter * length
mass = rho * volume

Cp = 900 #J/kgK
k = 170 #W/mK

g = 9.81 #m/s**2


#air
nu = 15.15E-6 #m**2/s
k_air = 25.95E-3 
thermal_expansion = 3.43E-3
Pr = 0.707


"""
---------------------
Empirical Correlation

"""
# Natural Convection

Gr = (g * thermal_expansion * (np.mean(nat_center_temp) - nat_T_inf) * diameter**3) / (nu**2)

def Nusselt(Nu):
    return Nu * np.exp(-2/Nu) - 0.6*((diameter/length)*Gr*Pr)**0.25

Nu = so.fsolve(Nusselt, 0.01)

h_nat_emp = Nu * k_air/diameter
h_nat_emp_err = h_nat_emp * ((0.00001/diameter) + (0.00001/diameter + 0.00001/length +
                                                    0.05/nat_T_inf + 0.00001/diameter)*0.25)


print("The Empirical Correlation gives an h of " + str(round(float(h_nat_emp), 2)) + " ± " +
      str(round(float(h_nat_emp_err), 2)))


# Forced Convection
velocity = 1.17 #m/s

Re = velocity * diameter / nu

if Re > 4:
    if Re > 40:
        if Re > 4000:
            if Re > 40000:
                if Re > 400000:
                    C = 0.027
                    m = 0.805
            else:
                C = 0.193
                m = 0.618
        else:
            C = 0.683
            m = 0.466
    else:
        C = 0.911
        m = 0.385
else:
    C = 0.989
    m = 0.330

h_empirical = k_air * C * (Re**m) * (Pr**(1/3)) / diameter
h_empirical_err = h_empirical * ((((0.01/velocity)+ (0.00001/diameter)) * m/(Re**m)) +
                                 (0.00001/diameter))

print("The Empirical Correlation gives an h of " + str(round(h_empirical, 2)) + " ± " +
      str(round(h_empirical_err, 2)))


"""
---------------------
Analytic Solution

"""


nat_LHS = np.log((nat_center_temp - nat_T_inf) / (nat_T_init - nat_T_inf))
plt.figure(0)
plt.plot(nat_time, nat_LHS, label = "Natural Convection")

for1_LHS = np.log((for1_center_temp - for1_T_inf) / (for1_T_init - for1_T_inf))
plt.plot(for1_time, for1_LHS, label = "Forced Convection")


def y(t, m):
    y = -t * m
    return y

nat_res = so.curve_fit(y, nat_time, nat_LHS)
nat_tau = float(1 / np.array(nat_res[0]))
for1_res = so.curve_fit(y, for1_time, for1_LHS)
for1_tau = float(1 / np.array(for1_res[0]))
 
plt.plot(nat_time, y(nat_time, np.array(nat_res[0])), c = 'k', label='Natural Tau = ' +
         str(round(nat_tau, 2)) + "$s^{-1}$") 
plt.plot(for1_time, y(for1_time, np.array(for1_res[0])), c = 'k', label='Forced Tau = ' +
         str(round(for1_tau, 2)) + "$s^{-1}$") 
plt.xlabel('Time (s)')
plt.ylabel('Log of Dimensionless Temperature')
plt.legend()

analytical_nat_h = mass * Cp / (nat_tau * sur_area)
analytical_for1_h = mass * Cp / (for1_tau * sur_area)

nat_delta_tau = 0.00008054 / nat_res[0]
for1_delta_tau = 0.00033873 / for1_res[0]

analytical_nat_h_err = float(analytical_nat_h * (nat_delta_tau + (0.00001/diameter) +
                                                 (0.00001/length)))
analytical_for1_h_err = float(analytical_for1_h * (for1_delta_tau + (0.00001/diameter) +
                                                 (0.00001/length)))

print("Analytical Solution for Natural Convection is " + str(round(analytical_nat_h, 2)) + " ± " +
      str(round(analytical_nat_h_err, 2)))

print("Analytical Solution for Forced Convection is " + str(round(analytical_for1_h, 2)) + " ± " +
      str(round(analytical_for1_h_err, 2)))



"""
---------------------
Numerical Solution

"""


#Forced Convection
alfa = k / (rho * Cp)
dt = 0.05
dr = diameter/ (2 * 4)
beta = mass * Cp / dt
s = alfa * dt / dr**2
print("s is equal to " + str(s))
if (s>0.5):
    print("bad bad change your dr dt")


num_for1_time = np.arange(0, int(for1_time[-1]), dt)

for1_T = np.zeros((len(num_for1_time), int(diameter/(2 * dr))))
for1_T[0, :] = for1_center_temp[0]

def F_sq_diff(time):
    exp_idx = np.where(for1_time == time)[0]
    idx = np.where(num_for1_time == time)[0]
    numerical = for1_T[idx, -1]
    experimental = for1_center_temp[int(exp_idx)]
    diff = numerical - experimental
    return diff**2

def F_num_method_first(h):
    for m in range(1, len(num_for1_time)):
        for i in range(1, int(diameter / (2 * dr)) - 1):
            a_i = s - s/(2*i)
            b_i = 1 - 2*s
            c_i = s + s/(2*i)
            for1_T[m, i] = a_i * for1_T[m-1, i-1] + b_i * for1_T[m-1, i] + c_i * for1_T[m-1, i+1]
            for1_T[m, 0] = for1_T[m, 1]
            for1_T[m, -1] = ((beta * for1_T[m-1, -1]) + (h * sur_area * for1_T_inf)) / (beta + (h *
                                                                                               sur_area))           
    print("Completed Solution")
    residuals = []
    for i in range(len(for1_time)):
        res = F_sq_diff(for1_time[i])
        residuals.append(res)
    residuals = np.array(residuals[:-1], dtype = float)
    return float(sum(residuals))

# guess = float(input("Guess the h value for forced convection "))
guess = 27

for1_ans = so.least_squares(F_num_method_first, guess, bounds = [guess-2, guess+2])
print("Numerical Solution for Forced Convection is " + str(round(float(for1_ans.x), 2)) + " ± " +
      str(0.0110))

numerical_for1_h = float(for1_ans.x)

def F_num_method_final(h):
    for m in range(1, len(num_for1_time)):
        for i in range(1, int(diameter/(2 * dr))-1):
            a_i = s - s/(2*i)
            b_i = 1 - 2*s
            c_i = s + s/(2*i)
            for1_T[m, i] = a_i * for1_T[m-1, i-1] + b_i * for1_T[m-1, i] + c_i * for1_T[m-1, i+1]
            for1_T[m, 0] = for1_T[m, 1]
            for1_T[m, -1] = ((beta * for1_T[m-1, -1]) + (h * sur_area * for1_T_inf)) / (beta + (h *
                                                                                               sur_area))

    print("Completed Solution")
    return(for1_T)

for1_T = F_num_method_final(for1_ans.x)

for1_heat_transfer_tot = numerical_for1_h * np.pi * diameter * length * si.simps(for1_T[:, -1] -
                                                                            for1_T_inf,
                                                                            num_for1_time)

print("Total Heat Transfer is " + str(round(for1_heat_transfer_tot, 2)) + " Joules")

plt.figure(1)
plt.plot(num_for1_time, for1_T[:, -1], 'k', ls='-', lw=2, alpha=0.5, label='Numerical Forced Convection Solution')


### Natural Convection
alfa = k / (rho * Cp)
dt = 0.05
dr = diameter/ (2 * 4)
beta = mass * Cp / dt
s = alfa * dt / dr**2
print("s is equal to " + str(s))
if (s>0.5):
    print("bad bad change your dr dt")

num_nat_time = np.arange(0, int(nat_time[-1]), dt)

nat_T = np.zeros((len(num_nat_time), int(diameter/(2 * dr))))
nat_T[0, :] = nat_center_temp[0]

def N_sq_diff(time):
    exp_idx = np.where(nat_time == time)[0]
    idx = np.where(num_nat_time == time)[0]
    numerical = nat_T[idx, -1]
    experimental = nat_center_temp[exp_idx]
    diff = numerical - experimental
    return diff**2

def N_num_method_first(h):
    for m in range(1, len(num_nat_time)):
        for i in range(1, int(diameter / (2 * dr)) - 1):
            a_i = s - s/(2*i)
            b_i = 1 - 2*s
            c_i = s + s/(2*i)
            nat_T[m, i] = a_i * nat_T[m-1, i-1] + b_i * nat_T[m-1, i] + c_i * nat_T[m-1, i+1]
            nat_T[m, 0] = nat_T[m, 1]
            nat_T[m, -1] = ((beta * nat_T[m-1, -1]) + (h * sur_area * nat_T_inf)) / (beta + (h *
                                                                                               sur_area))           
    print("Completed Solution")
    residuals = []
    for i in range(len(nat_time)):
        res = N_sq_diff(nat_time[i])
        residuals.append(res)
    residuals = np.array(residuals[:-1], dtype = float)
    return float(sum(residuals))

# guess = float(input("Guess the h value for Natural Convection: "))

guess = 11

num_ans = so.least_squares(N_num_method_first, guess, bounds = [guess-2, guess+2])
print("Numerical Solution for Natural Convection is " + str(round(float(num_ans.x), 2)) + " ± " +
      str(0.0047))
numerical_nat_h = float(num_ans.x)

def N_num_method_final(h):
    for m in range(1, len(num_nat_time)):
        for i in range(1, int(diameter / (2 * dr)) - 1):
            a_i = s - s/(2*i)
            b_i = 1 - 2*s
            c_i = s + s/(2*i)
            nat_T[m, i] = a_i * nat_T[m-1, i-1] + b_i * nat_T[m-1, i] + c_i * nat_T[m-1, i+1]
            nat_T[m, 0] = nat_T[m, 1]
            nat_T[m, -1] = ((beta * nat_T[m-1, -1]) + (h * sur_area * nat_T_inf)) / (beta + (h *
                                                                                               sur_area))           
    # print("Completed Solution")
    return nat_T

nat_T = N_num_method_final(num_ans.x)

nat_heat_transfer_tot = numerical_nat_h * np.pi * diameter * length * si.simps(nat_T[:, -1] -
                                                                            nat_T_inf,
                                                                            num_nat_time)

print("Total Heat Transfer is " + str(round(nat_heat_transfer_tot, 2)) + " Joules")

def cooling_to_point_1(t):
    return (0.1 / (nat_T_init - nat_T_inf) - np.exp(-numerical_nat_h * np.pi * diameter * length * t /
                                                    (mass * Cp)))

cooling_time = float(so.fsolve(cooling_to_point_1, 0.01))

print(f"It takes {cooling_time} seconds to cool the cylinder to within 0.1 deg C")

plt.plot(num_nat_time, nat_T[:, -1], 'r', ls='-', lw=2, alpha=0.5, label='Numerical Natural Convection Solution')


step = 20
new_for1_time = []
new_for1_center_temp = []
for i in range(int(len(for1_time)/step)):
    ft = for1_time[i * step]
    new_for1_time.append(ft)
    ftemp = for1_center_temp[i * step]
    new_for1_center_temp.append(ftemp)



step = 20
new_nat_time = []
new_nat_center_temp = []
for i in range(int(len(nat_time)/step)):
    ft = nat_time[i * step]
    new_nat_time.append(ft)
    ftemp = nat_center_temp[i * step]
    new_nat_center_temp.append(ftemp)

plt.scatter(new_nat_time, new_nat_center_temp, s=20, alpha=1, c='c', label='Natural Convection Data')
plt.scatter(new_for1_time, new_for1_center_temp, s=20, alpha=1, label='Forced Convection Data')

plt.plot(nat_time, (nat_T_init - nat_T_inf) * np.exp(-(nat_time / nat_tau)) + nat_T_inf, c='r',
         ls=':', label='Natural Analytical Solution')
plt.plot(for1_time, (for1_T_init - for1_T_inf) * np.exp(-(for1_time / for1_tau)) + for1_T_inf, c='k',
         ls=':', label='Forced Analytical Solution')

plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (˚C)')

plt.show()
