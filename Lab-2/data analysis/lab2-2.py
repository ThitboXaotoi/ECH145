import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as so
import scipy.stats as ss


# file_name = "/Users/Nguyen/Documents/ECH_145/Lab-2/data analysis/Lab2Data.xlsx"
file_name = "/Users/binhco/Documents/GitHub/ECH145/Lab-2/data analysis/Lab2Data.xlsx"

### Forced Convection 3
rawfile = pd.read_excel(file_name, "Forced Convection 2")
raw_time = np.array(rawfile)[2:, 4]
raw_time = np.array(raw_time, dtype=float)
raw_temps = np.array(rawfile)[2:, 2:4]
raw_temps = np.array(raw_temps, dtype=float)

raw_ambient_temp = raw_temps[:, 1]
raw_center_temp = raw_temps[:, 0]

i_start = np.argmax(raw_center_temp)
for2_time = np.array(raw_time[i_start:]) - np.array(raw_time[i_start])
for2_ambient_temp = np.array(raw_ambient_temp[i_start:])
for2_center_temp = np.array(raw_center_temp[i_start:])

for2_T_inf = np.mean(for2_ambient_temp)
for2_T_init = for2_center_temp[0]

### Forced Convection 3
rawfile = pd.read_excel(file_name, "Forced Convection backup 1")
raw_time = np.array(rawfile)[2:, 2]
raw_time = np.array(raw_time, dtype=float)
raw_temps = np.array(rawfile)[2:, 3:5]
raw_temps = np.array(raw_temps, dtype=float)

raw_ambient_temp = raw_temps[:, 1]
raw_center_temp = raw_temps[:, 0]

i_start = np.argmax(raw_center_temp)
for3_time = np.array(raw_time[i_start:]) - np.array(raw_time[i_start])
for3_ambient_temp = np.array(raw_ambient_temp[i_start:])
for3_center_temp = np.array(raw_center_temp[i_start:])

for3_T_inf = np.mean(for3_ambient_temp)
for3_T_init = for3_center_temp[0]

### Forced Convection 4
rawfile = pd.read_excel(file_name, "Forced Convection backup 2")
raw_time = np.array(rawfile)[2:, 0]
raw_time = np.array(raw_time, dtype=float)
raw_temps = np.array(rawfile)[2:, 1:3]
raw_temps = np.array(raw_temps, dtype=float)

raw_ambient_temp = raw_temps[:, 1]
raw_center_temp = raw_temps[:, 0]

i_start = np.argmax(raw_center_temp)
for4_time = np.array(raw_time[i_start:]) - np.array(raw_time[i_start])
for4_ambient_temp = np.array(raw_ambient_temp[i_start:])
for4_center_temp = np.array(raw_center_temp[i_start:])

for4_T_inf = np.mean(for4_ambient_temp)
for4_T_init = for4_center_temp[0]

### Parameters
rho = 2700 #kg/m**3
diameter = 0.02536 #m
length = 0.06018

volume = np.pi * diameter**2 * length / 4
sur_area = np.pi * diameter * length
mass = rho * volume

Cp = 900 #J/kgK
k = 170 #W/mK

#air
nu = 15.15E-6 #m**2/s
k_air = 25.95E-3 



###Empirical Correlation

velocity = np.array([1.88, 0.94, 1.78]) #m/s
Pr = 0.707

Re = velocity * diameter / nu

C = 0.683
m = 0.466

h_empirical = k_air * C * (Re**m) * (Pr**(1/3)) / diameter
h_empirical_err = h_empirical * ((((0.01/velocity) + (2 * 0.00001/diameter)) * m) + (0.00001/diameter))

for i in range(len(h_empirical)):
    print(f"The Empirical Correlation gives an h{i+2} of " + str(round(h_empirical[i], 2)) + " ± " +
          str(round(h_empirical_err[i], 2)))


###Analytic Solution

plt.figure(0)

for2_LHS = np.log((for2_center_temp - for2_T_inf) / (for2_T_init - for2_T_inf))
plt.plot(for2_time, for2_LHS, label = "Forced Convection 2")

for3_LHS = np.log((for3_center_temp - for3_T_inf) / (for3_T_init - for3_T_inf))
plt.plot(for3_time, for3_LHS, label = "Forced Convection 3")

for4_LHS = np.log((for4_center_temp - for4_T_inf) / (for4_T_init - for4_T_inf))
plt.plot(for4_time, for4_LHS, label = "Forced Convection 4")



def y(t, m):
    y = -t * m
    return y

for2_res = so.curve_fit(y, for2_time, for2_LHS)
for2_tau = float(1 / np.array(for2_res[0]))
 
for3_res = so.curve_fit(y, for3_time, for3_LHS)
for3_tau = float(1 / np.array(for3_res[0]))

for4_res = so.curve_fit(y, for4_time, for4_LHS)
for4_tau = float(1 / np.array(for4_res[0]))


plt.plot(for2_time, y(for2_time, np.array(for2_res[0])), label='Forced Tau = ' +
         str(round(for2_tau, 2)) + "$s^{-1}$") 
plt.plot(for3_time, y(for3_time, np.array(for3_res[0])), label='Forced Tau = ' +
         str(round(for3_tau, 2)) + "$s^{-1}$")
plt.plot(for4_time, y(for4_time, np.array(for4_res[0])), label='Forced Tau = ' +
         str(round(for4_tau, 2)) + "$s^{-1}$")

plt.xlabel('Time (s)')
plt.ylabel('Log of Dimensionless Temperature')
plt.legend()

analytical_for2_h = mass * Cp / (for2_tau * sur_area)
analytical_for3_h = mass * Cp / (for3_tau * sur_area)
analytical_for4_h = mass * Cp / (for4_tau * sur_area)

for2_delta_tau = 0.00033372 / for2_res[0]
for3_delta_tau = 0.00039077 / for3_res[0]
for4_delta_tau = 0.00027497 / for4_res[0]

analytical_for2_h_err = float(analytical_for2_h * (for2_delta_tau + (0.00001/diameter) +
                                                 (0.00001/length)))
analytical_for3_h_err = float(analytical_for3_h * (for3_delta_tau + (0.00001/diameter) +
                                                 (0.00001/length)))
analytical_for4_h_err = float(analytical_for4_h * (for4_delta_tau + (0.00001/diameter) +
                                                 (0.00001/length)))

print("Analytical Solution for Forced Convection is " + str(round(analytical_for2_h, 2)) + " ± " +
      str(round(analytical_for2_h_err, 2)))

print("Analytical Solution for Forced Convection is " + str(round(analytical_for3_h, 2)) + " ± " +
      str(round(analytical_for3_h_err, 2)))

print("Analytical Solution for Forced Convection is " + str(round(analytical_for4_h, 2)) + " ± " +
      str(round(analytical_for4_h_err, 2)))

###Numerical Solution

#Forced Convection 2
alfa = k / (rho * Cp)
dt = 0.05
dr = diameter/ (2 * 4)
beta = mass * Cp / dt
s = alfa * dt / dr**2
print("s is equal to " + str(s))
if (s>0.5):
    print("bad bad change your dr dt")


num_for2_time = np.arange(0, int(for2_time[-1]), dt)

for2_T = np.zeros((len(num_for2_time), int(diameter/(2 * dr))))
for2_T[0, :] = for2_center_temp[0]

def F2_sq_diff(time):
    exp_idx = np.where(for2_time == time)[0]
    idx = np.where(num_for2_time == time)[0]
    numerical = for2_T[idx, -1]
    experimental = for2_center_temp[int(exp_idx)]
    diff = numerical - experimental
    return diff**2

def F2_num_method_first(h):
    for m in range(1, len(num_for2_time)):
        for i in range(1, int(diameter / (2 * dr)) - 1):
            a_i = s - s/(2*i)
            b_i = 1 - 2*s
            c_i = s + s/(2*i)
            for2_T[m, i] = a_i * for2_T[m-1, i-1] + b_i * for2_T[m-1, i] + c_i * for2_T[m-1, i+1]
            for2_T[m, 0] = for2_T[m, 1]
            for2_T[m, -1] = ((beta * for2_T[m-1, -1]) + (h * sur_area * for2_T_inf)) / (beta + (h *
                                                                                               sur_area))           
    print("Completed Solution")
    residuals = []
    for i in range(len(for2_time)):
        res = F2_sq_diff(for2_time[i])
        residuals.append(res)
    residuals = np.array(residuals[:-1], dtype = float)
    return float(sum(residuals))

# guess = float(input("Guess the h value for forced convection "))
guess = 40

for2_ans = so.least_squares(F2_num_method_first, guess, bounds = [guess-2, guess+2])
print("Numerical Solution for Forced Convection is " + str(round(float(for2_ans.x), 2)) + " ± " +
      str(0.0163))

def F2_num_method_final(h):
    for m in range(1, len(num_for2_time)):
        for i in range(1, int(diameter/(2 * dr))-1):
            a_i = s - s/(2*i)
            b_i = 1 - 2*s
            c_i = s + s/(2*i)
            for2_T[m, i] = a_i * for2_T[m-1, i-1] + b_i * for2_T[m-1, i] + c_i * for2_T[m-1, i+1]
            for2_T[m, 0] = for2_T[m, 1]
            for2_T[m, -1] = ((beta * for2_T[m-1, -1]) + (h * sur_area * for2_T_inf)) / (beta + (h *
                                                                                               sur_area))

    print("Completed Solution")
    return(for2_T)

for2_T = F2_num_method_final(for2_ans.x)

plt.figure(1)
plt.plot(num_for2_time, for2_T[:, -1], 'b', ls='-', lw=2, alpha=0.5, label=f'Numerical Forced Convection (velocity = {velocity[0]}) Solution')


#Forced Convection 3
alfa = k / (rho * Cp)
dt = 0.05
dr = diameter/ (2 * 4)
beta = mass * Cp / dt
s = alfa * dt / dr**2
print("s is equal to " + str(s))
if (s>0.5):
    print("bad bad change your dr dt")


num_for3_time = np.arange(0, int(for3_time[-1]), dt)

for3_T = np.zeros((len(num_for3_time), int(diameter/(2 * dr))))
for3_T[0, :] = for3_center_temp[0]

def F3_sq_diff(time):
    exp_idx = np.where(for3_time == time)[0]
    idx = np.where(num_for3_time == time)[0]
    numerical = for3_T[idx, -1]
    experimental = for3_center_temp[int(exp_idx)]
    diff = numerical - experimental
    return diff**2

def F3_num_method_first(h):
    for m in range(1, len(num_for3_time)):
        for i in range(1, int(diameter / (2 * dr)) - 1):
            a_i = s - s/(2*i)
            b_i = 1 - 2*s
            c_i = s + s/(2*i)
            for3_T[m, i] = a_i * for3_T[m-1, i-1] + b_i * for3_T[m-1, i] + c_i * for3_T[m-1, i+1]
            for3_T[m, 0] = for3_T[m, 1]
            for3_T[m, -1] = ((beta * for3_T[m-1, -1]) + (h * sur_area * for3_T_inf)) / (beta + (h *
                                                                                               sur_area))           
    print("Completed Solution")
    residuals = []
    for i in range(len(for3_time)):
        res = F3_sq_diff(for3_time[i])
        residuals.append(res)
    residuals = np.array(residuals[:-1], dtype = float)
    return float(sum(residuals))

# guess = float(input("Guess the h value for forced convection "))
guess = 40

for3_ans = so.least_squares(F3_num_method_first, guess, bounds = [guess-2, guess+2])

print("Numerical Solution for Forced Convection is " + str(round(float(for3_ans.x), 2)) + " ± " +
      str(0.0161))


def F3_num_method_final(h):
    for m in range(1, len(num_for3_time)):
        for i in range(1, int(diameter/(2 * dr))-1):
            a_i = s - s/(2*i)
            b_i = 1 - 2*s
            c_i = s + s/(2*i)
            for3_T[m, i] = a_i * for3_T[m-1, i-1] + b_i * for3_T[m-1, i] + c_i * for3_T[m-1, i+1]
            for3_T[m, 0] = for3_T[m, 1]
            for3_T[m, -1] = ((beta * for3_T[m-1, -1]) + (h * sur_area * for3_T_inf)) / (beta + (h *
                                                                                               sur_area))

    print("Completed Solution")
    return(for3_T)

for3_T = F3_num_method_final(for3_ans.x)

plt.plot(num_for3_time, for3_T[:, -1], 'g', ls='-', lw=2, alpha=0.5, label=f'Numerical Forced Convection (velocity = {velocity[1]}) Solution')



#Forced Convection 4
alfa = k / (rho * Cp)
dt = 0.05
dr = diameter/ (2 * 4)
beta = mass * Cp / dt
s = alfa * dt / dr**2
print("s is equal to " + str(s))
if (s>0.5):
    print("bad bad change your dr dt")


num_for4_time = np.arange(0, int(for4_time[-1]), dt)

for4_T = np.zeros((len(num_for4_time), int(diameter/(2 * dr))))
for4_T[0, :] = for4_center_temp[0]

def F4_sq_diff(time):
    exp_idx = np.where(for4_time == time)[0]
    idx = np.where(num_for4_time == time)[0]
    numerical = for4_T[idx, -1]
    experimental = for4_center_temp[int(exp_idx)]
    diff = numerical - experimental
    return diff**2

def F4_num_method_first(h):
    for m in range(1, len(num_for4_time)):
        for i in range(1, int(diameter / (2 * dr)) - 1):
            a_i = s - s/(2*i)
            b_i = 1 - 2*s
            c_i = s + s/(2*i)
            for4_T[m, i] = a_i * for4_T[m-1, i-1] + b_i * for4_T[m-1, i] + c_i * for4_T[m-1, i+1]
            for4_T[m, 0] = for4_T[m, 1]
            for4_T[m, -1] = ((beta * for4_T[m-1, -1]) + (h * sur_area * for4_T_inf)) / (beta + (h *
                                                                                               sur_area))           
    print("Completed Solution")
    residuals = []
    for i in range(len(for4_time)):
        res = F4_sq_diff(for4_time[i])
        residuals.append(res)
    residuals = np.array(residuals[:-1], dtype = float)
    return float(sum(residuals))

# guess = float(input("Guess the h value for forced convection "))
guess = 30

for4_ans = so.least_squares(F4_num_method_first, guess, bounds = [guess-2, guess+2])

print("Numerical Solution for Forced Convection is " + str(round(float(for4_ans.x), 2)) + " ± " +
      str(0.0119))


def F4_num_method_final(h):
    for m in range(1, len(num_for4_time)):
        for i in range(1, int(diameter/(2 * dr))-1):
            a_i = s - s/(2*i)
            b_i = 1 - 2*s
            c_i = s + s/(2*i)
            for4_T[m, i] = a_i * for4_T[m-1, i-1] + b_i * for4_T[m-1, i] + c_i * for4_T[m-1, i+1]
            for4_T[m, 0] = for4_T[m, 1]
            for4_T[m, -1] = ((beta * for4_T[m-1, -1]) + (h * sur_area * for4_T_inf)) / (beta + (h *
                                                                                               sur_area))

    print("Completed Solution")
    return(for4_T)

for4_T = F4_num_method_final(for4_ans.x)

plt.plot(num_for4_time, for4_T[:, -1], 'm', ls='-', lw=2, alpha=0.5, label=f'Numerical Forced Convection (velocity = {velocity[2]}) Solution')






step = 20
new_for2_time = []
new_for2_center_temp = []
for i in range(int(len(for2_time)/step)):
    ft = for2_time[i * step]
    new_for2_time.append(ft)
    ftemp = for2_center_temp[i * step]
    new_for2_center_temp.append(ftemp)

step = 20
new_for3_time = []
new_for3_center_temp = []
for i in range(int(len(for3_time)/step)):
    ft = for3_time[i * step]
    new_for3_time.append(ft)
    ftemp = for3_center_temp[i * step]
    new_for3_center_temp.append(ftemp)

step = 20
new_for4_time = []
new_for4_center_temp = []
for i in range(int(len(for4_time)/step)):
    ft = for4_time[i * step]
    new_for4_time.append(ft)
    ftemp = for4_center_temp[i * step]
    new_for4_center_temp.append(ftemp)

plt.scatter(new_for2_time, new_for2_center_temp, c='b', s=20, alpha=1, label='Forced Convection Data 2')
plt.plot(for2_time, (for2_T_init - for2_T_inf) * np.exp(-(for2_time / for2_tau)) + for2_T_inf, c='b',
         ls=':', label='Forced Analytical Solution 2')

plt.scatter(new_for3_time, new_for3_center_temp, c='g', s=20, alpha=1, label='Forced Convection Data 3')
plt.plot(for3_time, (for3_T_init - for3_T_inf) * np.exp(-(for3_time / for3_tau)) + for3_T_inf, c='g',
         ls=':', label='Forced Analytical Solution 3')

plt.scatter(new_for4_time, new_for4_center_temp, c='m', s=20, alpha=1, label='Forced Convection Data 4')
plt.plot(for4_time, (for4_T_init - for4_T_inf) * np.exp(-(for4_time / for4_tau)) + for4_T_inf, c='m',
         ls=':', label='Forced Analytical Solution 4')

plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Temperature (˚C)')

plt.show()
