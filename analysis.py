import numpy as np 
import matplotlib.pyplot as plt 


labels = ['R1-D52', 'R2-D22', 'R2-D36', 'R2-D52', 'R3-D51', 'R4-D22', 'R4-D36', 'R4-D52', 'R5-D22', 'R5-D36', 'R5-D52', 'R6-D24', 'R6-D49', 'R7-D16', 'R7-D61']
avg_vel = [13.298, 14.986, 14.754, 15.450, 16.294, 15.269, 13.224, 12.584, 14.32, 12.62, 11.91, 8.203, 9.562, 20.164, 15.540]
r = [52, 22, 36, 52, 51, 22, 36, 52, 22, 36, 52, 24, 49, 16, 61]
colors = ['brown', 'black', 'black', 'black', 'green', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'orange', 'orange', 'pink', 'pink']

fig, ax = plt.subplots()
ax.scatter(r, avg_vel, s=10, edgecolor='black')
ax.plot(r[1:4], avg_vel[1:4], color='black')
ax.plot(r[5:8], avg_vel[5:8], color='red')
ax.plot(r[8:11], avg_vel[8:11], color='blue')
ax.plot(r[11:13], avg_vel[11:13], color='orange')
ax.plot(r[13:], avg_vel[13:], color='pink')

for i, txt in enumerate(labels):
    ax.annotate(txt, (r[i], avg_vel[i]), color=colors[i])

plt.title("Average velocity differences (all species) vs Distance from star")
plt.xlabel("Radius (Ny)")
plt.ylabel("Average velocity difference (m/s)")
plt.xlim(0, 80)
plt.show()


alpha = [10**(-4), 10**(-3), 10**(-3), 10**(-3), 0.01, 10**(-5), 10**(-5), 10**(-5), 10**(-6), 10**(-6), 10**(-6)]

fig, ax = plt.subplots()
ax.scatter(alpha, avg_vel[0:11], s=10, edgecolor='black')
ax.set_xscale('log')

for i, txt in enumerate(labels[0:11]):
    ax.annotate(txt, (alpha[i], avg_vel[i]), color=colors[i])

plt.title("Average velocity differences (all species) vs ALPHA")
plt.xlabel("Alpha Viscosity")
plt.ylabel("Average velocity difference (m/s)")
plt.xlim(10**-7, 0.1)
plt.show()


mass_labels = ['R6-D24', 'R6-D49', 'R1-D52', 'R7-D16', 'R7-D61']
mass = [0.0001482, 0.0001482, 0.0003, 0.001, 0.001]
vel_planet = [8.203, 9.562, 13.298, 20.164, 15.540]
p_col=['red', 'red', 'black', 'blue', 'blue']

fig, ax = plt.subplots()

ax.scatter(mass, vel_planet, s=10, edgecolor='black')
ax.set_xscale('log')

for i, txt in enumerate(mass_labels):
    ax.annotate(txt, (mass[i], vel_planet[i]), color=p_col[i])

# ax.tick_params(axis='x', rotation=70)
plt.axvline(0.0001482, color='red', alpha=0.5, label="Half Saturn")
plt.axvline(0.0003, color='black', alpha=0.5, label="Saturn")
plt.axvline(0.001, color='blue', alpha=0.5, label="Jupiter")

plt.title("Average velocity differences (all species) vs Planet Mass")
plt.xlabel("Planet Mass / Star Mass")
plt.ylabel("Average velocity difference (m/s)")
plt.xlim(10**-5, 10**-2)
plt.legend(loc='upper left')
plt.show()