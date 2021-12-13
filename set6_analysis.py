import numpy as np 
import matplotlib.pyplot as plt 

R = 5.2  # 5.2 AU

vel_avg = np.array([16.2457, 12.4237, 15.1385, 8.50613])
radii = np.array([1, 5.2, 10, 20])

fig, ax = plt.subplots()
ax.plot(radii, vel_avg, marker='o')
# ax.set_xscale('log')

plt.title("Average velocity differences (for all species) at different radii")
plt.xlabel("Radii (AU)")
plt.ylabel("Average velocity difference (m/s)")
plt.show()