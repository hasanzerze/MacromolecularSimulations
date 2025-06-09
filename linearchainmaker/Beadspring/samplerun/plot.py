import matplotlib.pyplot as plt
import numpy as np

# Add labels with Times New Roman font and specific font size
font = {'family': 'serif', 'size': 18}  # Adjust 'size' as needed

# Read data from file
filename = 'thermo_npt.dat'  # Replace with your actual file name
variables = np.loadtxt(filename, unpack=True, comments="#")

time = variables[0] / 1000.
temp = variables[1]
press = variables[2]
totenergy = variables[5]

# Create the plot
#plt.plot(time, temp, marker='o', ms=2.5, color = 'black', linestyle='none')
#plt.plot(time, temp, 'k-')

#plt.plot(time, press, marker='o', ms=2.5, color = 'black', linestyle='none')
#plt.plot(time, press, 'k-')

plt.plot(time, totenergy, marker='o', ms=2.5, color = 'black', linestyle='none')
plt.plot(time, totenergy, 'k-')

# Add labels and title
#plt.xlabel('Time ('+ r'$\sigma m^{1/2} \epsilon ^{-1/2}$'+')', fontdict = font)
#plt.ylabel('Temperature ('+ r'$\epsilon / k_{B}$'+')', fontdict = font)
#plt.title('y vs. x')

#plt.xlabel('Time ('+ r'$\sigma m^{1/2} \epsilon ^{-1/2}$'+')', fontdict = font)
#plt.ylabel('Pressure ('+ r'$\epsilon \sigma^{-3}$'+')', fontdict = font)

plt.xlabel('Time ('+ r'$\sigma m^{1/2} \epsilon ^{-1/2}$'+')', fontdict = font)
plt.ylabel('Total Energy ('+ r'$\epsilon$'+')', fontdict = font)

# Set font properties for tick labels
plt.xticks(fontfamily='serif', fontsize=16)  # Adjust fontsize as needed
plt.yticks(fontfamily='serif', fontsize=16)  # Adjust fontsize as needed

plt.xlim(xmin = 0)
#plt.ylim(ymin = 0.92, ymax=1.33)

# Auto-adjust layout to avoid overlap
plt.tight_layout()

# Save plot to a PDF file
#plt.savefig('Temp_vs_time_NPT.pdf')
#plt.savefig('Press_vs_time_NPT.pdf')
plt.savefig('Energy_vs_time_NPT.pdf')

# Optionally show the plot
#plt.show()

