import matplotlib.pyplot as plt

# Load the chain data (e.g., from Rosenbluth sampling)
# Format: each line in 'chain_coords.txt' should contain "x y"
positions = []
with open("chain_coords.txt") as f:
    for line in f:
        x, y = map(int, line.strip().split())
        positions.append((x, y))

# Separate x and y coordinates
x_coords = [p[0] for p in positions]
y_coords = [p[1] for p in positions]

# Create the plot
fig, ax = plt.subplots(figsize=(6, 6))

# Plot the polymer chain
ax.plot(x_coords, y_coords, '-o', color='black', markersize=7, linewidth=2)

# Mark the starting monomer with a white square (black border)
ax.plot(x_coords[0], y_coords[0], marker='s', markerfacecolor='white', markeredgecolor='black',
        markersize=10, label='Start')

# Mark the ending monomer with a black square
ax.plot(x_coords[-1], y_coords[-1], marker='s', color='black', markersize=10, label='End')

# Set axis limits with margin
margin = 2
ax.set_xlim(min(x_coords) - margin, max(x_coords) + margin)
ax.set_ylim(min(y_coords) - margin, max(y_coords) + margin)

# Keep axis labels and ticks
ax.set_xlabel("X (lattice units)")
ax.set_ylabel("Y (lattice units)")

# Add lattice grid lines
ax.set_xticks(range(min(x_coords) - margin, max(x_coords) + margin + 1))
ax.set_yticks(range(min(y_coords) - margin, max(y_coords) + margin + 1))
ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')

# Equal aspect ratio
ax.set_aspect('equal')

# Optional: add legend
ax.legend(loc='upper left')

# Show and/or save
plt.tight_layout()
plt.savefig("polymer_chain_eps-1.png", dpi=300)
#plt.show()


