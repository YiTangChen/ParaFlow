import matplotlib.pyplot as plt
import numpy as np

latResolution = 1200
lonResolution = 2400
nblocks = 256

with open(f"drawSubdomain/{latResolution}_{lonResolution}.bin", "rb") as f:
    data = np.fromfile(f, dtype=np.float64)
    data = data.reshape((latResolution, lonResolution, 1))

data = np.clip(data, 0.0, 1.0)  # Ensure data is within [0, 1] for visualization
plt.figure(figsize=(10, 5))
plt.imshow(np.clip(data, 0.0, 1.0), cmap='ocean', origin='lower', vmin=0.0, vmax=1.0)
plt.title(f"Subdomain Visualization ({latResolution}x{lonResolution})")
plt.savefig(f"png/subdomain_{latResolution}_{lonResolution}.png")
plt.close()

for i in range(nblocks):
    with open(f"drawSubdomain/{latResolution}_{lonResolution}_{i}.bin", "rb") as f:
        data = np.fromfile(f, dtype=np.float64)
        data = data.reshape((latResolution, lonResolution, 1))
        data = np.clip(data, 0.0, 1.0)  # Ensure data is within [0, 1] for visualization
        plt.figure(figsize=(10, 5))
        plt.imshow(np.clip(data, 0.0, 1.0), cmap='ocean', origin='lower', vmin=0.0, vmax=1.0)
        plt.title(f"Subdomain Visualization {i} ({latResolution}x{lonResolution})")
        plt.savefig(f"png/subdomain_{latResolution}_{lonResolution}_{i}.png")
        plt.close()