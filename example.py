from time import time
from typing import List
import numpy as np
import matplotlib.pyplot as plt
from gillespy import gillespie


def lotka(y: List[int]) -> List[float]:
    return [y[0], y[0] * y[1], y[1]]


updater = [[1, 0], [-1, 1], [0, -1]]

start = time()
t, species = gillespie([1000, 1000], [10.0, 0.01, 10.0], updater, lotka, 30)
end = time()
print(f"Simulation took {end - start:.2f} seconds")

t = np.array(t)
species = np.array(species)

fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(9, 3), gridspec_kw=dict(wspace=0.4))

ax0.plot(t, species[:, 0], lw=0.5, label="Species 1")
ax0.plot(t, species[:, 1], lw=0.5, label="Species 2")
ax0.set_xlabel("Time")
ax0.set_ylabel("Abundance")
ax0.legend()

ax1.plot(species[:, 0], species[:, 1], lw=0.5)
ax1.set_xlabel("Species 1")
ax1.set_ylabel("Species 2")

plt.savefig("figs/lokia.png", dpi=300, bbox_inches="tight")
