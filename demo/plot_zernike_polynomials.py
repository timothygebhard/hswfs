"""
Create plot of first 15 Zernike polynomials.
"""

# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np

from hswfs.utils import get_unit_disk_meshgrid
from hswfs.zernike import eval_cartesian, ZernikePolynomial


# -----------------------------------------------------------------------------
# MAIN CODE
# -----------------------------------------------------------------------------

if __name__ == "__main__":

    # Define the maximum value of n
    n_max = 4

    # Set up a new figure
    fig = plt.figure(figsize=(n_max, n_max))

    # Loop over the indices m and n of the Zernike polynomials
    for row, n in enumerate(range(0, n_max + 1)):
        for col, m in enumerate(range(-n, n + 1, 2)):

            print(f"Running for Z^{m}_{n}...", end=" ", flush=True)

            # Create new Zernike polynomial
            zernike = ZernikePolynomial(m=m, n=n)

            # Add a new subplot to the plot
            ax = plt.subplot2grid(
                (n_max + 1, 2 * (n_max + 1)),
                (n, n_max - n + 2 * col),
                colspan=2,
            )

            # Evaluate the Zernike polynomial on the unit disk
            x_0, y_0 = get_unit_disk_meshgrid(resolution=512)
            img = eval_cartesian(zernike.cartesian, x_0, y_0)
            limit = np.nanmax(np.abs(img))

            # Actually plot the Zernike polynomial
            ax.imshow(
                img, origin="lower", cmap="RdBu_r", vmin=-limit, vmax=limit
            )

            # Add plot options (title, disable ticks etc.)
            ax.set_title(rf"$Z^{{{m}}}_{{{n}}}$", fontsize=4, pad=4)
            ax.axis("off")

            print("Done!", flush=True)

    # Final plot options; save as PNG
    plt.tight_layout(h_pad=0.1, w_pad=0.1)
    plt.savefig("zernike_polynomials.png", dpi=300, pad=0)
