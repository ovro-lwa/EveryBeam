from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np


def read_oskar_beams():
    """
    Read oskar beam and store into 

    Returns
    -------
    np.array
        (N, N, 2, 2) array, where N number of image pixels
    """
    A = None
    for i, pol1 in enumerate(["X", "Y"]):
        for j, pol2 in enumerate(["X", "Y"]):
            filename_amp = "basefunctions_050_MHz_S0000_TIME_SEP_CHAN_SEP_AMP_{}{}.fits".format(
                pol1, pol2
            )
            filename_phase = "basefunctions_050_MHz_S0000_TIME_SEP_CHAN_SEP_PHASE_{}{}.fits".format(
                pol1, pol2
            )
            d = (
                fits.getdata(filename_amp) * np.exp(1j * fits.getdata(filename_phase))
            )[0, 0, :, :].T
            if A is None:
                N = d.shape[-1]
                A = np.zeros((N, N, 2, 2), dtype=np.complex128)
            A[:, :, i, j] = d[::-1, :]
        N = A.shape[0]

    for i in range(N):
        x = 2.0 * i / (N - 1) - 1.0
        for j in range(N):
            y = 2.0 * j / (N - 1) - 1.0
            theta = np.arcsin(np.sqrt(x * x + y * y))
            phi = np.arctan2(y, x)

            # Need to swap the sign of phi to get the transformation right
            # Still need to figure out why
            # phi = -phi

            e_theta = np.array([[np.cos(phi)], [np.sin(phi)]])
            e_phi = np.array([[-np.sin(phi)], [np.cos(phi)]])

            # This matrix applies the transformation defined in
            # https://github.com/OxfordSKA/OSKAR/blob/master/oskar/convert/define_convert_ludwig3_to_theta_phi_components.h
            T = np.concatenate((e_theta, e_phi), axis=1)

            # Transformation is applied to the right side of A
            # so the rows a of A are tranformed
            A[i, j, :, :] = np.dot(A[i, j, :, :], T)
    return A


if __name__ == "__main__":

    A = read_oskar_beams()

    plt.subplot(2, 2, 1)
    plt.imshow(np.abs(A[:, :, 0, 0]).T, origin="lower")
    plt.colorbar()
    plt.title("abs(Etheta)")

    plt.subplot(2, 2, 2)
    plt.imshow(np.abs(A[:, :, 0, 1]).T, origin="lower")
    plt.colorbar()
    plt.title("abs(Ephi)")

    plt.subplot(2, 2, 3)
    plt.imshow(
        np.angle(A[:, :, 0, 0]).T, clim=(-np.pi, np.pi), cmap="twilight", origin="lower"
    )
    plt.colorbar()
    plt.title("angle(Etheta)")

    plt.subplot(2, 2, 4)
    plt.imshow(
        np.angle(A[:, :, 0, 1]).T, clim=(-np.pi, np.pi), cmap="twilight", origin="lower"
    )
    plt.colorbar()
    plt.title("angle(Ephi)")

    plt.show()
