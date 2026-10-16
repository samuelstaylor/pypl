import numpy as np
from scipy import constants

# Energy of one unit of each supported SOC unit, in J.
SOC_UNITS = {
    "GHz": constants.h * 1e9,
    "meV": constants.eV * 1e-3,
    "ueV": constants.eV * 1e-6,
    "cm-1": constants.h * constants.c * 1e2,
}


def compute_isc_rate(soc, vof, soc_unit="GHz"):
    r"""
    Intersystem crossing rate from Fermi's golden rule.

    .. math::

        $\Gamma(\Delta) = \frac{4\pi}{\hbar} \, |\lambda_\perp|^2 \, F(\Delta)$

    where :math:`\lambda_\perp` is the spin-orbit coupling matrix element between the two
    electronic states and :math:`F(\Delta)` the vibrational overlap function (VOF), i.e. the
    Franck-Condon lineshape of the transition evaluated at the electronic energy gap
    :math:`\Delta` (e.g. ``fc_lineshape`` of ``config_coord_1d_solver.bulid_fc_lsp`` or
    ``config_coord_1d_vof_solver.build_fc_lsp_from_levels``). The prefactor :math:`4\pi/\hbar` follows
    the convention used for the :math:`{}^3E \rightarrow {}^1A_1` transition of the NV center
    (Goldman et al., Phys. Rev. Lett. 114, 145502 (2015); Jin et al., Phys. Rev. Lett. 135, 036401
    (2025)); scale the result if your definition of :math:`\lambda_\perp` uses a different one.

    Parameters
    ----------
    soc : float
        Spin-orbit coupling matrix element :math:`\lambda_\perp`, in ``soc_unit``.
    vof : float or ndarray
        Vibrational overlap function :math:`F(\Delta)` in 1/meV.
    soc_unit : {'GHz', 'meV', 'ueV', 'cm-1'}, optional
        Unit of ``soc``; GHz means :math:`h \nu` with :math:`\nu` in GHz. Default is 'GHz'.

    Returns
    -------
    float or ndarray
        :math:`\Gamma` in 1/s (Hz), with the shape of ``vof``.
    """

    if soc_unit not in SOC_UNITS:
        raise ValueError("soc_unit must be one of %s, got %r" % (", ".join(SOC_UNITS), soc_unit))
    soc_j = soc * SOC_UNITS[soc_unit]
    vof_per_j = np.asarray(vof, dtype=float) / (constants.eV * 1e-3)
    return 4 * np.pi / constants.hbar * soc_j**2 * vof_per_j


def compute_isc_lifetime(rate):
    r"""
    Lifetime :math:`\tau = 1 / \Gamma` of an intersystem crossing channel.

    Parameters
    ----------
    rate : float or ndarray
        Rate :math:`\Gamma` in 1/s, e.g. from ``compute_isc_rate``.

    Returns
    -------
    float or ndarray
        :math:`\tau` in s; ``inf`` where :math:`\Gamma \le 0`.
    """

    rate = np.asarray(rate, dtype=float)
    with np.errstate(divide="ignore"):
        tau = np.where(rate > 0, 1.0 / np.where(rate > 0, rate, 1.0), np.inf)
    return tau if tau.ndim else float(tau)
