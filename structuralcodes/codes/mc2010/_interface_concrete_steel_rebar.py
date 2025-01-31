"""
Provides functions to determine parameters of the bond stress-slip relations as shown in Figure 6.1-1 of the fib
ModelCode 2010.

The BondStressSlip class includes methods to calculate:
- Maximum bond stress in case of pull-out failure (tau_bmax).
- Maximum bond stress in case of splitting failure (tau_bu_split).
- Slip corresponding to tau_bu_split (s_tau_bu_split).
- Maximum bond stress due to reinforcement yield (tau_yield).
"""

import warnings


class BondStressSlip:
    def __init__(self, f_cm: float, bond: str):
        """
        Initialize the BondStressSlip class with concrete strength and bond condition.

        Args:
            f_cm: Mean concrete cylindrical compressive strength in MPa.
            bond: "Good" or "Other" bond conditions, according to MC2010 section 6.1.3.2.
        """

        if bond not in ["Good", "Other"]:
            raise ValueError("Invalid input for bond. Use 'Good' or 'Other'")

        self.f_cm = f_cm
        self.bond = bond
        self._tau_bmax = None  # Initialize as None
        self._tau_bu_split = None  # Initialize as None
        self._s_tau_bu_split = None  # Initialize as None

    def tau_bmax(self) -> float:
        """
        Calculate maximum bond stress in case of pull-out failure, according to MC2010 Table 6.1-1

        Returns:
            tau_bmax in MPa.
        """
        if self._tau_bmax is None:  # Calculate only if not already calculated
            if self.bond == 'Good':
                self._tau_bmax = 2.5 * self.f_cm**0.5
            if self.bond == 'Other':
                self._tau_bmax = 1.25 * self.f_cm**0.5
        return self._tau_bmax

    def tau_bu_split(
        self,
        phi: float,
        c_min: float,
        c_max: float,
        k_m: float,
        K_tr: float,
    ) -> float:
        """
        Calculate maximum bond stress in case of splitting failure, Eq.(6.1-5) of MC2010.
        Raises warnings for invalid parameter input values according to MC2010 Eq.(6.1-19), from which Eq.(6.1-5) is
        derived.
        Ensures that tau_bu_split cannot exceed tau_bmax.
        Note: Power of (c_min / phi) deviates from MC2010 Eq. 6.1-5, instead it follows MC2010 Eq. 6.1-19 and fib
        Bulletin 72, Appendix A.

        Args:
            phi: Nominal bar diameter in mm;
            c_min: Parameter according to MC2010 Figure 6.1-2;
            c_max: Parameter according to MC2010 Figure 6.1-2;
            k_m: Parameter according to MC2010 Figure 6.1-3;
            K_tr: To be calculated with MC2010 Eq.(6.1-6);

        Returns:
            tau_bu_split in MPa.
        """
        if (
            self._tau_bu_split is None
        ):  # Calculate only if not already calculated
            # Raise warnings
            if not 15 < self.f_cm < 110:
                warnings.warn(
                    "Warning: Eq.(6.1-19) is valid for 15 MPa < f_cm < 110 MPa",
                    UserWarning,
                )
            if not 0.5 < c_min / phi < 3.5:
                warnings.warn(
                    "Warning: Eq.(6.1-19) is valid for 0.5 < c_min / phi < 3.5",
                    UserWarning,
                )
            if not 1.0 < c_max / c_min < 5.0:
                warnings.warn(
                    "Warning: Eq.(6.1-19) is valid for 1.0 < c_max / c_min < 5.0",
                    UserWarning,
                )
            if not K_tr <= 0.05:
                warnings.warn(
                    "Warning: Eq.(6.1-19) is valid for K_tr <= 0.05",
                    UserWarning,
                )

            if self.bond == 'Good':
                eta_2 = 1.0
            if self.bond == 'Other':
                eta_2 = 0.7

            self._tau_bu_split = (
                eta_2
                * 6.5
                * (self.f_cm / 25) ** 0.25
                * (25 / phi) ** 0.2
                * (
                    (c_min / phi)
                    ** 0.25  # <-- Deviates from MC2010 Eq. 6.1-5, but follows MC2010 Eq. 6.1-19 and fib Bulletin 72, Appendix A.
                    * (c_max / c_min) ** 0.1
                    + k_m * K_tr
                )
            )

            return self._tau_bu_split

    def s_tau_bu_split(
        self,
        tau_bu_split: float,
    ) -> float:
        """
        Calculate the slip corresponding to tau_bu_split, MC2010 Eq.(6.1-1).

        Args:
            tau_bu_split: Concrete bond stress at splitting failure.

        Returns:
            The slip corresponding to tau_bu_split.
        """
        if (
            self._s_tau_bu_split is None
        ):  # Calculate only if not already calculated
            if self.bond == 'Good':
                s1_PO = 1.0
            elif self.bond == 'Other':
                s1_PO = 1.8

            tau_bmax = self.tau_bmax()
            self._s_tau_bu_split = s1_PO * (tau_bu_split / tau_bmax) ** (
                1 / 0.4
            )

        return self._s_tau_bu_split

    @staticmethod
    def tau_yield(
        f_ym: float,
        l_b: float = None,
        phi: float = None,
    ) -> float:
        """
        Calculate the maximum bond stress due to reinforcement yield, based on Eq.(6.1-19).
        Both l_b and phi must be provided together or both must have no input. In the latter case, a uniform bond stress
        over 5 times the bar diameter is assumed, following MC2010 Eq.(6.1-5)

        Args:
            f_ym: Reinforcement mean yield stress.
            l_b: Assumed length with uniform bond stress.
            phi: Reinforcement diameter.

        Returns:
            Bond stress at reinforcement yield.
        """
        if (l_b is None and phi is not None) or (
            l_b is not None and phi is None
        ):
            raise ValueError(
                "Both l_b and phi must be provided together or both must be excluded."
            )

        # If both are excluded, assume uniform bond stress over l_b = 5 * phi
        if l_b is None and phi is None:
            l_b = 5
            phi = 1

        # tau_yield derived from the equilibrium: f_ym * 0.25 * pi * phi ** 2 = pi * phi * l_b * tau_yield
        return f_ym * phi / 4 / l_b
