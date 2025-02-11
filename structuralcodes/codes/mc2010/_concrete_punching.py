"""Covers punching in Model code 2010, 7.3.5.1 to 7.3.5.4."""

import warnings
from math import cos, pi, sin
from dataclasses import dataclass


@dataclass
class ColumnPosition:
    """Column position configuration."""
    inner: bool
    edge_par: bool  # edge column with tension reinf. parallel to edge
    edge_per: bool  # edge column with tension reinf. perpendicular to edge
    corner: bool

    def validate(self) -> None:
        """Validate that only one position is True."""
        positions = [self.inner, self.edge_par, self.edge_per, self.corner]
        if sum(positions) != 1:
            raise ValueError('Exactly one column position must be True')


class ConcretePunching:
    """
    The goal of the class is to calculate the punching shear resistance
    calculations according to MC2010.
    """

    def __init__(
        self,
        l_x: float,
        l_y: float,
        f_yd: float,
        d: float,
        e_s: float,
        v_ed: float,
        e_u: float,
        inner: bool,
        edge_par: bool,
        edge_per: bool,
        corner: bool,
        m_rd: float,
        m_pd: float,
        x_direction: bool,
        gamma_c: float = 1.5,
        gamma_s: float = 1.15,
    ):
        """Initialize ConcretePunching calculator.

        Args:
            l_x (float): The distance between two columns in x direction.
            l_y (float): The distance between two columns in y direction.
            f_yd (float): Design strength of reinforment steel in MPa.
            d (float): The mean value of the effective depth in mm.
            e_s (float): The E_modulus for steel in MPa.
            v_ed (float): The acting shear force from the columns.
            e_u (float): Refers to the eccentricity of the resultant of shear
                forces with respect to the centroid.
            inner (bool): Is true only if the column is a inner column.
            edge_par (bool): Is true only if the column is a edge column with
                tension reinforcement parallel to the edge.
            edge_per (bool): Is true only if the column is a edge column with
                tension reinforcement perpendicular to the edge.
            corner (bool): Is true only if the column is a corner column.
            m_rd (float): The design average strength per unit length in MPa.
            m_pd (float): The average decompresstion moment due to
                prestressing.
            x_direction (bool): True if x direction is considered.
            gamma_c (float): Safety factor for concrete. Defaults to 1.5.
            gamma_s (float): Safety factor for reinforcement. Defaults to 1.15.
        """
        # Geometric parameters
        self.l_x = l_x
        self.l_y = l_y
        self.d = d

        # Material properties
        self.f_yd = f_yd
        self.e_s = e_s

        # Loading parameters
        self.v_ed = v_ed
        self.e_u = e_u

        # Column position configuration
        self.position = ColumnPosition(inner, edge_par, edge_per, corner)
        self.position.validate()

        # Design parameters
        self.m_rd = m_rd
        self.m_pd = m_pd
        self.x_direction = x_direction

        # Safety factors
        self.gamma_c = gamma_c
        self.gamma_s = gamma_s

    def b_0(self, v_perp_d_max: float) -> float:
        """Calculate shear-resisting control perimeter b_0.

        fib Model Code 2010, eq. (7.3-57).

        Args:
            v_perp_d_max (float): The maximum shear force per unit length
                perpendicular to the basic control parameter (Figure 7.3-24).

        Returns:
            float: The shear-resisting control perimeter, b_0.
        """
        return self.v_ed / v_perp_d_max

    def m_ed(self) -> float:
        """Calculate average bending moment in support strip.

        fib Model Code 2010, eq. (7.3-76), (7.3-71), (7.3-72), (7.3-73)
        and (7.3-74).

        Returns:
            float: The bending moment acting in the support strip.
        """
        r_sx = 0.22 * self.l_x
        r_sy = 0.22 * self.l_y
        l_min = min(self.l_x, self.l_y)
        b_s = min(1.5 * (r_sx * r_sy) ** 0.5, l_min)

        if self.position.inner:
            return self.v_ed * ((1 / 8) + abs(self.e_u) / (2 * b_s))
        if self.position.edge_par:
            return max(
                self.v_ed * ((1 / 8) + self.e_u / (2 * b_s)),
                self.v_ed / 4
            )
        if self.position.edge_per:
            return self.v_ed * ((1 / 8) + self.e_u / (b_s))
        if self.position.corner:
            return max(
                self.v_ed * ((1 / 8) + self.e_u / (b_s)),
                self.v_ed / 2
            )
        raise ValueError('Placement is not defined, only one needs to be True')

    def psi_punching(self, approx_lvl_p: float) -> float:
        """Calculate rotation of slab around supported area.

        fib Model Code 2010, eq. (7.3-70), (7.3-75) and (7.3-77).

        Args:
            approx_lvl_p (float): The approx level for punching.

        Returns:
            float: psi for the chosen approx level in punching.
        """
        r_s = max(0.22 * self.l_x, 0.22 * self.l_y)
        if approx_lvl_p == 1:
            psi = 1.5 * r_s * self.f_yd / (self.d * self.e_s)
        elif approx_lvl_p == 2:
            r_s = 0.22 * self.l_x if self.x_direction else 0.22 * self.l_y
            psi = (1.5 * r_s * self.f_yd / (self.d * self.e_s)) * (
                    self.m_ed() / self.m_rd
                ) ** 1.5
        return psi

    def v_rdc_punching(
        self,
        approx_lvl_p: float,
        dg: float,
        f_ck: float,
        d_v: float,
        v_perp_d_max: float,
    ) -> float:
        """Calculate punching resistance from concrete.

        fib Model Code 2010, eq. (7.3-61), (7.3-62) and (7.3-63).

        Args:
            approx_lvl_p (float): The approx level for punching.
            dg (float): Maximum size of aggregate.
            f_ck (float): Characteristic strength in MPa.
            d_v (float): The effective depth considering support in mm.
            v_perp_d_max (float): The maximum shear force per unit length
                perpendicular to the basic control parameter (Figure 7.3-24).

        Returns:
            float: v_rdc for punching with the right approx level.
        """
        k_dg = max(32 / (16 + dg), 0.75)
        k_psi = min(
            1 / (1.5 + 0.9 * k_dg * self.d * self.psi_punching(approx_lvl_p)),
            0.6,
        )
        return (
            k_psi * self.b_0(v_perp_d_max) * d_v * (f_ck**0.5) / self.gamma_c
        )

    def v_rds_punching(
        self,
        b_u: float,
        approx_lvl_p: float,
        alpha: float,
        f_bd: float,
        f_ywk: float,
        phi_w: float,
        a_sw: float,
    ) -> float:
        """Calculate punching resistance from shear reinforcement.

        fib Model Code 2010, eq. (7.3-64) and (7.3-65).

        Args:
            b_u (float): The diamter of a circle with same surface as the
                region inside the basic control perimeter (Figure 7.3-27b).
            approx_lvl_p (float): The approx level for punching.
            alpha (float): Inclination of the stirrups in degrees.
            f_bd (float): The design bond strength in MPa.
            f_ywk (float): Characteristic yield strength of shear
                reinforcement.
            phi_w (float): The diameter of the shear reinforcement.
            a_sw (float): The area of the shear reinforcement in mm^2.

        Returns:
            float: Punching resistance that comes from reinforcement.
        """
        f_ywd = f_ywk / self.gamma_s
        k_e = 1 / (1 + self.e_u / b_u)
        sigma_swd = min(
            (self.e_s * self.psi_punching(approx_lvl_p) / 6)
            * (sin(alpha * pi / 180) + cos(alpha * pi / 180))
            * (sin(alpha * pi / 180) + f_bd * self.d / (f_ywd * phi_w)),
            f_ywd,
        )

        if (a_sw * k_e * f_ywd) < 0.5 * self.v_ed:
            warnings.warn(
                "Consider increasing punching shear reinforcement for "
                "sufficient deformation capacity"
            )
        return a_sw * k_e * sigma_swd * sin(alpha * pi / 180)

    def v_rd_max_punching(
        self,
        approx_lvl_p: float,
        dg: float,
        v_perp_d_max: float,
        d_v: float,
        f_ck: float,
        d_head: bool,
        stirrups_compression: bool,
    ) -> float:
        """Calculate maximum value for v_rd_punching.

        fib Model Code 2010, eq. (7.3-68) and (7.3-69).

        Args:
            approx_lvl_p (float): The approx level for punching.
            dg (float): Maximum size of aggregate.
            v_perp_d_max (float): The maximum shear force per unit length
                perpendicular to the basic control parameter (Figure 7.3-24).
            d_v (float): The effective depth considering support in mm.
            f_ck (float): Characteristic strength in MPa.
            d_head (bool): True if diameter of heads is three times larger.
            stirrups_compression: (bool): Stirrups with sufficient length at
                compression face, and bent on tension face.

        Return:
            float: The maximum allowed punching resistance.
        """
        if d_head:
            k_sys = 2.8
        elif stirrups_compression:
            k_sys = 2.4
        else:
            k_sys = 2

        k_dg = max(32 / (16 + dg), 0.75)
        k_psi = min(
            1 / (1.5 + 0.9 * k_dg * self.d * self.psi_punching(approx_lvl_p)),
            0.6,
        )
        return min(
            (k_sys * k_psi * self.b_0(v_perp_d_max) * d_v * f_ck**0.5)
            / self.gamma_c,
            (self.b_0(v_perp_d_max) * d_v * f_ck**0.5) / self.gamma_c,
        )

    def v_rd_punching(
        self,
        b_u: float,
        approx_lvl_p: float,
        alpha: float,
        f_bd: float,
        f_ywk: float,
        phi_w: float,
        a_sw: float,
        dg: float,
        f_ck: float,
        d_v: float,
        v_perp_d_max: float,
        d_head: bool,
        stirrups_compression: bool,
    ) -> float:
        """Calculate total punching resistance (Vrd,c + Vrd,s).

        fib Model Code 2010, eq. (7.3-60).

        Args:
            b_u (float): The diamter of a circle with same surface as the
                region inside the basic control perimeter (Figure 7.3-27b).
            approx_lvl_p (float): The approx level for punching.
            alpha (float): Inclination of the stirrups in degrees.
            f_bd (float): The design bond strength in MPa.
            f_ywk (float): Characteristic yield strength of shear
                reinforcement.
            phi_w (float): The diameter of the shear reinforcement.
            a_sw (float): The area of the shear reinforcement in mm^2.
            dg (float): Maximum size of aggregate.
            f_ck (float): Characteristic strength in MPa.
            d_v (float): The effective depth considering support in mm.
            v_perp_d_max (float): The maximum shear force per unit length
                perpendicular to the basic control parameter (Figure 7.3-24).
            d_head (bool): True if diameter of heads is three times larger.
            stirrups_compression (bool): Stirrups with sufficient length at
                compression face, and bent on tension face.

        Return:
            float: The maximum allowed punching resistance, regardless of
            values from v_rdc and v_rds.
        """
        return min(
            self.v_rdc_punching(
                approx_lvl_p,
                dg,
                f_ck,
                d_v,
                v_perp_d_max,
            )
            + self.v_rds_punching(
                b_u,
                approx_lvl_p,
                alpha,
                f_bd,
                f_ywk,
                phi_w,
                a_sw,
            ),
            self.v_rd_max_punching(
                approx_lvl_p,
                dg,
                v_perp_d_max,
                d_v,
                f_ck,
                d_head,
                stirrups_compression,
            ),
        )
