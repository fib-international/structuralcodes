"""The concrete class for EC2 2023 Concrete Material."""

import typing as t
import warnings

from structuralcodes.codes import ec2_2023

from ..constitutive_laws import ConstitutiveLaw, create_constitutive_law
from ._concrete import Concrete


class ConcreteEC2_2023(Concrete):  # noqa: N801
    """Concrete implementation for EC2 2023 Concrete."""

    _kE: t.Optional[float] = None  # noqa: N815
    _strength_dev_class: t.Optional[str] = None
    _fcm: t.Optional[float] = None
    _fctm: t.Optional[float] = None
    _Ecm: t.Optional[float] = None
    _fctk_5: t.Optional[float] = None
    _fctk_95: t.Optional[float] = None
    _eps_c1: t.Optional[float] = None
    _eps_cu1: t.Optional[float] = None
    _k_sargin: t.Optional[float] = None
    _eps_c2: t.Optional[float] = None
    _eps_cu2: t.Optional[float] = None
    _n_parabolic_rectangular: t.Optional[float] = None

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400.0,
        kE: float = 9500,
        strength_dev_class: t.Literal[
            'CS', 'CN', 'CR', 'slow', 'normal', 'rapid'
        ] = 'CN',
        gamma_c: t.Optional[float] = None,
        constitutive_law: t.Optional[
            t.Union[
                t.Literal[
                    'elastic',
                    'parabolarectangle',
                    'bilinearcompression',
                    'sargin',
                    'popovics',
                ],
                ConstitutiveLaw,
            ]
        ] = 'parabolarectangle',
        fcm: t.Optional[float] = None,
        fctm: t.Optional[float] = None,
        fctk_5: t.Optional[float] = None,
        fctk_95: t.Optional[float] = None,
        Ecm: t.Optional[float] = None,
        eps_c1: t.Optional[float] = None,
        eps_cu1: t.Optional[float] = None,
        k_sargin: t.Optional[float] = None,
        eps_c2: t.Optional[float] = None,
        eps_cu2: t.Optional[float] = None,
        n_parabolic_rectangular: t.Optional[float] = None,
        **kwargs,
    ):
        """Initializes a new instance of Concrete for EC2 2023.

        Arguments:
            fck (float): Characteristic strength in MPa if concrete is not
                existing.

        Keyword Arguments:
            name (str): A descriptive name for concrete.
            density (float): Density of material in kg/m3 (default: 2400).
            kE (float): Coefficient relating aggregates.
            strength_dev_class (str, optional): Default is CN. Possible values:
                CS, CN, CR, slow, normal or rapid.
            gamma_c (float, optional): Partial factor of concrete (default is
                1.5).
            consitutive_law (ConstitutiveLaw | str): A valid ConstitutiveLaw
                object for concrete or a string defining a valid constitutive
                law type for concrete. (valid options for string: 'elastic',
                'parabolarectangle', 'bilinearcompression', 'sargin',
                'popovics').
            fcm (float, optional): The mean compressive strength.
            fctm (float, optional): The mean tensile strength.
            fctk_5 (float, optional): The 5% fractile for the tensile strength.
            fctk_95 (float, optional): The 95% fractile for the tensile
                strength.
            Ecm (float, optional): The mean secant Young's modulus.
            eps_c1 (float, optional): The strain at peak stress for the Sargin
                constitutive law.
            eps_cu1 (float, optional): The ultimate strain for the Sargin
                constitutive law.
            k_sargin (float, optional): The coefficient for the Sargin
                constitutive law.
            eps_c2 (float, optional): The strain at peak stress for the
                parabolic rectangular constitutive law.
            eps_cu2 (float, optional): The ultimate strain for the parabolic
                rectangular constitutive law.
            n_parabolic_rectangular (float, optional): The coefficient for the
                parabolic rectangular constitutive law.

        Raises:
            ValueError: If the provided strength development class is invalid.
            ValueError: If gamma_c is less than zero.
            ValueError: If fcm is less than fck.
            ValueError: If k_sargin is negative.
            ValueError: If n_parabolic_rectangular is negative.
            ValueError: If the constitutive law name is not available for the
                material.
            ValueError: If the provided constitutive law is not valid for
                concrete.
            Warning: If eps_c1 is larger than 0.1.
            Warning: If eps_cu1 is larger than 0.1.
            Warning: If eps_c2 is larger than 0.1.
            Warning: If eps_cu2 is larger than 0.1.
        """
        del kwargs
        if name is None:
            name = f'C{round(fck):d}'

        super().__init__(
            fck=fck,
            name=name,
            density=density,
            existing=False,
            gamma_c=gamma_c,
        )
        self._kE = kE
        self._strength_dev_class = strength_dev_class.strip().lower()
        self._fcm = abs(fcm) if fcm is not None else None
        self._fctm = abs(fctm) if fctm is not None else None
        self._fctk_5 = abs(fctk_5) if fctk_5 is not None else None
        self._fctk_95 = abs(fctk_95) if fctk_95 is not None else None
        self._Ecm = abs(Ecm) if Ecm is not None else None
        self._eps_c1 = abs(eps_c1) if eps_c1 is not None else None
        self._eps_cu1 = abs(eps_cu1) if eps_cu1 is not None else None
        self._k_sargin = k_sargin if k_sargin is not None else None
        self._eps_c2 = abs(eps_c2) if eps_c2 is not None else None
        self._eps_cu2 = abs(eps_cu2) if eps_cu2 is not None else None
        self._n_parabolic_rectangular = (
            n_parabolic_rectangular
            if n_parabolic_rectangular is not None
            else None
        )

        self.__post_init__()

        # The constitutive law requires valid attributes, so it should be set
        # after validation
        self._constitutive_law = (
            constitutive_law
            if isinstance(constitutive_law, ConstitutiveLaw)
            else create_constitutive_law(
                constitutive_law_name=constitutive_law, material=self
            )
        )
        if 'concrete' not in self._constitutive_law.__materials__:
            raise ValueError(
                'The provided constitutive law is not valid for concrete.'
            )

    def __post_init__(self):
        """Validator for the attributes that are set in the constructor."""
        # Check if strength_dev_class is valid
        valid_dev_classes = ('cs', 'cn', 'cr', 'slow', 'normal', 'rapid')
        if self._strength_dev_class not in valid_dev_classes:
            raise ValueError(
                'strength_dev_class not valid. '
                + f'{self._strength_dev_class} not {valid_dev_classes}'
            )

        # Check for valid partial factor
        if self._gamma_c is not None and self._gamma_c <= 0:
            raise ValueError(
                f'gamma_c cannot be less than zero. Current: {self._gamma_c}'
            )

        # Check the mean compressive strength
        if self._fcm is not None and self._fcm <= self._fck:
            raise ValueError(
                (
                    'Mean compressive strength cannot be lower than',
                    'characteristic strength.\n',
                    'Current characteristing strength: ',
                    f'fck = {self._fck}.',
                    f'Current value: value = {self._fcm}',
                )
            )

        # Check eps_c1
        if self._eps_c1 is not None and self._eps_c1 >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_c1 that should be a pure'
                f' number without units. Please check ({self._eps_c1} given).'
            )

        # Check eps_cu1
        if self._eps_cu1 is not None and self._eps_cu1 >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_cu1 that should be a pure'
                f' number without units. Please check ({self._eps_cu1} given).'
            )

        # Check k_sargin
        if self._k_sargin is not None and self._k_sargin < 0:
            raise ValueError(
                f'k_sargin should be a positive value ({self._k_sargin} given)'
            )

        # Check eps_c2
        if self._eps_c2 is not None and self._eps_c2 >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_c2 that should be a pure'
                f' number without units. Please check ({self._eps_c2} given).'
            )

        # Check eps_cu2
        if self._eps_cu2 is not None and self._eps_cu2 >= 0.1:
            warnings.warn(
                'A suspect value is input for eps_cu2 that should be a pure'
                f' number without units. Please check ({self._eps_cu2} given).'
            )

        # Check n_parabolic_rectangular
        if (
            self._n_parabolic_rectangular is not None
            and self._n_parabolic_rectangular < 0
        ):
            raise ValueError(
                'n should be a positive value '
                f'({self._n_parabolic_rectangular} given)'
            )
        if (
            self._n_parabolic_rectangular is not None
            and self._n_parabolic_rectangular >= 5
        ):
            warnings.warn(
                'A suspect value is input for eps_cu2 that should be a pure'
                ' number without units. Please check '
                f'({self._n_parabolic_rectangular} given).'
            )

    @property
    def fcm(self) -> float:
        """Returns the mean strength of concrete.

        Returns:
            float: The mean compressive strength in MPa.

        Note:
            The returned value is derived from fck if fcm is not manually
            provided when initializing the object.
        """
        if self._fcm is None:
            return ec2_2023.fcm(self.fck)
        return self._fcm

    @property
    def fctm(self) -> None:
        """Returns the mean concrete tensile strength.

        Returns:
            float: The mean concrete tensile strength.

        Note:
            The returned value is derived from fck if fctm is not manually
            provided when initializing the object.
        """
        if self._fctm is None:
            return ec2_2023.fctm(self.fck)
        return self._fctm

    @property
    def fctk_5(self) -> float:
        """Returns the 5% fractile for the concrete tensile strength.

        Returns:
            float: The 5% fractile for the concrete tensile strength in MPa.

        Note:
            The returned value is derived from fctm if fctk_5 is not manually
            provided when initializing the object.
        """
        if self._fctk_5 is None:
            return ec2_2023.fctk_5(self.fctm)
        return self._fctk_5

    @property
    def fctk_95(self) -> float:
        """Returns the 95% fractile for the concrete tensile strength.

        Returns:
            float: The 5% fractile for the concrete tensile strength in MPa.

        Note:
            The returned value is derived from fctm if fctk_95 is not manually
            provided when initializing the object.
        """
        if self._fctk_95 is None:
            return ec2_2023.fctk_95(self.fctm)
        return self._fctk_95

    @property
    def Ecm(self) -> float:
        """Returns the secant modulus.

        Returns:
            float: The secant concrete modulus in MPa.

        Note:
            The returned value is derived from fcm and kE if Ecm is not
            manually provided when initializing the object.
        """
        if self._Ecm is None:
            return ec2_2023.Ecm(self.fcm, self._kE)
        return self._Ecm

    def fcd(
        self, t_ref: float = 28, t0: float = 91, fck_ref: float = 40
    ) -> float:
        """Computes the value of the design compressive strength of concrete.

        Arguments:
            t_ref (float,optional): The reference time in days (default is 28
                days).
            t0 (float, optional): Age at loading in days (default is 91 days).
            fck_ref (float, optional): The reference compressive strength in
                MPa (default is 40 MPa).

        Returns:
            float: The design compressive strength of concrete in MPa.

        Raises:
            ValueError: If fck_ref is less or equal to 0.
            ValueError: If t_ref is less than 0.
            ValueError: If t0 is less than 0.
        """
        eta_cc = ec2_2023.eta_cc(self.fck, fck_ref=fck_ref)
        k_tc = ec2_2023.k_tc(t_ref, t0, self._strength_dev_class)
        return ec2_2023.fcd(self.fck, eta_cc, k_tc, self.gamma_c)

    def fctd(self, t_ref: float = 28) -> float:
        """Computes the value of the design tensile strength of concrete.

        Arguments:
            t_ref (float, optional): The reference time in days (default is 28
                days).

        Returns:
            float: The design tensile strength of concrete in MPa.

        Raises:
            ValueError: If t_ref is less than 0.
        """
        k_tt = ec2_2023.k_tt(
            t_ref=t_ref, strength_dev_class=self._strength_dev_class
        )
        return ec2_2023.fctd(self.fctk_5, k_tt, self.gamma_c)

    @property
    def eps_c1(self) -> float:
        """Returns the strain at maximum compressive strength of concrete (fcm)
        for the Sargin constitutive law.

        Returns:
            float: The strain at maximum compressive strength of concrete.

        Note:
            The returned value is derived from fcm if eps_c1 is not manually
            provided when initializing the object.
        """
        if self._eps_c1 is None:
            return ec2_2023.eps_c1(self.fcm)
        return self._eps_c1

    @property
    def eps_cu1(self) -> float:
        """Returns the strain at concrete failure of concrete.

        Returns:
            float: The maximum strength at failure of concrete.

        Note:
            The returned value is derived from fcm if eps_cu1 is not manually
            provided when initializing the object.
        """
        if self._eps_cu1 is None:
            return ec2_2023.eps_cu1(self.fcm)
        return self._eps_cu1

    @property
    def k_sargin(self) -> float:
        """Returns the k coefficient for Sargin constitutive law.

        Returns:
            float: k coefficient for Sargin constitutive law.

        Note:
            The returned value is derived from Ecm, fcm, and eps_c1 if k_sargin
            is not manually provided when initializing the object.
        """
        if self._k_sargin is None:
            return ec2_2023.k_sargin(
                Ecm=self.Ecm,
                fcm=self.fcm,
                eps_c1=self.eps_c1,
            )
        return self._k_sargin

    @property
    def eps_c2(self) -> float:
        """Returns the strain at maximum compressive strength of concrete (fcd)
        for the Parabola-rectangle constitutive law.

        Returns:
            float: The strain at maximum compressive strength of concrete.

        Note:
            The returned value is the default value in the design code if not
            manually provided when initializing the object.
        """
        if self._eps_c2 is None:
            return ec2_2023.eps_c2()
        return self._eps_c2

    @property
    def eps_cu2(self) -> float:
        """Returns the strain at concrete failure of concrete for the
        Parabola-rectangle constitutive law.

        Returns:
            float: The maximum strain at failure of concrete.

        Note:
            The returned value is the default value in the design code if not
            manually provided when initializing the object.
        """
        if self._eps_cu2 is None:
            return ec2_2023.eps_cu2()
        return self._eps_cu2

    @property
    def n_parabolic_rectangular(self) -> float:
        """Returns the coefficient for Parabola-rectangle constitutive law.

        Returns:
            float: The exponent for Parabola-recangle law.

        Note:
            The returned value is the default value in the design code if not
            manually provided when initializing the object.
        """
        if self._n_parabolic_rectangular is None:
            return ec2_2023.n_parabolic_rectangular()
        return self._n_parabolic_rectangular

    @property
    def gamma_c(self) -> float:
        """The partial factor for concrete."""
        # Here we should implement the interaction with the globally set
        # national annex. For now, we simply return the default value.
        return self._gamma_c or 1.5

    def __elastic__(self) -> dict:
        """Returns kwargs for creating an elastic constitutive law."""
        return {'E': self.Ecm}

    def __parabolarectangle__(self) -> dict:
        """Returns kwargs for creating a parabola rectangle const law."""
        return {
            'fc': self.fcd(),
            'eps_0': self.eps_c2,
            'eps_u': self.eps_cu2,
            'n': self.n_parabolic_rectangular,
        }

    def __sargin__(self) -> dict:
        """Returns kwargs for creating a Sargin const law."""
        return {
            'fc': self.fcd(),
            'eps_c1': self.eps_c1,
            'eps_cu1': self.eps_cu1,
            'k': self.k_sargin,
        }

    def __popovics__(self) -> dict:
        """Returns kwargs for creating a Sargin const law."""
        return {
            'fc': self.fcd(),
            'eps_c': self.eps_c1,
            'eps_cu': self.eps_cu1,
            'Ec': self.Ecm,
        }
