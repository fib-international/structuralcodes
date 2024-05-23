"""The concrete class for EC2 2023 Concrete Material."""

import typing as t

from structuralcodes.codes import ec2_2023

from ._concrete import Concrete


class ConcreteEC2_2023(Concrete):  # noqa: N801
    """Concrete implementation for EC2 2023 Concrete."""

    # Inherent concrete properties
    _kE: t.Optional[float] = None  # noqa: N815
    _strength_dev_class: t.Optional[str] = None

    # Computed attributes
    _fcm: t.Optional[float] = None
    _fctm: t.Optional[float] = None
    _Ecm: t.Optional[float] = None
    _fctk_5: t.Optional[float] = None
    _fctk_95: t.Optional[float] = None
    _eps_c1: t.Optional[float] = None
    _eps_cu1: t.Optional[float] = None

    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400.0,
        kE: float = 9500,
        strength_dev_class: str = 'CN',
        gamma_c: t.Optional[float] = None,
        existing: bool = False,
    ):
        """Initializes a new instance of Concrete for EC2 2023.

        Args:
            fck (float): Characteristic strength in MPa if concrete is not
                existing.

        Keyword Args:
            name (str): A descriptive name for concrete
            density (float): Density of material in kg/m3 (default: 2400)
            kE (float): coefficient relating aggregates.
            strength_dev_class (str, optional): default is CN.
                Possible values: CS, CN, CR, slow, normal or rapid
            gamma_c (float, optional): partial factor of concrete
                (default is 1.5)
            existing (bool, optional): The material is of an existing structure
                (default: False)
        """
        if name is None:
            name = f'C{round(fck):d}'

        # Check if strength_dev_class is valid
        strength_dev_class = strength_dev_class.strip()
        valid_dev_classes = ('cs', 'cn', 'cr', 'slow', 'normal', 'rapid')
        if strength_dev_class.lower() not in valid_dev_classes:
            raise ValueError(
                'strength_dev_class not valid. '
                + f'{strength_dev_class} not {valid_dev_classes}'
            )

        # Check for valid security coeff
        if gamma_c is not None and gamma_c <= 0:
            raise ValueError(
                f'gamma_c cannot be less than zero. Current: {gamma_c}'
            )

        super().__init__(
            fck=fck,
            name=name,
            density=density,
            existing=existing,
            gamma_c=gamma_c,
        )
        self._kE = kE
        self._strength_dev_class = strength_dev_class

    def _reset_attributes(self):
        """Reset computed properties."""
        # We only need to reset computed properties
        self._fcm = None
        self._fctm = None
        self._Ecm = None
        self._fctk_5 = None
        self._fctk_95 = None
        self._eps_c1 = None
        self._eps_cu1 = None

    @property
    def fcm(self) -> float:
        """Returns the mean strength of concrete.

        Returns:
            float: The mean compressive strength in MPa.
        """
        self._fcm = self._fcm or ec2_2023.fcm(self.fck)
        return self._fcm

    @fcm.setter
    def fcm(self, value: float):
        """Sets a user defined value for the mean strength of concrete.

        Args:
            value (float): the value of the mean strength of
                concrete in MPa.

        Raises:
            ValueError: if value is less or equal than the value of fck
        """
        if abs(value) <= self._fck:
            raise ValueError(
                (
                    'Mean compressive strength cannot be lower than',
                    'characteristic strength.\n',
                    'Current characteristing strength: ',
                    f'fck = {self._fck}.',
                    f'Current value: value = {value}',
                )
            )
        self._fcm = value

    @property
    def fctm(self) -> None:
        """Returns the mean concrete tensile strength.

        Returns:
            float: the mean concrete tensile strength.
        """
        self._fctm = self._fctm or ec2_2023.fctm(self.fck)
        return self._fctm

    @fctm.setter
    def fctm(self, value: float):
        """Sets a custom user defined value for the concrete tensile
            strength for the concrete.

        Args:
            value (float): the new value for fctm in MPa.
        """
        self._fctm = value

    @property
    def fctk_5(self) -> float:
        """Returns the 5% mean concrete tensile strength fractile.

        Returns:
            float: the 5% mean concrete tensile strength fractile in MPa
        """
        self._fctk_5 = self._fctk_5 or ec2_2023.fctk_5(self.fctm)
        return self._fctk_5

    @property
    def fctk_95(self) -> float:
        """Returns the 95% mean concrete tensile strength fractile.

        Returns:
            float: the 5% mean concrete tensile strength fractile in MPa
        """
        self._fctk_95 = self._fctk_95 or ec2_2023.fctk_95(self.fctm)
        return self._fctk_95

    @property
    def Ecm(self) -> float:
        """Returns the secant modulus.

        Returns:
            float: the secant concrete modulus in MPa
        """
        self._Ecm = self._Ecm or ec2_2023.Ecm(self.fcm, self._kE)
        return self._Ecm

    @Ecm.setter
    def Ecm(self, value: float) -> None:
        """Sets the secant modulus.

        Args:
            float: the secand modulus value in MPa.
        """
        self._Ecm = value

    def fcd(
        self, t_ref: float = 28, t0: float = 91, fck_ref: float = 40
    ) -> float:
        """Computes the value of the design compressive strength of
        concrete.

        Args:
            t_ref (float,optional): the reference time in days
                (default is 28 days).
            t0 (float, optional): age at loading in days (default is 91 days).
            fck_ref (float, optional): the reference compressive strength
                MPa (default is 40 MPa).

        Returns:
            float: the design compressive strength of concrete in MPa

        Raises:
            ValueError: if fck_ref is less or equal to 0
            ValueError: if t_ref is less than 0
            ValueError: if t0 is less than 0
        """
        eta_cc = ec2_2023.eta_cc(self.fck, fck_ref=fck_ref)
        k_tc = ec2_2023.k_tc(t_ref, t0, self._strength_dev_class)
        return ec2_2023.fcd(self.fck, eta_cc, k_tc, self.gamma_c)

    def fctd(self, t_ref: float = 28) -> float:
        """Computes the value of the design tensile strength of concrete.

        Args:
            t_ref (float, optional): the reference time in days
                (default is 28 days).

        Returns:
            float: the design tensile strength of concrete in MPa

        Raises:
            ValueError: if t_ref is less than 0
        """
        k_tt = ec2_2023.k_tt(
            t_ref=t_ref, strength_dev_class=self._strength_dev_class
        )
        return ec2_2023.fctd(self.fctk_5, k_tt, self.gamma_c)

    @property
    def eps_c1(self) -> float:
        """Returns the strain at maximum compressive strength
            of concrete (fcm) for the Sargin constitutive law.

        Returns:
            float: the strain at maximum compressive strength of concrete
        """
        self._eps_c1 = self._eps_c1 or ec2_2023.eps_c1(self.fcm)
        return self._eps_c1

    @property
    def eps_cu1(self) -> float:
        """Returns the strain at concrete failure of concrete.

        Returns:
            float: the maximum strength at failure of concrete
        """
        self._eps_cu1 = self._eps_cu1 or ec2_2023.eps_cu1(self.fcm)
        return self._eps_cu1

    def sigma_c(self, eps_c: float) -> float:
        """Computes the compressive stress of concrete given a
        strain eps_c under short term uniaxial compression.

        Args:
            eps_c (float): the strain of concrete

        Returns:
            float: the compressive stress of concrete in MPa

        Raises:
            ValueError: if eps_c is less than 0
        """
        return ec2_2023.sigma_c(
            self.Ecm, self.fcm, eps_c, self.eps_c1, self.eps_cu1
        )

    @property
    def gamma_c(self) -> float:
        """The partial factor for concrete."""
        # Here we should implement the interaction with the globally set
        # national annex. For now, we simply return the default value.
        return self._gamma_c or 1.5
