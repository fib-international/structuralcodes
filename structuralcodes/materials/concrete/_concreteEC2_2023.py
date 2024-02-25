"""The concrete class for EC2 2023 Concrete Material."""

import typing as t

from structuralcodes.codes import ec2_2023
from ._concrete import Concrete


# TODO: complete the documentation and discuss if properties/functions
# are considered proper for the framework
class ConcreteEC2_2023(Concrete):  # noqa: N801
    """Concrete implementation for EC2 2023 Concrete."""

    # Inherent concrete properties
    _kE: t.Optional[float] = None  # noqa: N815
    _concrete_dev_class: t.Optional[str] = None
    _gamma_C: t.Optional[float] = None  # noqa: N815

    # Computed attributes
    _fcm: t.Optional[float] = None
    _fctm: t.Optional[float] = None
    _Ecm: t.Optional[float] = None

    # TODO: do we agree with the optional keyword arguments here?
    def __init__(
        self,
        fck: float,
        name: t.Optional[str] = None,
        density: float = 2400.0,
        kE: float = 9500,
        concrete_dev_class: str = 'CS',
        gamma_C: float = 1.5,
        existing: bool = False,
    ):
        """Initializes a new instance of Concrete for MC 2010.

        Args:
            fck (float): Characteristic strength in MPa if concrete is not
                existing.

        Keyword Args:
            name (str): A descriptive name for concrete
            density (float): Density of material in kg/m3 (default: 2400)
            kE (float): coefficient relating aggregates.
            concrete_dev_class (str, optional): default is CS
            gamma_C (float, optioanl): partial factor of concrete
                (default is 1.5)
            existing (bool, optional): The material is of an existing structure
                (default: False)
        """
        if name is None:
            name = f'C{round(fck):d}'
        super().__init__(
            fck=fck, name=name, density=density, existing=existing
        )
        self._kE = kE
        self._concrete_dev_class = concrete_dev_class
        self._gamma_C = gamma_C

        # TODO: implement this method
        def _reset_attributes(self):
            # We only need to reset computed properties
            self._fcm = None
            self._fctm = None
            self._Ecm = None

        @property
        def fcm(self) -> float:
            """Returns the mean strength of concrete.

            Returns:
                float: The mean compressive strength in MPa.
            """
            if self._fcm is not None:
                return self._fcm
            self._fcm = ec2_2023.fcm(self.fck)
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
            self._reset_attributes()
            self._fcm = value

        @property
        def fctm(self) -> None:
            """Returns the mean concrete tensile strength.

            Returns:
                float: the mean concrete tensile strength.
            """
            if self._fctm is not None:
                return self._fctm
            self._fctm = ec2_2023.fctm(self.fck)
            return self._fctm

        @fctm.setter
        def fctm(self, value: float):
            """Sets a custom user defined value for the concrete tensile
                strength for the concrete.

            Args:
                value (float): the new value for fctm in MPa.
            """
            self._reset_attributes()
            self._fctm = value

        @property
        def fctk_5(self) -> float:
            return ec2_2023.fctk_5(self.fctm)

        @property
        def fctk_95(self) -> float:
            return ec2_2023.fctk_95(self.fctm)

        @property
        def Ecm(self) -> float:
            if self._Ecm is not None:
                return self._Ecm
            self._Ecm = ec2_2023.Ecm(self.fcm, self._kE)
            return self._Ecm

        @Ecm.setter
        def Ecm(self, value: float) -> None:
            self._reset_attributes()
            self._Ecm = value

        def fcd(self, t_ref: float, t0: float, fck_ref: float = 40) -> float:
            """Computes the value of the design compressive strength of
            concrete.

            Args:
                t_ref (float): the reference time in days
                t0 (float): age at loading in days
                fck_ref (float, optional): the reference compressive strength
                    MPa (default is 40 MPa).

            Returns:
                float: the design compressive strength of concrete in MPa

            Raises:
                ValueError: if fkc_ref is less or equal to 0
                ValueError: if t_ref is less than 0
                ValueError: if t0 is less than 0
            """
            eta_cc = ec2_2023.eta_cc(self.fck, fck_ref=fck_ref)
            k_tc = ec2_2023.k_tc(t_ref, t0, self.concrete_dev_class)
            return ec2_2023.fcd(self.fck, eta_cc, k_tc, self.gamma_C)

        def fctd(self, t_ref: float) -> float:
            """Computes the value of the design tensile strength of concrete.

            Args:
                t_ref (float): the reference time in days

            Returns:
                float: the design tensile strength of concrete in MPa

            Raises:
                ValueError: if t_ref is less than 0
            """
            k_tt = ec2_2023.k_tt(t_ref=t_ref)
            return ec2_2023.fctd(self.fctk_5, k_tt, self.gamma_C)

        @property
        def eps_c1(self) -> float:
            return ec2_2023.eps_c1(self.fcm)

        @property
        def eps_cu1(self) -> float:
            return ec2_2023.eps_cu1(self.fcm)

        @property
        def sigma_c(self) -> float:
            return ec2_2023.sigma_c(
                self.Ecm, self.fcm, self.eps_c, self.eps_c1, self.eps_cu1
            )
