"""Abstract base classes"""
import abc


class Material(abc.ABC):
    """Abstract base class for materials."""


class DesignCode(abc.ABC):
    """Abstract base class for design codes."""


class ConcreteDesignCode(DesignCode):
    """Abstract class for concrete design codes. 
    (e.g. MC2010, EC2, EC8-3 for existing buildings etc.)"""