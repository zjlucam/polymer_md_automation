from dataclasses import dataclass
from enum import Enum

AVOGADROS_NUMBER = 6.022e23
CM3_TO_NM3 = 1e21
NM3_TO_CM3 = 1e-21
NM3_TO_M3 = 1e-27
ANGSTROM_TO_CM = 1e-8
ANGSTROM3_TO_M3 = 1e-30
KG_TO_G = 1e3
DENSITY_TOLERANCE_PERCENTAGE = 5.0
DEFAULT_BOX_SIZE_NM = [5, 5, 5]
ATOMISTIC_TRUNCATED_POLYMER_LENGTH = 5


class LengthUnits(Enum):
    ANGSTROM: float = 1e-10
    NANOMETER: float = 1e-9
    MICROMETER: float = 1e-6
    MILLIMETER: float = 1e-3
    CENTIMETER: float = 1e-2
    METER: float = 1e0
    KILOMETER: float = 1e3


CONVERSION_FACTORS_TO_M = {
    LengthUnits.ANGSTROM: 1e-10,
    LengthUnits.NANOMETER: 1e-9,
    LengthUnits.MICROMETER: 1e-6,
    LengthUnits.MILLIMETER: 1e-3,
    LengthUnits.CENTIMETER: 1e-2,
    LengthUnits.METER: 1,
    LengthUnits.KILOMETER: 1e3,
}


class MassUnits(Enum):
    AMU: float = 1.66053906660e-27
    GRAM: float = 1e-3
    KILOGRAM: float = 1


CONVERSION_FACTORS_TO_KG = {
    MassUnits.AMU: 1.66053906660e-27,
    MassUnits.GRAM: 1e-3,
    MassUnits.KILOGRAM: 1,
}
