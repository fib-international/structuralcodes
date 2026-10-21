from shapely import Polygon

from structuralcodes import set_design_code
from structuralcodes.core.base import ConstitutiveLaw
from structuralcodes.geometry import SurfaceGeometry, add_reinforcement
from structuralcodes.materials.concrete import create_concrete
from structuralcodes.materials.constitutive_laws import (
    ElasticFragileTension,
    ParabolaRectangle,
    Parallel,
)
from structuralcodes.materials.reinforcement import create_reinforcement
from structuralcodes.sections import BeamSection

set_design_code('ec2_2004')

# Parametres beton et armatures (repris de quickstart.py)
fck = 45
Ec = 36000
fct = 2.7
fyk = 500
ftk = 550
Es = 200000
epsuk = 0.07

# Lois constitutives : compression seule, traction seule, et les deux (Parallel)
compression: ConstitutiveLaw = ParabolaRectangle(fc=fck)
tension: ConstitutiveLaw = ElasticFragileTension(fct=fct, Ec=Ec)
compression_tension: ConstitutiveLaw = Parallel(
    constitutive_laws=[compression, tension]
)

concrete_c = create_concrete(fck=fck, constitutive_law=compression)
concrete_t = create_concrete(fck=fck, constitutive_law=tension)
concrete_ct = create_concrete(fck=fck, constitutive_law=compression_tension)

reinforcement = create_reinforcement(fyk=fyk, Es=Es, ftk=ftk, epsuk=epsuk)

# Geometrie rectangulaire avec armatures (reprise de quickstart.py)
width = 250
height = 500
diameter_reinf = 25
cover = 50


def build_section(
    concrete, reinforced: bool, integrator: str = 'fiber'
) -> BeamSection:
    """Construit une section rectangulaire pour un beton donne.

    Reprend la geometrie de quickstart.py : section 250x500. Si
    `reinforced` est vrai, ajoute deux barres de diametre 25 en bas,
    positionnees avec un enrobage de 50.
    """
    polygon = Polygon(
        [
            (-width / 2, -height / 2),
            (width / 2, -height / 2),
            (width / 2, height / 2),
            (-width / 2, height / 2),
        ]
    )
    geometry = SurfaceGeometry(poly=polygon, material=concrete)
    if reinforced:
        geometry = add_reinforcement(
            geometry,
            (
                -width / 2 + cover + diameter_reinf / 2,
                -height / 2 + cover + diameter_reinf / 2,
            ),
            diameter_reinf,
            reinforcement,
        )
        geometry = add_reinforcement(
            geometry,
            (
                width / 2 - cover - diameter_reinf / 2,
                -height / 2 + cover + diameter_reinf / 2,
            ),
            diameter_reinf,
            reinforcement,
        )
    return BeamSection(geometry, integrator=integrator)


def m_y(concrete, reinforced: bool, integrator: str) -> float:
    """Retourne le moment resistant m_y pour un beton, un ferraillage et un integrateur donnes."""
    section = build_section(
        concrete, reinforced=reinforced, integrator=integrator
    )
    return section.section_calculator.calculate_bending_strength().m_y


concretes = {
    'C  (compression only)    ': concrete_c,
    'T  (tension only)        ': concrete_t,
    'CT (compression+tension) ': concrete_ct,
}


def diff_pct(a: float, b: float) -> float:
    """Ecart relatif de b par rapport a a, en pourcent (0 si a est nul)."""
    return 0.0 if a == 0 else (b - a) / a * 100


print('Results (m_y), comparaison fiber vs marin :')
for reinforced in (True, False):
    etat = 'avec armatures' if reinforced else 'sans armatures'
    print(f'--- {etat} ---')
    for label, concrete in concretes.items():
        fiber = m_y(concrete, reinforced, 'fiber')
        marin = m_y(concrete, reinforced, 'marin')
        print(
            f'Section {label} fiber: {fiber:>20.6f}  marin: {marin:>20.6f}'
            f'  diff: {diff_pct(fiber, marin):+.3f} %'
        )
