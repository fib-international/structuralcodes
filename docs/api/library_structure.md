(api-structure)=
# Library structure

The content of StructuralCodes is organized in a hierarchy:

1. On the lowest level, we have the {ref}`code equations <api-design-codes>` from the various design codes. For example, the code equations from Eurocode 2 (2004) can be imported by typing `from structuralcodes.codes import ec2_2004`.
2. Next, we have the classes that represent {ref}`materials <api-materials>` and {ref}`constitutive laws <api-constitutive-laws>` from the design codes. For example, import concrete from Eurocode 2 (2004) by typing `from structuralcodes.materials.concrete import ConcreteEC2_2004`. The materials in turn refer to code equations from the lower level to calculate material properties.
3. Next we have the classes that represent {ref}`sections <api-sections>` and {ref}`geometry <api-geometry>`. We rely on the [Shapely](https://shapely.readthedocs.io/en/stable/) library for creating geometries. The sections in turn refer to materials and code equations from the lower levels. Section calculators are responsible for calculating the response in sections, and section integrators are responsible for integrating stress resultants from strains.

:::{seealso}
{ref}`Design codes reference <api-design-codes>`.

{ref}`Materials reference <api-materials>`

{ref}`Geometry reference <api-geometry>`

{ref}`Sections reference <api-sections>`

The [Shapely](https://shapely.readthedocs.io/en/stable/) library.

:::