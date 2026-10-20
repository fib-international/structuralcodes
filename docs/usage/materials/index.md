(usage-materials)=
# Materials

## General

The [material classes](#api-materials) in StructuralCodes are the natural starting point when building any structural design workflow. The material classes are specific to each design code, and they contain attributes that represent common properties from the design code. A material object also contains a [constitutive law](#api-constitutive-laws). The constitutive law is accessed from the [section calculator](#api-section-calculator) when integrating the strain response in a section.

StructuralCodes contains material classes for representing [concrete](#usage-concrete-materials) and [reinforcement](#usage-reinforcement-materials), in addition to a collection of [basic materials](#usage-basic-materials) for representing arbitrary materials with given constitutive laws.

To create concrete and reinforcement objects, we can either directly import specific classes as shown in the sections below, or we can set an active design code and use the factory functions.

(code-usage-materials-factory)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_materials_factory.py
:lines: 3-22
:caption: Create material objects by setting a design code and calling factory functions.
:::
::::

Notice how the factory functions accept the same keyword arguments as the constructors of the specific classes. This allows us to pass e.g. `gamma_c` and `alpha_cc` to `create_concrete`. Also notice how we are using {func}`reinforcement_duct_props() <structuralcodes.codes.ec2_2004.reinforcement_duct_props>` to get a dictionary with properties related to the specified ductility class of reinforcement.

After creating a material object, its properties are easily accessed, e.g. `concrete.Ecm`. Most notably, the `.constitutive_law` attribute reveals the [constitutive law](#usage-constitutive-laws) of the material.

:::{tip}

All material classes are subclasses of the {class}`Material <structuralcodes.core.base.Material>` base class, and all constitutive laws are subclasses of the {class}`ConstitutiveLaw <structuralcodes.core.base.ConstitutiveLaw>` base class. These base classes define interfaces. This means that we can create our own custom classes as subclasses of these base classes, and as long as we implement the necessary methods defined by the interface, we can swap the builtin classes with our custom ones. This opens for a lot of flexibility inside StructuralCodes.
:::

(usage-materials-are-immutable)=
:::{attention}

The material objects in StructuralCodes are immutable. This means that as soon as a material object is created, the attributes of the object cannot be changed.
:::

(usage-concrete-materials)=
## Concrete materials

There are classes for representing concrete according to the [current version of Eurocode 2](#api-concrete-ec2-2004), [_fib_ Model Code 2010](#api-concrete-mc2010), and the [next generation of Eurocode 2](#api-concrete-ec2-2023).

Import {class}`ConcreteEC2_2004 <structuralcodes.materials.concrete.ConcreteEC2_2004>` to create an object representing concrete according to Eurocode 2 (2004). 

(code-usage-concrete-ec2-2004)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_specific_concrete.py
   :lines: 3, 5-6
   :caption: Create a concrete object according to Eurocode 2 (2004).
:::
::::

Notice how it is possible to override the default values of the properties by passing these to the constructor.

(code-usage-concrete-ec2-2004-override-defaults)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_specific_concrete.py
:lines: 8
:caption: Override default values in a concrete object.
:::
::::

To override the constitutive law, we can also pass a [constitutive law object](#api-constitutive-laws) valid for concrete.

(code-usage-concrete-ec2-2004-override-constitutive-law)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_specific_concrete.py
:lines: 3-4, 9-11
:caption: Override the constitutive law in a concrete object.
:::
::::

(usage-reinforcement-materials)=
## Reinforcement materials

Reinforcement objects are created and used in a similar manner as the concrete objects, and there are reinforcement classes according to the [current version of Eurocode 2](#api-reinforcement-ec2-2004), [_fib_ Model Code 2010](#api-reinforcement-mc2010), and the [next generation of Eurocode 2](#api-reinforcement-ec2-2023).

Import {class}`ReinforcementEC2_2004 <structuralcodes.materials.reinforcement.ReinforcementEC2_2004>` to create a reinforcement object according to Eurocode 2 (2004).

(code-usage-reinforcement-ec2-2004)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_specific_reinforcement.py
:lines: 3-5
:caption: Create a reinforcement object according to Eurocode 2 (2004).
:::
::::

(usage-constitutive-laws)=
## Constitutive laws

When initializing a material object it will include a [constitutive law](#api-constitutive-laws) object as well. The purpose of the constitutive law object is to provide a relation between stresses and strains, or modulus and strain. Therefore all constitutive law objects implement the `.get_stress()`, `.get_tangent()`, and `.get_secant()` methods. As shown above, the [concrete](#usage-concrete-materials) and [reinforcement](#usage-reinforcement-materials) classes define default constitutive laws, that can be overridden by the constructor.

The example below shows how to access the ultimate strains of the constitutive law, create an array of linearly spaced strain values, and passing these strain values to the `.get_stress()` method to get corresponding stresses. Note that only compressive strains and stresses are computed in this example. In this example we are passing an array of strain values to `.get_stress()`, but the method also accepts a scalar argument.

(code-usage-compute-stress)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_compute_stress_strain.py
:lines: 4-14
:caption: Compute stresses with a constitutive law.
:::
::::

:::::{tip}

Use your favourite plotting library to visualize the computed stresses. The code below shows how to produce the figure [below](#fig-usage-parabola-rectangle) with [Matplotlib](https://matplotlib.org/).

(code-usage-visualize-stress-strain)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_compute_stress_strain.py
:lines: 3, 15-24
:caption: Visualize stress-strain relation.
:::
::::

(fig-usage-parabola-rectangle)=
:::{figure} parabola_rectangle.png
:width: 75%

The parabola-rectangle distribution computed with the code [above](#code-usage-compute-stress).
:::

:::::

To create a constitutive law, you can either use the [constitutive law factory](#api-constitutive-law-factory), or directly import a specific constitutive law class as shown in the [example above](#code-usage-concrete-ec2-2004-override-constitutive-law). You can call {func}`get_constitutive_laws_list() <structuralcodes.materials.constitutive_laws.get_constitutive_laws_list>` to get a list of all the available constitutive laws.

The constitutive law factory takes a name and a material object as arguments, and returns a constitutive law object. For it to work, the material object which is used as argument needs to implement a dunder method related to the specified constitutive law, in order to provide the relevant parameters for the constitutive law.

The example below shows how we can use the constitutive law factory to create a {class}`BilinearCompression <structuralcodes.materials.constitutive_laws.BilinearCompression>` object based on a {class}`ConcreteEC2_2004 <structuralcodes.materials.concrete.ConcreteEC2_2004>` object. Since the {class}`ConcreteEC2_2004 <structuralcodes.materials.concrete.ConcreteEC2_2004>` object implements the {meth}`__bilinearcompression__ <structuralcodes.materials.concrete.ConcreteEC2_2004.__bilinearcompression__>` method, the constitutive law factory can call this method to get a dictionary with the relevant parameters to create the constitutive law. This workflow can be used if you wish to compare several constitutive laws that are based on the same material.

(code-usage-constitutive-law-factory)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_constitutive_law_factory.py
:lines: 3-
:caption: Create a constitutive law using the constitutive law factory.
:::
::::

:::{admonition} Customize the material behaviour
:class: tip

There are a few ways to customize constitutive laws in StructuralCodes.

The first is to simply provide user values to the parameters of the available constitutive laws.

The second is to use the special {class}`UserDefined <structuralcodes.materials.constitutive_laws.UserDefined>` constitutive laws which allows you to provide two arrays in the constructor: one for strains and one for stresses, i.e. representing any arbitrary uniaxial relation between stresses and strains.

The third is to create a new subclass of {class}`ConstitutiveLaw <structuralcodes.core.base.ConstitutiveLaw>` that represents your custom constitutive law. As long as you implement the `.get_stress()`, `.get_tangent()`, and `.get_ultimate_strain()` methods, you can use it in a material as basis for a geometry in a section calculation.

:::

(usage-basic-materials)=
## Basic materials

To simplify creating materials with simple constitutive laws, StructuralCodes contains a set of [basic materials](#api-basic-materials), including {class}`ElasticMaterial <structuralcodes.materials.basic.ElasticMaterial>`, {class}`ElasticPlasticMaterial <structuralcodes.materials.basic.ElasticPlasticMaterial>`, and {class}`GenericMaterial <structuralcodes.materials.basic.GenericMaterial>`.

The {class}`GenericMaterial <structuralcodes.materials.basic.GenericMaterial>` is a material class that accepts any constitutive law.

The {class}`ElasticMaterial <structuralcodes.materials.basic.ElasticMaterial>` and the {class}`ElasticPlasticMaterial <structuralcodes.materials.basic.ElasticPlasticMaterial>` are materials that are initialized with just the necessary data to also initialize an {class}`Elastic <structuralcodes.materials.constitutive_laws.Elastic>` and an {class}`ElasticPlastic <structuralcodes.materials.constitutive_laws.ElasticPlastic>` constitutive law, respectively. The example below shows how to create an {class}`ElasticPlasticMaterial <structuralcodes.materials.basic.ElasticPlasticMaterial>`.

(code-usage-create-elastic-plastic-material)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_elasticplastic_material.py
:lines: 3-
:caption: Create a material with an elastic-plastic constitutive law.
:::
::::

Note how you can use {func}`.from_material() <structuralcodes.materials.basic.ElasticMaterial.from_material>` on {class}`ElasticMaterial <structuralcodes.materials.basic.ElasticMaterial>` to create an elastic representation of any other material on the fly. This method makes use of the constitutive law of the other material to get the tangent modulus at zero strain, and passes this value to the constructor of the elastic material. See the example code below. This feature is used for example in {func}`calculate_elastic_cracked_properties() <structuralcodes.sections.calculate_elastic_cracked_properties>` to represent the materials in a reinforced concrete geometry with elastic materials.

(code-usage-create-elastic-material-from-other)=
::::{dropdown-syntax}
:::{literalinclude} ../../_example_code/usage_create_elastic_material_from_other.py
:lines: 3-
:caption: Create a material with an elastic constitutive law based on another material.
:::
::::