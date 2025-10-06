# Material

In Structuralcodes, materials are objects that serve as interfcaces between material properties according to Design Codes and the Sections. This includes constitutive behavior, that may be uniaxial stress-strain laws (useful for beam section analyses) or 2D constitutive laws (useful for shell section analysis).
In structural codes we have a library of Materials that can be already used. Furthermore users can define their own materials, provided that they adhere to the interface specified by base [`Material`](api-base-material-class) class.

Available materials in Structuralcodes can be grouped according to the following:
1. [Basic materials](theory-basic-materials)
2. [Concrete materials](theory-concrete-materials)
3. [Reinforcement steel materials](theory-reinforcement-materials).

(theory-basic-materials)=
## Basic materials

Base materials include some standard materials for general-purpose analysis like [`ElasticMaterial`](api-elastic-material), [`ElasticPlasticMaterial`](api-elastic-plastic-material) or [`GenericMaterial`](api-generic-material).

[`ElasticMaterial`](api-elastic-material) represents a simple elastic constitutive law and is defined by the elastic modulus.

[`ElasticPlasticMaterial`](api-elastic-plastic-material) represents a material characterized by an elastic-plastic behavior, optionally with an hardening branch.

[`GenericMaterial`](api-generic-material) permits to define a custom constitutive law characterized by multi-linear branch law.


(theory-concrete-materials)=
## Concrete materials

Concrete materials is a collection of several materials representing concrete according to different Design codes (e.g. Model Code 2010, EC2 2004, EC2 2023). For more information refer to the [api](api-concrete-materials).

(theory-reinforcement-materials)=
## Reinforcement steel materials

Reinforcement steel materials is a collection of several materials representing reinforcement steel according to different Design codes (e.g. Model Code 2010, EC2 2004, EC2 2023). For more information refer to the [api](api-reinforcement-steel-materials).

