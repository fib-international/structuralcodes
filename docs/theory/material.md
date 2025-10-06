# Material

In Structuralcodes, materials are objects that serve as interfcaces between material properties according to Design Codes and the Sections. This includes constitutive behavior, that may be uniaxial stress-strain laws (useful for beam section analyses) or 2D constitutive laws (useful for shell section analysis).
In structural codes we have a library of Materials that can be already used. Furthermore users can define their own materials, provided that they adhere to the interface specified by base [`Material`](api-base-material-class) class.

Available materials in Structuralcodes can be grouped according to the following:
1. Base materials
2. Concrete materials
3. Reinforcement steel materials.

## Base materials

Base materials include some standard material for general-purpose analysis like `ElasticMaterial`, `ElasticPlasticMaterial`, 