(versioning)=
# Versioning

## General

StructuralCodes does not yet have a stable API. This does not mean that you should not use the library for solving your design problems, but that there are more features to come, and that some of these new features might be implemented in ways that break existing behaviour.

Until then, StructuralCodes will use a custom versioning scheme that uses the __minor__ version number for breaking changes, and the __patch__ version number for other changes.

## Version changes

A __patch__ version number increase, for example from `0.0.1` to `0.0.2`, occurs when:
- Bugs are fixed.
- New features are added in a backwards compatible manner.

A __minor__ version number increase, for example from `0.0.2` to `0.1.0`, occurs when:
- Changes are introduced in a backwards incompatible manner.

A __major__ version number increase to `1.0.0` indicates that the API is stable, and that [semantic versioning](https://semver.org/) applies from here.
