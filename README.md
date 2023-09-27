# Structural Design Codes

*A Python package with models from structural design codes.*

[![Build](https://github.com/fib-international/structuralcodes/actions/workflows/build.yaml/badge.svg)](https://github.com/fib-international/structuralcodes/actions/workflows/build.yaml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**:construction: Note! This package is work in progress. Feel free to join the discussions in the pull requests. :construction:**

## How to install :grey_question:

The package is still on a premature level. As soon as we start tagging versions, the package will be pip-installable from source. Eventually, the package will be published on PyPI.

## Package structure :card_file_box:

The package will be based on a hierachical structure. On the lowest level, models from different structural design codes are implemented in sub-packages. The work is coordinated by the _fib_, and the _fib_ Model Code and Eurocode 2 are used as cases for establishing the structure.

On the higher levels there will be material, cross section and element classes related to specific design codes, e.g. `ConcreteMC2010` and `RCRectangularMC2010`.

## Contribution guidelines :handshake:

For the moment, the core team has completed an initial structure for the two lowest levels, i.e. how to structure the models from the structural design codes, and how to create classes that represent materials. Contributions that seek to add content to the lowest level are highly appreciated :heart:

We have defined a set of [milestones](https://github.com/fib-international/structuralcodes/milestones) that indicate our priorities. Issues are linked to each milesstone, and some of the issues include a to-do list and have linked pull requests. Contributions related to the milestones are highly appreciated. For the moment, our first priority is to add code related to _fib_ Model Code 2020, but other contributions are of course also welcome!

Here is how to get started.

First:

- Take a look at the milestones, find one that suits your interests, and navigate to a related issue. Create a comment in the issue indicating the topic that you wish to contribute. Make sure that the topic is not already started by someone else.
- If you didn't find a suitable milestone or issue, please create a new issue where you describe your request or feature idea. Make sure that you get a GO from the core team before continuing your contribution. We promise to provide feedback as soon as we can.

Then:
- Create a fork of the repo.
- Set up your local development environment, see below :arrow_heading_down:.
- Start developing the feature in a feature branch with a suitable title :rocket:
- Please follow common best practice when committing, provide a meaningful commit message, and split different sub-topics in separate commits :thumbsup:
- Feel free to open a draft pull request as soon as you start committing so that the community can follow the development. The core team will make sure the pull request is linked to the relevant issue, and that we [keep track of the progress](https://github.com/orgs/fib-international/projects/1).
- We use `black`, `flake8` and `pylint` for code formatting and linting, and we have a test suite based on `pytest`. Please make sure that these checks pass before committing code :heavy_check_mark:
- Remember to add tests for the new features that you implement :man_scientist:
- Remember to add docstrings to your functions and classes. We use the [google format](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).
- When you feel that your contribution is ready, mark your pull request as ready for review, and we will assign some suitable reviewers :nerd_face:

The core team continues working on the structure for the higher levels of the package. Feel free to join the discussions in the pull requests :raising_hand:

## Set up your local development environment :woman_technologist:

For the moment, there are no third party dependencies to the package, but it is always useful to work in a virtual environment. To create one on Windows, type `py -m venv venv`, and activate it by typing `venv\scripts\activate`.

We use `black` for code formatting, `flake8` and `pylint` for linting, and have set up our test suite using `pytest`. Make sure you have these installed (`make deps`). We have set up VS Code to use these tools automatically. Open VS Code inside the repo (`code .`), view the tests in the test pane on the left side, and notice how the linters highlight parts of your code that should be rewritten.

We have supplied a makefile that you may use to run formatting (`make form`), linting (`make lint`) and testing (`make test`). Don't have `make`? You can get it from [chocolatey](https://community.chocolatey.org/packages/make), `choco install make`. If you don't want to use `make`, feel free to e.g. `py -m black structuralcodes` and `py -m pytest`.
