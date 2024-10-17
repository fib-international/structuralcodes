(contributing)=
# Contributing

[![Code style: Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/format.json)](https://github.com/charliermarsh/ruff)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/charliermarsh/ruff)

[![Tests](https://github.com/fib-international/structuralcodes/actions/workflows/build.yaml/badge.svg)](https://github.com/fib-international/structuralcodes/actions/workflows/build.yaml)
[![Publish](https://github.com/fib-international/structuralcodes/actions/workflows/create-release.yml/badge.svg)](https://github.com/fib-international/structuralcodes/actions/workflows/create-release.yml)
[![Docs](https://github.com/fib-international/structuralcodes/actions/workflows/sphinx.yml/badge.svg)](https://github.com/fib-international/structuralcodes/actions/workflows/sphinx.yml)

## Contribution guidelines

Here is how to get started.

_First_

- Take a look at the [milestones](https://github.com/fib-international/structuralcodes/milestones), find one that suits your interests, and navigate to a related issue. Create a comment in the issue indicating the topic that you wish to contribute. Make sure that the topic is not already started by someone else.
- If you didn't find a suitable milestone or issue, please create a new issue where you describe your request or feature idea. Make sure that you get a GO from the core team before continuing your contribution, so that you do not waste your time duplicating the work of someone else. We promise to provide feedback as soon as we can.

_Then_

- Create a fork of the repo.
- Set up your local development environment, see [below](local-dev-env).
- Start developing the feature in a feature branch with a suitable title.
- We use `ruff` for code formatting and linting, and we have a test suite based on `pytest`. Please make sure that these checks pass before committing code.
- Please follow common best practice when committing, provide a meaningful commit message, and split different sub-topics in separate commits.
- Feel free to open a draft pull request as soon as you start committing so that the community can follow the development. The core team will make sure the pull request is linked to the relevant issue, and that we [keep track of the progress](https://github.com/orgs/fib-international/projects/1).
- Remember to add tests for the new features that you implement.
- Remember to add docstrings to your functions and classes. We use the [google format](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).
- When you feel that your contribution is ready, mark your pull request as ready for review, and we will assign some suitable reviewers.

The core team continues working on the structure for the higher levels of the package. Feel free to join the discussions in the pull requests.

(local-dev-env)=
## Set up your local development environment

To keep track of the library dependencies, we recommend to work in a virtual environment. To create one on Windows, type `py -m venv venv`, and activate it by typing `venv\scripts\activate`.

We use `ruff` for linting and code formatting, and have set up our test suite using `pytest`. Make sure you have these installed (`make deps`). We have set up VS Code to use testing. Open VS Code inside the repo (`code .`) and view the tests in the test pane on the left side.

We are using the [Ruff extension for VS Code](https://github.com/astral-sh/ruff-vscode?tab=readme-ov-file), and we have by default enabled all automatic fixes. If you wish to disable these fixes locally, e.g. auto-format or auto-sort imports on save, you should have a look at `.vscode/settings.json`.

We have supplied a makefile that you may use to run formatting (`make form`), linting (`make lint`) and testing (`make test`). Don't have `make`? You can get it from [chocolatey](https://community.chocolatey.org/packages/make), `choco install make`. If you don't want to use `make`, feel free to manually type the commands from the makefile e.g. `py -m ruff format structuralcodes` and `py -m pytest`.