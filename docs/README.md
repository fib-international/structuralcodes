# About the structuralcodes docs

## General

The docs are built using [sphinx](https://pypi.org/project/Sphinx/). We rely on [sphinx-design](https://pypi.org/project/sphinx_design/) and [furo](https://pypi.org/project/furo/) for the design. The content is written in MyST Markdown using the [MyST-Parser](https://pypi.org/project/myst-parser/4.0.0/).

## How to build

Navigate to the docs folder.

```
cd docs
```

If you have not created the virtual environment, you should do that. Note that the docs require a different virtual environment than the main library, such that there will be one venv-folder in the root directory and one venv-folder in the docs directory.

```
py -m venv venv
```

Activate the virtual environment.

```
venv\scripts\activate
```

Install the dependencies.

```
py -m pip install -r requirements.txt
```

The docs are built using [sphinx-autobuild](https://pypi.org/project/sphinx-autobuild/). Type the following to build the docs, and host a test server locally which automatically rebuilds the docs and reloads the browser when you make changes. By default, the docs are available in your browser on `http://127.0.0.1:8000`, but the port might be different from `8000` if you have other instances running. Make sure to review the output in the terminal.

```
sphinx-autobuild . _build
```

You can stop the test server by hitting `ctrl+c`.

Sometimes, you might need to force a complete rebuild of the docs. You accomplish this through the following steps:

1. Stop the test server if it is running, i.e. hit `ctrl+c` in the terminal.
2. Delete the `_build`-folder, i.e. typing `rmdir .\_build\` in powershell.
3. Restart the test server with `sphinx-autobuild` as described above.
