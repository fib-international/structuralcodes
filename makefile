.DEFAULT_GOAL := deps
.PHONY: deps form lint

PACKAGE_NAME = structuralcodes

deps:  ## Install dependencies
	python -m pip install --upgrade pip
	python -m pip install --upgrade wheel
	python -m pip install -r requirements.txt
	python -m pip install --upgrade black
	python -m pip install --upgrade flake8 mccabe pylint
	python -m pip install --upgrade flit
	python -m pip install --upgrade pytest pytest-cov

form:  ## Code formatting
	python -m black $(PACKAGE_NAME)

lint:  ## Linting and static type checking
	python -m flake8 $(PACKAGE_NAME)
	python -m pylint $(PACKAGE_NAME)

test:  ## Run tests and output reports
	python -m pytest --junitxml=junit/test-results.xml --cov=$(PACKAGE_NAME) --cov-report=term-missing --cov-report=xml
