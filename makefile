.DEFAULT_GOAL := deps
.PHONY: deps form lint

PACKAGE_NAME = structuralcodes

deps:  ## Install dependencies
	python -m pip install --upgrade pip
	python -m pip install --upgrade wheel
	python -m pip install -r requirements.txt
	python -m pip install ruff==0.3.3
	python -m pip install --upgrade flit
	python -m pip install --upgrade pytest pytest-cov

form:  ## Code formatting
	python -m ruff format $(PACKAGE_NAME)
	python -m ruff format tests

lint:  ## Linting and static type checking
	python -m ruff check $(PACKAGE_NAME)
	python -m ruff check tests

test:  ## Run tests and output reports
	python -m pytest --junitxml=junit/test-results.xml --cov=$(PACKAGE_NAME) --cov-report=term-missing --cov-report=xml