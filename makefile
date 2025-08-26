.DEFAULT_GOAL := deps
.PHONY: deps form lint test

PACKAGE_NAME = structuralcodes

deps:  ## Install dependencies
	python -m pip install --upgrade uv
	python -m uv pip install -r requirements.txt
	python -m uv pip install ruff==0.12.10
	python -m uv pip install --upgrade flit
	python -m uv pip install --upgrade pytest pytest-cov

form:  ## Code formatting
	python -m ruff format $(PACKAGE_NAME)
	python -m ruff format tests

lint:  ## Linting and static type checking
	python -m ruff check $(PACKAGE_NAME)
	python -m ruff check tests

test:  ## Run tests and output reports
	python -m pytest --junitxml=junit/test-results.xml --cov=$(PACKAGE_NAME) --cov-report=term-missing --cov-report=xml