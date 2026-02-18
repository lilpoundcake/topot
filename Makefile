.PHONY: help install install-dev clean test lint format type build docs

help:
	@echo "TOPOT Development Commands"
	@echo ""
	@echo "Installation:"
	@echo "  make install           Install package"
	@echo "  make install-dev       Install with development dependencies"
	@echo ""
	@echo "Testing & Quality:"
	@echo "  make test              Run tests"
	@echo "  make test-cov          Run tests with coverage report"
	@echo "  make lint              Run linters (black, flake8, isort)"
	@echo "  make format            Auto-format code"
	@echo "  make type              Run type checker (mypy)"
	@echo "  make quality           Run all quality checks"
	@echo ""
	@echo "Build & Distribution:"
	@echo "  make build             Build wheel and sdist packages"
	@echo "  make docs              Build documentation"
	@echo "  make clean             Clean build artifacts"
	@echo ""

install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

clean:
	rm -rf build/ dist/ *.egg-info/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
	find . -type d -name ".mypy_cache" -exec rm -rf {} +
	find . -type d -name ".tox" -exec rm -rf {} +

test:
	pytest tests/

test-cov:
	pytest --cov=src --cov-report=html --cov-report=term-missing tests/
	@echo "Coverage report generated in htmlcov/index.html"

lint:
	@echo "Running black..."
	black --check src/ tests/
	@echo "Running isort..."
	isort --check-only src/ tests/
	@echo "Running flake8..."
	flake8 src/ tests/

format:
	@echo "Formatting with black..."
	black src/ tests/
	@echo "Sorting imports with isort..."
	isort src/ tests/
	@echo "Done!"

type:
	mypy src/

quality: lint type

build: clean
	@echo "Building distribution packages..."
	python -m build
	@echo "Build complete! Packages in dist/"
	@ls -lh dist/

docs:
	@echo "Building documentation..."
	sphinx-build -W -b html docs docs/_build/html
	@echo "Documentation built in docs/_build/html/"

.DEFAULT_GOAL := help
