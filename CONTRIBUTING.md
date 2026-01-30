# Contributing to FragMentor

Thank you for your interest in contributing to FragMentor! 🧬

## Getting Started

### Development Setup

```bash
# Clone the repository
git clone https://github.com/fragmentomics/fragmentomics.git
cd fragmentomics

# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=fragmentomics --cov-report=html

# Run specific test file
pytest tests/test_sizes.py -v
```

### Code Style

We use:
- **Black** for code formatting
- **Ruff** for linting
- **MyPy** for type checking

```bash
# Format code
black src/ tests/

# Lint
ruff check src/

# Type check
mypy src/fragmentomics
```

## How to Contribute

### Reporting Bugs

1. Check existing issues first
2. Include Python version and OS
3. Provide minimal reproducible example
4. Include full error traceback

### Suggesting Features

1. Open an issue with `[Feature]` prefix
2. Describe the use case
3. Provide examples of expected behavior

### Pull Requests

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass
6. Update documentation if needed
7. Submit PR with clear description

### Commit Messages

We follow conventional commits:

```
feat: add new feature
fix: bug fix
docs: documentation changes
test: add or modify tests
refactor: code refactoring
chore: maintenance tasks
```

## Code Guidelines

### Python Style

- Use type hints for all public functions
- Write docstrings in NumPy format
- Keep functions focused and small
- Prefer composition over inheritance

### Documentation

- All public APIs must have docstrings
- Include examples in docstrings
- Update README for significant changes

### Testing

- Write tests for all new functionality
- Maintain or improve coverage
- Use descriptive test names
- Use fixtures for common test data

## Project Structure

```
fragmentomics/
├── src/fragmentomics/
│   ├── io/           # I/O utilities (BAM reading, etc.)
│   ├── features/     # Feature extraction modules
│   ├── viz/          # Visualization
│   ├── ml/           # Machine learning (future)
│   └── utils/        # Utilities
├── tests/            # Test suite
├── docs/             # Documentation
└── notebooks/        # Tutorial notebooks
```

## Questions?

- Open a GitHub issue
- Tag with `[Question]`

---

**Thank you for helping make FragMentor better!** 🙏
