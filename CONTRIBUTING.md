# Contributing to MolSwapper

Thank you for your interest in contributing to this project. This document provides guidelines for contributions.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/yourusername/MolSwapper.git
   cd MolSwapper
   ```
3. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```
4. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Development Workflow

1. Create a feature branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes with clear, well-commented code

3. Test thoroughly:
   - Run the CLI with various options and inputs
   - Verify output files (CSV, SDF, images) are correct
   - Test edge cases (single molecule, large libraries, etc.)

4. Commit with descriptive messages:
   ```bash
   git commit -m "Add feature: description of changes"
   ```

5. Push to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

6. Submit a Pull Request with:
   - Clear description of changes
   - Relevant issue references (e.g., "Fixes #123")
   - Explanation of the problem solved or feature added

## Code Standards

- Follow PEP 8 guidelines
- Use descriptive variable and function names
- Include docstrings for all functions and classes
- Keep functions focused and modular
- Ensure code is testable and maintainable

## Types of Contributions

### Bug Reports
- File reports using GitHub Issues
- Include: Python version, RDKit version, exact command, error message
- Provide a minimal reproducible example

### Feature Requests
- Clearly describe the proposed feature
- Explain the use case and benefits
- Suggest implementation approach if possible

### Documentation
- Fix typos and improve clarity
- Add examples for new features
- Enhance API documentation

### Code Improvements
- Optimize performance
- Refactor complex sections
- Add error handling and validation
- Improve modularity

## Questions

Open a GitHub issue or discussion for questions about contributing.

---

Thank you for contributing to MolSwapper.
