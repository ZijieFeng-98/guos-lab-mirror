# Contributing to qPCR Cleaner and Analyzer

Thank you for your interest in contributing to the qPCR Cleaner and Analyzer! This document provides guidelines for contributing to the project.

## 🚀 How to Contribute

### Reporting Bugs
- Use the [bug report template](.github/ISSUE_TEMPLATE/bug_report.md)
- Include detailed steps to reproduce the issue
- Provide sample data if possible
- Specify your environment (OS, R version, browser)

### Suggesting Features
- Use the [feature request template](.github/ISSUE_TEMPLATE/feature_request.md)
- Describe the use case and benefits
- Consider implementation complexity
- Check if similar features already exist

### Code Contributions
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Test thoroughly with sample data
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

## 🛠️ Development Setup

### Prerequisites
- R 4.0.0 or higher
- RStudio (recommended)
- Git

### Local Development
1. Clone the repository:
   ```bash
   git clone https://github.com/ZijieFeng-98/qPCR_Cleaner_Analyzer.git
   cd qPCR_Cleaner_Analyzer
   ```

2. Install required packages:
   ```r
   source("install_packages.R")
   ```

3. Run the app locally:
   ```r
   shiny::runApp()
   ```

4. Make your changes and test thoroughly

## 📋 Code Style Guidelines

### R Code
- Use meaningful variable and function names
- Add comments for complex logic
- Follow R style guidelines
- Use consistent indentation (2 spaces)

### Shiny App Structure
- Keep UI and server logic organized
- Use reactive expressions for complex calculations
- Handle errors gracefully
- Provide user feedback

### Documentation
- Update README.md for new features
- Add inline comments for complex functions
- Update sample data if needed

## 🧪 Testing

### Before Submitting
- Test with the included sample data
- Test with different file formats (CSV, Excel)
- Test with various data sizes
- Test error conditions
- Test on different browsers

### Test Cases
- Valid data with 3 replicates
- Valid data with 4 replicates
- Data with missing values
- Data with outliers
- Invalid file formats
- Empty files

## 📝 Pull Request Guidelines

### Before Submitting a PR
- [ ] Code follows style guidelines
- [ ] All tests pass
- [ ] Documentation is updated
- [ ] No sensitive information is included
- [ ] Changes are focused and minimal

### PR Description
- Describe the changes made
- Explain the rationale
- Include screenshots if UI changes
- Reference related issues

## 🏷️ Issue Labels

- `bug` - Something isn't working
- `enhancement` - New feature or request
- `documentation` - Improvements to documentation
- `good first issue` - Good for newcomers
- `help wanted` - Extra attention is needed

## 📞 Getting Help

- Check existing issues and discussions
- Join our community discussions
- Contact the maintainers directly

## 📄 License

By contributing, you agree that your contributions will be licensed under the same license as the project.

---

Thank you for contributing to the qPCR Cleaner and Analyzer! 🎉
