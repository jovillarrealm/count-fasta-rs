---
description: Check the rust project for errors, tests, and lints
---
# Rust Check Workflow

This workflow ensures the Rust project is correct by running `check`, `test`, and `clippy`. It stops if any step fails.

1. Run cargo check to verify compilation
// turbo
2. Run cargo test to verify functionality
// turbo
3. Run cargo clippy to catch lints
