# Product Guidelines: count-fasta-rs

## Prose & Tone
- **User-Centric Friendly:** Use language that is warm, helpful, and accessible. Avoid unnecessary jargon while maintaining technical accuracy.
- **Direct & Helpful:** Provide clear instructions and actionable feedback in error messages and documentation.

## User Experience (UX)
- **Dual Mode Output:**
    - **Human-Readable Mode:** Visually rich, formatted, and easy to read for quick manual checks in the terminal.
    - **CSV-Friendly Mode:** Standardized, machine-parsable output (via `--csv`) for seamless integration into pipelines and large-scale data analysis.
- **Interactive Feedback:** When running interactively, provide helpful progress indicators and clear summaries of results.

## Documentation Standards
- **Practicality Above All:** Focus on clear, copy-pasteable command examples that demonstrate real-world usage.
- **Scenario-Based Guides:** Provide guides for common tasks, such as processing large directories or integrating with other tools.

## Design Principles
- **Predictable Behavior:** Ensure the tool behaves consistently across different platforms and file formats.
- **Robust Error Handling:** Catch and report errors gracefully, providing hints for resolution where possible.
