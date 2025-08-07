# Failing Tests Archive

This directory contains test files that were written during initial development but fail due to incomplete understanding of the Finesse 3 API.

## Files

- `test_finesse_interop.jl` - Complex Julia→Finesse tests that fail due to API issues
- `test_finesse_interop.py` - Complex Python/Finesse→Julia tests that fail due to API issues

## Purpose

These files are preserved to document:

1. **What doesn't work** with the current Finesse 3.0a33 API
2. **Why it fails** (detailed post-mortem analysis in file headers)
3. **Lessons learned** for future development
4. **API discovery insights** for when we better understand Finesse 3

## Key Lessons from These Failures

1. **Finesse 3 API differs significantly** from documentation and examples
2. **Component constructors require careful parameter handling** due to strict validation
3. **Model building is complex** and requires explicit port management
4. **Start with simple components** before attempting complex models
5. **Empirical API discovery** is needed rather than assumption-based development

## Status

These tests remain **FAILING** but contain valuable analysis. They should not be run as part of the main test suite.

## Future Work

When we gain better understanding of the Finesse 3 API through:
- Official documentation study
- Working example analysis
- Incremental component exploration

These tests can be revisited, corrected, and potentially moved back to the main test suite.