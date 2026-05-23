#pragma once
//
// Audit the printed form of registered models for suspicious tokens.
//
// Diagnostic complement to testmodelequations: rather than comparing
// against a locked canonical string, this test scans the printed
// expression for substrings that should NEVER appear given the model's
// source (e.g. trig functions of spatial coordinates in a model that
// never invokes trig).  Useful for catching cases where a macro like
// e(x) silently means something different from what the model author
// intended.
//
// Build:  -DCHOSEN_TESTS=testequationaudit
//
// Returns 0 on success, otherwise nonzero count of failing audits.

int testequationaudit();
