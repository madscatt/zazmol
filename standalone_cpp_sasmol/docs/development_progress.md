# Development Progress Notes

This file records design follow-ups that should not be lost while the
standalone C++ SasMol parity work proceeds in small implementation slices.

## Selection / Eval Scope

The first C++ selection parser intentionally covers only a bounded subset of
Python `get_subset_mask` expressions. A future pass should survey real Python
SasMol/SASSIE selection usage and determine the maximum practical scope that
can be supported safely.

That future expansion must keep the current failure contract:

- unsupported or ambiguous expressions fail with a clear error
- failures return no partial indices
- the parser must stop rather than guessing
- arbitrary Python-style `eval` behavior must not become the default C++ path
