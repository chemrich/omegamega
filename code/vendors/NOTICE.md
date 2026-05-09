# Vendored from construct_compiler

Files in this directory are copied verbatim from
[chemrich/construct_compiler](https://github.com/chemrich/construct-compiler)
at commit `9f42e8b01ae7f30680cf8f71eb1300c23327af82` (2026-05-09).

| File | Source path |
|---|---|
| `base.py` | `src/construct_compiler/vendors/base.py` |
| `twist.py` | `src/construct_compiler/vendors/twist.py` |

These are kept as a vendored copy rather than a Python dependency to keep
omegamega self-contained. If upstream changes are needed, re-copy the files
manually and update the commit SHA above.

`twist.py` is used by `code/pricing.py` for the optional `--twist-quote` live
pricing path. Authentication uses the env vars documented at the top of
`twist.py`.
