# Security Policy

## Supported Versions

jwst-l2-coupled-dynamics is in early-release phase (v0.1.x). Only the latest
release tag on the `main` branch receives security attention.

| Version       | Supported |
|---------------|-----------|
| 0.1.x (HEAD)  | yes       |
| earlier       | no        |

## Scope

jwst-l2-coupled-dynamics is a numerical reference implementation of
coupled rigid-body dynamics for two extended bodies under mutual Newtonian
gravity. It is **not intended for production attitude-control or mission
planning** in its current form. It produces data files (NDJSON, JSON) and
plots; it does not accept network input, does not authenticate, does not
write outside the working directory.

Known intentional limitations (not treated as vulnerabilities):
- No input sanitisation on parameter dictionaries passed to integrators
  (callers are assumed to be the author or a contributor running locally)
- No bounds on integration time or step count (caller responsibility)
- Output files (NDJSON trajectories, JSON fixtures) may be large; no
  size guard

## Reporting a Vulnerability

If you discover a security issue that would affect a user running this code
in a networked or shared environment, please report it privately before
opening a public issue.

**Contact:** petr.yamyang@gmail.com
**Subject line:** `[jwst-l2-coupled-dynamics SECURITY] <brief description>`

Expected response time: within 7 days. If the issue is confirmed, a fix or
mitigation note will be published in the next commit and credited to the
reporter (unless anonymity is requested).

Do not open a public GitHub issue for security-sensitive findings until a
fix or documented mitigation is in place.
