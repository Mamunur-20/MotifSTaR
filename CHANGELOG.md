# Changelog

All notable changes to MotifSTaR will be recorded here.

## 1.4.7 — 2026

- Public release branding changed from the development name **STR-Interruption** to **MotifSTaR**.
- Compound-panel equality is determined from the complete ordered vector of displayed repeat-block sizes, including companion blocks that are not part of a primary-only report.
- Panel-specific calls follow SVG vertical order: the top panel is allele 1 and the bottom panel is allele 2.
- Supported same-size calls are retained without adding a review solely because their top/bottom order cannot be established from pooled evidence.
- Primary-repeat-only behaviour remains available through `--component-mode primary`; this is the recommended mode for the biological target in the supplied catalogue.

The branding update does not alter the calling algorithm relative to the supplied 1.4.7 development script.
