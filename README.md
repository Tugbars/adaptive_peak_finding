# Adaptive Peak Finding and Analysis

**Author:** Tugbars Heptaskin
**Date:** 10/06/2025
**Company:** Aminic Aps

---

## Overview

This module provides an **adaptive peak detection and analysis algorithm** for noisy signals.
It dynamically adjusts its detection criteria based on **local signal statistics** such as mean and standard deviation, allowing it to perform reliably across varying noise levels and signal amplitudes.

Unlike simple threshold-based peak detection, this implementation incorporates **adaptive thresholds**, **hysteresis**, and **topological prominence** (MATLAB-compatible) to do peak identification.
It also computes **width metrics** (Full Width at Half Maximum – FWHM, and width at 10% height) to characterize peak shapes.

---

## Features

* **Adaptive thresholding:** Uses local mean and standard deviation to dynamically adjust detection sensitivity.
* **Topological prominence:** Compatible with MATLAB’s `findpeaks` prominence calculation.
* **Hysteresis-based width estimation:** Improves robustness in noisy data.
* **Flat-top handling:** Optionally detects plateaus and returns their center.
* **Safe numerical handling:**

  * Welford’s algorithm for local statistics (avoids floating-point instability).
  * Safe midpoint computation in recursion to prevent integer overflow.
* **Error handling:** Uses a return code enum (`PeakResult`) for precise error reporting.
* **Memory safety:**

  * Dynamic allocation with failure handling.
  * Guard checks for array boundaries and invalid configurations.

---

---

## Key Functions

### `PeakResult processPeak(MqsRawDataPoint_t a[], int size, uint16_t* peakIndex, bool* isEdgeCase, const PeakFinderConfig* user_config)`

**Description:**
Main entry point for the peak finding process. Finds and verifies the **widest valid peak** based on adaptive criteria.

**Parameters:**

* `a[]` — Input signal array (`MqsRawDataPoint_t`).
* `size` — Length of the signal array.
* `peakIndex` — Output index of the detected peak.
* `isEdgeCase` — Set `true` if the peak lies near the boundary and still climbing.
* `user_config` — Optional configuration pointer (use `NULL` for defaults).

**Returns:**
A `PeakResult` code indicating success or failure.

---

## Debugging

If `DEBUG_PRINT` is defined, verbose messages will be printed for:

* Local statistics and adaptive thresholds
* Peak indices, values, and prominences
* Width measurements and climbing detection

> **Note:** `DEBUG_PRINT` is **not thread-safe**. For multi-threaded use, replace `printf` with a callback-based logging mechanism.

---

## Error Handling

The algorithm gracefully handles:

* Invalid configurations (negative or zero window sizes)
* Memory allocation failures
* Out-of-bounds access during recursive or adaptive detection
* Missing or insignificant peaks

Return codes (`PeakResult`) clearly indicate the type of failure.

---

## Performance Notes

* Complexity: O(N) per window range (gradient-based detection).
* Memory usage: Minimal dynamic allocation for temporary buffers.
* Suitable for **embedded or real-time environments** with moderate sampling rates.
* Designed for **deterministic operation** — no randomization or adaptive learning loops.

---

## Future Improvements

* Add **configurable logging callback** for multi-threaded environments.
* Introduce **peak merging** for overlapping peaks.
* Add **moving-window continuous mode** for real-time signal streams.
* Optimize for **fixed-point arithmetic** in low-resource MCUs.

---

## License

Proprietary – © 2025 Aminic Aps.
All rights reserved.

---

Would you like me to include a short **“Signal flow diagram”** (showing how each step leads to the next: local stats → adaptive threshold → gradient → prominence → width → validation)? It’d make the README more visual and helpful for onboarding new engineers.
