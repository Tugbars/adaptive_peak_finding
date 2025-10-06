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

1. **Second-order-derivative core**
   Peaks are located at zero-crossings of the **second difference** (ΔΔ signal).
   This converts the search into a *sign-change* problem, making the detector largely insensitive to slow baseline drift and linear trends while keeping computational cost **O(N)**.

2. **Adaptive thresholding**
   Every candidate zero-crossing is validated against a **local** amplitude gate:
   `threshold = μ_local + k · σ_local`
   where μ and σ are computed with **Welford’s** stable on-line algorithm inside a user-defined window. Sensitivity therefore **tracks** non-stationary noise and signal power in real time.

3. **Topological (MATLAB-compatible) prominence**
   Prominence is calculated with the **topographic rule** used by MATLAB’s `findpeaks`: descend left and right until a higher summit or the data edge is met; the prominence is the peak height minus the **higher** of the two minima encountered. No Gaussian priors, no fitting – just the internationally accepted definition.

4. **Hysteresis-based width estimation**
   FWHM and 10%-height widths are measured with **hysteresis**: after the signal drops below the target level it must rise back by a configurable percentage (default **20%**) before the search stops. This suppresses “hair-line” widths caused by single outlier samples and gives repeatable results on noisy data.

5. **Flat-top plateau handling**
   When enabled, the divide-and-conquer peak search continues until the **entire plateau** is traversed; the **centre index** is returned. Useful for rectangular pulses, clipped peaks, or over-sampled data.

6. **Safe numerical handling**

   * Welford’s on-line algorithm for mean/variance – no catastrophic cancellation even when the window is large and the signal almost constant.
   * Mid-point calculation uses `l + (r - l) / 2` to avoid 32-bit overflow.
   * All floating work is done in `double`; user data may be `float`.

7. **Robust error handling**
   Functions return the `PeakResult` enum (`PEAK_OK`, `PEAK_NO_MEMORY`, `PEAK_INVALID_CONFIG`, …). Callers can switch on the exact cause instead of guessing from a magic `-1`.

8. **Memory safety & boundary guards**

   * Dynamic scratch buffers are **always** bounded by `peak_detection_window_size`; `malloc` failures are caught and reported.
   * Window indices are clipped to `[0 … length-1]`; the code never accesses out-of-range samples even when the peak sits at the very edge of the buffer.

9. **Zero-dependency, single-file ANSI-C**
   No third-party libraries, no C++, no dynamic memory on the hot path (except the optional scratch buffers), making the library suitable for bare-metal MCUs, DSPs, and FPGA soft-cores.

10. **Thread-aware design**

    * No global mutable state (after making `default_config` `const`).
    * Optional user-supplied log callback can be added without breaking the ABI.


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
