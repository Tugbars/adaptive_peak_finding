/*!
 * Adaptive Peak Finding and Analysis Implementation.
 *
 * Author: Tugbars Heptaskin
 * Date: 06/18/2024
 * Company: Aminic Aps
 *
 * This implementation provides algorithms for detecting and analyzing peaks in noisy signals.
 * The approach is based on adaptive peak detection, which dynamically adjusts detection thresholds
 * and criteria based on local signal properties, such as mean and standard deviation.
 *
 * The algorithms utilize local mean and standard deviation to find peaks that adhere to the overall
 * trend of the dataset in a given window. The implementation considers both the full width at half
 * maximum (FWHM) and the width at 10% of the peak height, ensuring robust peak detection and analysis.
 *
 * Improvements Addressed:
 * - Replaced magic numbers with a PeakFinderConfig struct for configurability.
 * - Fixed boundary bugs: Guarded against window_size > length in calculate_local_stats; used dynamic
 *   allocation in adaptive_gradient_find_peaks.
 * - Added hysteresis to width estimation in calculate_peak_widths for robustness to noise.
 * - Added optional flat-top handling in findPeakRec (enabled via config).
 * - Used safe midpoint calculation in findPeakRec to prevent overflow.
 * - Improved numerical precision in calculate_local_stats using Welford's algorithm.
 * - Added return code enum (PeakResult) instead of bare bool for better error handling.
 * - Replaced original prominence calculation with MATLAB-compatible topological prominence,
 *   which uses the higher of the left and right minima as the reference level.
 * - DEBUG_PRINT retained but noted as not thread-safe; consider callback for multi-threading.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include "adaptive_peak_finding.h" // Defines MqsRawDataPoint_t and necessary types

// Define return code enum for error handling
typedef enum {
    PEAK_OK,              // Successful operation
    PEAK_NO_MEMORY,       // Memory allocation failure
    PEAK_INVALID_CONFIG,  // Invalid configuration parameters
    PEAK_NO_PEAK_FOUND,   // No valid peak found
    PEAK_RANGE_TOO_LARGE  // Range exceeds allowed size
} PeakResult;

// Configuration struct to replace hard-coded constants
typedef struct {
    float adaptive_threshold_multiplier; // Multiplier for adaptive threshold (e.g., 0.7)
    float low_height_threshold;          // Threshold for low width (e.g., 0.1 for 10%)
    float half_height_threshold;         // Threshold for FWHM (e.g., 0.5 for 50%)
    float hysteresis_factor;             // Factor for hysteresis (e.g., 1.2 for 20% rise)
    float prominence_threshold;          // Minimum prominence for valid peaks
    float noise_tolerance;               // Noise tolerance for climbing check
    int window_size;                     // Window size for local statistics
    int peak_detection_window_size;      // Max range for peak detection
    int peak_threshold;                  // Threshold for edge case detection
    bool flat_top_mode;                  // Enable flat-top handling (return plateau center)
    // void (*log_callback)(const char* msg); // Optional logging callback (not implemented)
} PeakFinderConfig;

// Default configuration values (overridable by user)
static PeakFinderConfig default_config = {
    .adaptive_threshold_multiplier = 0.7f,
    .low_height_threshold = 0.1f,
    .half_height_threshold = 0.5f,
    .hysteresis_factor = 1.2f,
    .prominence_threshold = 1.0f,    // Default prominence threshold
    .noise_tolerance = 0.01f,        // Default noise tolerance
    .window_size = 20,               // Default window size for local stats
    .peak_detection_window_size = 100, // Default max range for peak detection
    .peak_threshold = 5,             // Default edge threshold
    .flat_top_mode = true            // Enable flat-top handling by default
};

/*!
 * @brief Determines if a peak is still climbing at the end of a dataset.
 *
 * This function assesses whether the identified peak in a dataset is still rising
 * as it reaches the end of the dataset. This is important in peak finding algorithms,
 * particularly when analyzing segments of data where the peak might extend beyond
 * the current dataset's boundary.
 *
 * The function iterates from the peak index to the end of the dataset, calculating
 * the derivative (rate of change) at each point. It checks if this derivative is
 * less than or equal to a specified noise tolerance. If the condition fails more than once,
 * it indicates that the peak is no longer climbing.
 *
 * @param b The array of data points (MqsRawDataPoint_t) containing the peak.
 * @param sizeB The size of the array.
 * @param peakIndex The index of the peak within the array.
 * @param config The configuration struct containing noise_tolerance.
 * @return True if the peak is still climbing; false otherwise.
 */
static bool isPeakClimbing(const MqsRawDataPoint_t b[], int sizeB, int peakIndex, const PeakFinderConfig* config) {
    // Validate peak index to ensure we can compute derivatives
    if (peakIndex <= 0 || peakIndex >= sizeB - 1) {
        return false; // Cannot climb if at or near array bounds
    }
    int failCount = 0; // Track consecutive points where derivative <= noise tolerance
    for (int i = peakIndex; i < sizeB - 1; i++) {
        // Compute derivative as difference between consecutive points
        float derivativeAfter = b[i + 1].phaseAngle - b[i].phaseAngle;
        if (derivativeAfter <= config->noise_tolerance) {
            failCount++;
            if (failCount >= 2) {
                return false; // Peak stops climbing after two small derivatives
            }
        }
    }
    return failCount < 2; // Peak is climbing if fewer than two failures
}

/*!
 * @brief Finds the index of the maximum value in a column of a 1D array.
 *
 * @param a The array of data points (MqsRawDataPoint_t) to search through.
 * @param l The starting index.
 * @param r The ending index.
 * @param max_val A pointer to store the maximum value found.
 * @param max_index A pointer to store the index of the maximum value.
 * @return The index of the maximum value found in the specified column.
 */
static inline int maxrow(const MqsRawDataPoint_t a[], int l, int r, float* max_val, int* max_index) {
    // Initialize max value to negative infinity
    *max_val = -INFINITY;
    // Iterate through the specified range to find the maximum
    for (int i = l; i <= r; i++) {
        if (*max_val < a[i].phaseAngle) {
            *max_val = a[i].phaseAngle;
            *max_index = i;
        }
    }
    return *max_index;
}

/*!
 * @brief Recursively finds a peak in a dataset using a divide-and-conquer approach.
 *
 * This function implements a recursive peak finding algorithm. It divides the dataset
 * into two halves at each recursive step and determines the direction (left or right)
 * to continue the search based on the comparison of adjacent values. If flat_top_mode
 * is enabled, it detects plateaus and returns the center index.
 *
 * @param a The array of data points (MqsRawDataPoint_t) to search through for a peak.
 * @param l The starting index of the current search window.
 * @param r The ending index of the current search window.
 * @param peakIndex A pointer to store the index of the found peak.
 * @param config The configuration struct (for flat_top_mode).
 * @return The value of the peak found, or -1 if no peak is found.
 */
static double findPeakRec(const MqsRawDataPoint_t a[], int l, int r, uint16_t* peakIndex, const PeakFinderConfig* config) {
    // Base case: invalid range
    if (l > r) {
        return -1;
    }
    // Compute midpoint safely to avoid integer overflow
    int mid = l + (r - l) / 2;
    float mid_value = a[mid].phaseAngle;
    // Handle edge cases: if mid is at the boundary, return it
    if (mid == l || mid == r) {
        *peakIndex = mid;
        return mid_value;
    }
    // Compare with neighbors to decide search direction
    if (a[mid - 1].phaseAngle > mid_value) {
        // Left neighbor is larger; search left half
        return findPeakRec(a, l, mid - 1, peakIndex, config);
    } else if (a[mid + 1].phaseAngle > mid_value) {
        // Right neighbor is larger; search right half
        return findPeakRec(a, mid + 1, r, peakIndex, config);
    } else {
        // Current point is a peak (or plateau)
        *peakIndex = mid;
        if (config->flat_top_mode) {
            // Detect plateau by expanding left and right to find equal values
            int left = mid, right = mid;
            while (left > l && a[left - 1].phaseAngle == mid_value) left--;
            while (right < r && a[right + 1].phaseAngle == mid_value) right++;
            *peakIndex = (left + right) / 2; // Return center of plateau
        }
        return mid_value;
    }
}

/*!
 * @brief Calculate local mean and standard deviation for a given index in the signal.
 *
 * This function computes the local mean and standard deviation of the signal within a specified
 * window size centered at a given index. Uses Welford's online algorithm for numerical stability.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param config The configuration struct containing window_size.
 * @param index The index at which to compute the local statistics.
 * @param local_mean Pointer to store the computed local mean.
 * @param local_std Pointer to store the computed local standard deviation.
 * @return PEAK_OK on success, PEAK_INVALID_CONFIG if window_size invalid.
 */
static PeakResult calculate_local_stats(const MqsRawDataPoint_t* signal, int length, const PeakFinderConfig* config, int index, double* local_mean, double* local_std) {
    // Clamp window size to prevent invalid access
    int window_size = config->window_size;
    if (window_size <= 0 || window_size > length) {
        window_size = length; // Use entire signal if window_size is invalid
    }
    // Compute window boundaries, ensuring they stay within signal bounds
    int start = fmax(0, index - window_size / 2);
    int end = fmin(length - 1, index + window_size / 2);
    // Initialize variables for Welford's algorithm
    double mean = 0.0, M2 = 0.0;
    int count = 0;
    // Single-pass computation of mean and variance
    for (int i = start; i <= end; ++i) {
        double value = signal[i].phaseAngle;
        count++;
        double delta = value - mean;
        mean += delta / count;
        double delta2 = value - mean;
        M2 += delta * delta2; // Accumulate sum of squared differences
    }
    // Compute final mean and sample standard deviation
    *local_mean = mean;
    *local_std = (count > 1) ? sqrt(M2 / (count - 1)) : 0.0; // Use count-1 for sample variance
    return PEAK_OK;
}

/*!
 * @brief Detects peaks in the signal using an adaptive gradient method.
 *
 * This function identifies peaks in the given signal by analyzing the gradient within a specified range
 * and using a sliding window to compute local statistics. Uses dynamic allocation for scratch arrays.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param start The starting index for peak detection.
 * @param end The ending index for peak detection.
 * @param config The configuration struct.
 * @param peaks Array to store the detected peak indices (caller must allocate).
 * @param num_peaks Pointer to store the number of detected peaks.
 * @return PEAK_OK on success, PEAK_NO_MEMORY on allocation failure.
 */
static PeakResult adaptive_gradient_find_peaks(const MqsRawDataPoint_t* signal, int length, int start, int end, const PeakFinderConfig* config, int* peaks, int* num_peaks) {
    // Validate range size
    int range_length = end - start + 1;
    if (range_length <= 3) {
        *num_peaks = 0; // Too small for derivative calculations
        return PEAK_OK;
    }
    // Allocate dynamic arrays for gradient calculations
    double* dY = (double*)malloc((range_length - 2) * sizeof(double));
    int* S = (int*)malloc((range_length - 2) * sizeof(int));
    int* ddS = (int*)malloc((range_length - 4) * sizeof(int));
    if (!dY || !S || !ddS) {
        free(dY);
        free(S);
        free(ddS);
        return PEAK_NO_MEMORY; // Handle allocation failure
    }
    // Compute first derivative (central difference)
    double local_mean, local_std;
    for (int i = start + 1; i < end; ++i) {
        dY[i - start - 1] = signal[i + 1].phaseAngle - signal[i - 1].phaseAngle;
    }
    // Compute sign of first derivative (-1, 0, or 1)
    for (int i = 0; i < range_length - 2; ++i) {
        S[i] = (dY[i] > 0) - (dY[i] < 0);
    }
    // Compute second derivative as difference of signs
    for (int i = 1; i < range_length - 3; ++i) {
        ddS[i - 1] = S[i] - S[i - 1];
    }
    // Detect peaks where second derivative indicates a peak (ddS = -2)
    *num_peaks = 0;
    for (int i = 1; i < range_length - 3; ++i) {
        if (ddS[i - 1] == -2) {
            int peak_index = start + i + 1; // Adjust index for convolution offset
            // Verify peak using local statistics
            calculate_local_stats(signal, length, config, peak_index, &local_mean, &local_std);
            double adaptive_threshold = local_mean + config->adaptive_threshold_multiplier * local_std;
            if (signal[peak_index].phaseAngle > adaptive_threshold) {
                peaks[*num_peaks] = peak_index;
                (*num_peaks)++;
            }
        }
    }
    // Clean up allocated memory
    free(dY);
    free(S);
    free(ddS);
    return PEAK_OK;
}

/*!
 * @brief Calculate the width of each detected peak.
 *
 * This function computes the width of each detected peak in the signal based on both the full width
 * at half maximum (FWHM) and the width at low_height_threshold (e.g., 10%) of the peak height.
 * Includes hysteresis to handle noise.
 *
 * @param signal The input signal array.
 * @param peaks Array of detected peak indices.
 * @param num_peaks The number of detected peaks.
 * @param widths Array to store the computed widths of the peaks.
 * @param length The length of the input signal array.
 * @param config The configuration struct.
 */
static void calculate_peak_widths(const MqsRawDataPoint_t* signal, const int* peaks, int num_peaks, double* widths, int length, const PeakFinderConfig* config) {
    for (int i = 0; i < num_peaks; i++) {
        int peak = peaks[i];
        double peak_height = signal[peak].phaseAngle;
        // Compute thresholds for low width (e.g., 10%) and FWHM (e.g., 50%)
        double low_threshold = peak_height * config->low_height_threshold;
        double half_threshold = peak_height * config->half_height_threshold;
        double low_hyst = low_threshold * config->hysteresis_factor; // Hysteresis for low threshold
        double half_hyst = half_threshold * config->hysteresis_factor; // Hysteresis for FWHM
        // Left side for low threshold
        int start = peak;
        bool crossed_low = false;
        while (start > 0) {
            if (signal[start].phaseAngle <= low_threshold) {
                crossed_low = true; // Signal dropped below threshold
            }
            if (crossed_low && signal[start].phaseAngle > low_hyst) {
                break; // Stop if signal rises above hysteresis level
            }
            start--;
        }
        // Right side for low threshold
        int stop = peak;
        crossed_low = false;
        while (stop < length) {
            if (signal[stop].phaseAngle <= low_threshold) {
                crossed_low = true;
            }
            if (crossed_low && signal[stop].phaseAngle > low_hyst) {
                break;
            }
            stop++;
        }
        // Left side for half threshold (FWHM)
        int left_half = peak;
        bool crossed_half = false;
        while (left_half > 0) {
            if (signal[left_half].phaseAngle <= half_threshold) {
                crossed_half = true;
            }
            if (crossed_half && signal[left_half].phaseAngle > half_hyst) {
                break;
            }
            left_half--;
        }
        // Right side for half threshold (FWHM)
        int right_half = peak;
        crossed_half = false;
        while (right_half < length) {
            if (signal[right_half].phaseAngle <= half_threshold) {
                crossed_half = true;
            }
            if (crossed_half && signal[right_half].phaseAngle > half_hyst) {
                break;
            }
            right_half++;
        }
        // Compute weighted average of FWHM and low width
        double fwhm_width = right_half - left_half;
        double low_width = stop - start;
        widths[i] = fwhm_width * config->half_height_threshold + low_width * (1.0 - config->half_height_threshold);
    }
}

/*!
 * @brief Calculate the topological prominence of a given peak (MATLAB-compatible).
 *
 * This function computes the prominence of a peak as the difference between the peak value
 * and the higher of the minima on either side, following MATLAB's findpeaks approach.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param peakIndex The index of the peak for which to calculate prominence.
 * @param refLevelOut Optional pointer to store the reference level (higher minimum).
 * @return The computed prominence of the peak.
 */
static double topological_prominence(const MqsRawDataPoint_t* signal, int length, int peakIndex, double* refLevelOut) {
    double pv = signal[peakIndex].phaseAngle; // Peak value
    // Left contour: find minimum or higher peak
    int leftRef = 0;
    double leftMin = pv;
    for (int i = peakIndex - 1; i >= 0; --i) {
        if (signal[i].phaseAngle >= pv) {
            leftRef = i; // Stop at higher or equal peak
            break;
        }
        if (signal[i].phaseAngle < leftMin) {
            leftMin = signal[i].phaseAngle; // Update minimum
        }
    }
    // Right contour: find minimum or higher peak
    int rightRef = length - 1;
    double rightMin = pv;
    for (int i = peakIndex + 1; i < length; ++i) {
        if (signal[i].phaseAngle >= pv) {
            rightRef = i;
            break;
        }
        if (signal[i].phaseAngle < rightMin) {
            rightMin = signal[i].phaseAngle;
        }
    }
    // Use higher of the two minima as reference level
    double refLevel = (leftMin > rightMin) ? leftMin : rightMin;
    if (refLevelOut) *refLevelOut = refLevel;
#ifdef DEBUG_PRINT
    // Note: DEBUG_PRINT is not thread-safe; consider logging callback for multi-threaded use
    printf("Topological prominence for peak at %d: %f\n", peakIndex, pv - refLevel);
#endif
    return pv - refLevel;
}

/*!
 * @brief Wrapper for topological prominence with original prototype.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param peakIndex The index of the peak for which to calculate prominence.
 * @return The computed prominence of the peak.
 */
static double find_prominence(const MqsRawDataPoint_t* signal, int length, int peakIndex) {
    // Delegate to topological_prominence without reference level output
    return topological_prominence(signal, length, peakIndex, NULL);
}

/*!
 * @brief Calculate the prominences of all detected peaks.
 *
 * This function computes the topological prominence of each detected peak in the signal.
 *
 * @param signal The input signal array.
 * @param peaks Array of detected peak indices.
 * @param num_peaks The number of detected peaks.
 * @param prominences Array to store the computed prominences of the peaks.
 * @param length The length of the input signal array.
 */
static void calculate_peak_prominences(const MqsRawDataPoint_t* signal, const int* peaks, int num_peaks, double* prominences, int length) {
    for (int i = 0; i < num_peaks; ++i) {
        // Compute prominence for each peak
        prominences[i] = topological_prominence(signal, length, peaks[i], NULL);
#ifdef DEBUG_PRINT
        // Note: DEBUG_PRINT is not thread-safe
        printf("Peak %d (idx %d) prominence %f\n", i + 1, peaks[i], prominences[i]);
#endif
    }
}

/*!
 * @brief Find the primary peak in the signal.
 *
 * This function identifies the primary peak in the signal by recursively searching for the maximum
 * value.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param config The configuration struct.
 * @return The index of the primary peak.
 */
static uint16_t find_primary_peak(const MqsRawDataPoint_t* signal, int length, const PeakFinderConfig* config) {
    uint16_t primary_peak = 0;
    // Use recursive peak finding to locate primary peak
    float peak_value = findPeakRec(signal, 0, length - 1, &primary_peak, config);
#ifdef DEBUG_PRINT
    // Note: DEBUG_PRINT is not thread-safe
    printf("Primary peak at index %d with value %f\n", primary_peak, peak_value);
#endif
    return primary_peak;
}

/*!
 * @brief Detect peaks within a specified range in the signal.
 *
 * This function identifies peaks within a specified range in the signal using an adaptive gradient method.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param range_start The starting index of the range for peak detection.
 * @param range_end The ending index of the range for peak detection.
 * @param config The configuration struct.
 * @param peaks Array to store the detected peak indices (caller allocates).
 * @return The number of detected peaks, or negative on error.
 */
static int detect_peaks_in_range(const MqsRawDataPoint_t* signal, int length, int range_start, int range_end, const PeakFinderConfig* config, int* peaks) {
    int num_peaks = 0;
    // Delegate to adaptive gradient method for peak detection
    PeakResult res = adaptive_gradient_find_peaks(signal, length, range_start, range_end, config, peaks, &num_peaks);
    if (res != PEAK_OK) {
#ifdef DEBUG_PRINT
        printf("Error in peak detection: %d\n", res);
#endif
        return -1; // Propagate error
    }
#ifdef DEBUG_PRINT
    // Note: DEBUG_PRINT is not thread-safe
    printf("Detected Peaks:\n");
    for (int i = 0; i < num_peaks; i++) {
        printf("Peak at index %d, value = %f\n", peaks[i], signal[peaks[i]].phaseAngle);
    }
#endif
    return num_peaks;
}

/*!
 * @brief Detects the peak with the widest width in the signal.
 *
 * This function identifies the peak with the widest width in the signal using the adaptive gradient method
 * and verifies its prominence.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param peaks Detected peaks.
 * @param num_peaks Number of peaks.
 * @param widths Array for widths (caller allocates).
 * @param widest_peak_index Pointer to store the index of the widest peak.
 * @param config The configuration struct.
 * @return PEAK_OK if found, PEAK_NO_PEAK_FOUND if no valid peak, PEAK_NO_MEMORY on alloc fail.
 */
static PeakResult find_widest_peak(const MqsRawDataPoint_t* signal, int length, int* peaks, int num_peaks, double* widths, uint16_t* widest_peak_index, const PeakFinderConfig* config) {
    if (num_peaks == 0) {
        return PEAK_NO_PEAK_FOUND; // No peaks to process
    }
    // Compute widths for all detected peaks
    calculate_peak_widths(signal, peaks, num_peaks, widths, length, config);
    // Allocate array for prominences
    double* prominences = (double*)malloc(num_peaks * sizeof(double));
    if (!prominences) {
        return PEAK_NO_MEMORY;
    }
    // Compute prominences for all peaks
    calculate_peak_prominences(signal, peaks, num_peaks, prominences, length);
#ifdef DEBUG_PRINT
    // Note: DEBUG_PRINT is not thread-safe
    printf("Peak Widths and Prominences:\n");
    for (int i = 0; i < num_peaks; i++) {
        printf("Peak at index %d has width %f and prominence %f\n", peaks[i], widths[i], prominences[i]);
    }
#endif
    // Find peak with maximum width that meets prominence threshold
    double max_width = 0.0;
    int max_width_index = -1;
    for (int i = 0; i < num_peaks; i++) {
        if (prominences[i] > config->prominence_threshold && widths[i] > max_width) {
            max_width = widths[i];
            max_width_index = i;
        }
    }
    free(prominences);
    if (max_width_index != -1) {
        *widest_peak_index = peaks[max_width_index];
        return PEAK_OK;
    }
    return PEAK_NO_PEAK_FOUND;
}

/*!
 * @brief Detects the peak with the widest width in the signal.
 *
 * This function identifies the peak with the widest width in the signal using the adaptive gradient method
 * and verifies its prominence.
 *
 * @param signal The input signal array.
 * @param length The length of the input signal array.
 * @param peak_range The range for detecting peaks.
 * @param widest_peak_index Pointer to store the index of the peak with the widest width.
 * @param config The configuration struct.
 * @return PEAK_OK if a valid peak is found; error code otherwise.
 */
static PeakResult detect_peak_with_width(const MqsRawDataPoint_t* signal, int length, int peak_range, uint16_t* widest_peak_index, const PeakFinderConfig* config) {
    // Find primary peak to center the search range
    uint16_t primary_peak = find_primary_peak(signal, length, config);
    // Compute search range around primary peak
    int range_start = primary_peak - peak_range / 2;
    int range_end = primary_peak + peak_range / 2;
    if (range_start < 0) range_start = 0; // Clamp to signal bounds
    if (range_end >= length) range_end = length - 1;
    // Allocate arrays for peaks and widths
    int max_possible_peaks = range_end - range_start + 1;
    int* peaks = (int*)malloc(max_possible_peaks * sizeof(int));
    double* widths = (double*)malloc(max_possible_peaks * sizeof(double));
    if (!peaks || !widths) {
        free(peaks);
        free(widths);
        return PEAK_NO_MEMORY;
    }
    // Detect peaks in the specified range
    int num_peaks = detect_peaks_in_range(signal, length, range_start, range_end, config, peaks);
    if (num_peaks < 0) {
        free(peaks);
        free(widths);
        return PEAK_NO_PEAK_FOUND; // Error in detection
    }
    // Find the widest peak
    PeakResult res = find_widest_peak(signal, length, peaks, num_peaks, widths, widest_peak_index, config);
    // Clean up allocated memory
    free(peaks);
    free(widths);
    return res;
}

/*!
 * @brief Process the input data to find and verify the peak with the widest width.
 *
 * This function processes the input data to identify and verify the peak with the widest width,
 * printing its index and magnitude if found.
 *
 * @param a The input signal array.
 * @param size The size of the input signal array.
 * @param peakIndex Pointer to store the index of the peak with the widest width.
 * @param isEdgeCase Pointer to store if it's an edge case (climbing).
 * @param user_config Optional user-provided config; uses default if NULL.
 * @return PEAK_OK if a valid peak with sufficient width and prominence is found; error code otherwise.
 */
PeakResult processPeak(MqsRawDataPoint_t a[], int size, uint16_t* peakIndex, bool* isEdgeCase, const PeakFinderConfig* user_config) {
    // Use user config if provided, otherwise default
    const PeakFinderConfig* config = user_config ? user_config : &default_config;
    // Validate configuration
    if (config->window_size <= 0 || config->peak_detection_window_size <= 0) {
        return PEAK_INVALID_CONFIG;
    }
    int peak_range = config->peak_detection_window_size;
    // Find the widest peak in the signal
    PeakResult res = detect_peak_with_width(a, size, peak_range, peakIndex, config);
    if (res != PEAK_OK) {
#ifdef DEBUG_PRINT
        // Note: DEBUG_PRINT is not thread-safe
        printf("No peak with width detected.\n");
#endif
        return res;
    }
    // Check if peak is near the end and still climbing
    if (*peakIndex >= size - config->peak_threshold) {
        *isEdgeCase = isPeakClimbing(a, size, *peakIndex, config);
    } else {
        *isEdgeCase = false;
    }
    // Report peak details
    double peak_magnitude = a[*peakIndex].phaseAngle;
#ifdef DEBUG_PRINT
    // Note: DEBUG_PRINT is not thread-safe
    printf("Peak found at index: %d\n", *peakIndex);
    printf("Peak magnitude: %f\n", peak_magnitude);
#endif
    return PEAK_OK;
}
