package itu840;

import java.io.*;
import java.util.*;

/**
 * <p><b>ITU-R P.840-9 Cloud and Fog Attenuation Calculator</b></p>
 *
 * @author Mehmet Sait SEVER
 * @version 1.0.0, 26.09.2025
 *
 * <p><b>Mathematical Notation</b></p>
 * <ul>
 *   <li>L – Integrated cloud liquid water content (kg/m<sup>2</sup> or mm)</li>
 *   <li>p – Exceedance probability (CCDF) (percent)</li>
 *   <li>L(p) – Integrated cloud liquid water content at exceedance probability p (kg/m<sup>2</sup> or mm)</li>
 *   <li>m<sub>L</sub> and &sigma;<sub>L</sub> – Log-normal mean and standard deviation parameters (natural log scale)</li>
 *   <li>P<sub>L</sub> – Probability of cloud (percent)</li>
 * </ul>
 *
 * <p><b>Statistical Definitions</b></p>
 * <p><b>CDF</b> (cumulative distribution function) gives the probability that a random variable X is
 *    less than or equal to a value x: Pr[X &le; x].</p>
 * <p><b>CCDF</b> (complementary cumulative distribution function) gives the probability that X exceeds x:
 *    Pr[X &gt; x]. In this Recommendation, CCDF is used to express <b>exceedance probability</b>.</p>
 *
 * <p><b>Data Structures</b></p>
 * <p>Grids: Provided as digital maps with 721 × 1441 points at 0.25° resolution, starting from latitude −90° and longitude −180°.
 *    See the README file for details.</p>
 * <p>A <b>grid</b> is a global, regular latitude-longitude matrix that represents the Earth’s surface. Each cell (grid point)
 *    stores the value of a climatological parameter (e.g., L, m<sub>L</sub>, &sigma;<sub>L</sub>, or P<sub>L</sub>).</p>
 *
 * <p><b>This implementation follows Recommendation ITU-R P.840-9:</b></p>
 * <ul>
 *   <li>Equation (1): &gamma;<sub>c</sub>, the specific attenuation within the cloud or fog (Rayleigh approximation)</li>
 *   <li>Equation (2): K<sub>l</sub>, the cloud liquid water specific attenuation coefficient</li>
 *   <li>Equations (3)–(10): Double-Debye dielectric model for liquid water,
 *       including permittivity parameters (ε<sub>0</sub>, ε<sub>1</sub>, ε<sub>2</sub>),
 *       relaxation frequencies (f<sub>p</sub>, f<sub>s</sub>), and derived functions (ε&prime;, ε&Prime;, &eta;)</li>
 *   <li>Equation (11): A<sub>C</sub>, the slant path instantaneous cloud attenuation</li>
 *   <li>Equations (12), (14), (16): K<sub>L</sub>, the cloud liquid mass absorption coefficient</li>
 *   <li>Equation (13): A<sub>C</sub>, the slant path statistical cloud attenuation</li>
 *   <li>Equation (15): A<sub>C</sub>, the log-normal approximation to the slant path statistical cloud attenuation</li>
 *   <li>Section 4.2.1 and 4.2.2: Spatial and statistical interpolation rules</li>
 * </ul>
 */
public class Itu840 {

    /** <p>Utility class, no instances.</p> */
    private Itu840() {
        throw new AssertionError("Itu840 is a utility class and cannot be instantiated.");
    }

    // ===================== Input parameter ranges =====================
    // Latitude in degrees: [-90, +90]
    // Longitude in degrees: [-180, +180]
    // p in percent: [0.01, 100] for the annual statistics and [0.1, 100] for the monthly statistics
    // Frequency in GHz: [1–200]
    // Elevation angle (θ) in degrees: (0, 90]
    private static final double LATITUDE_MIN_DEGREES = -90.0;
    private static final double LATITUDE_MAX_DEGREES = 90.0;
    private static final double LONGITUDE_MIN_DEGREES = -180.0;
    private static final double LONGITUDE_MAX_DEGREES = 180.0;
    private static final double P_MIN_PERCENT_ANNUAL = 0.01;
    private static final double P_MIN_PERCENT_MONTHLY = 0.1;
    private static final double P_MAX_PERCENT = 100.0;
    private static final double FREQUENCY_MIN_GHZ = 1.0;
    private static final double FREQUENCY_MAX_GHZ = 200.0;
    private static final double ELEVATION_ANGLE_MIN_DEGREES = 0.0;
    private static final double ELEVATION_ANGLE_MAX_DEGREES = 90.0;

    // Reference temperature (273.75 K) used by ITU-R P.840-9
    private static final double REFERENCE_TEMPERATURE_KELVIN = 273.75;

    // Equation 15 NOTE threshold: P_L ≤ 0.02% ⇒ A_C = 0
    private static final double CLOUD_PROBABILITY_THRESHOLD_PERCENT = 0.02;

    // ===================== Gaussian parameters for K_L (Equations 12, 14, and 16) =====================
    private static final double GAUSSIAN_A1 = 0.1522;
    private static final double GAUSSIAN_A2 = 11.51;
    private static final double GAUSSIAN_A3 = -10.4912;
    private static final double GAUSSIAN_F1 = -23.9589;
    private static final double GAUSSIAN_F2 = 219.2096;
    private static final double GAUSSIAN_SIGMA1 = 3.2991e3;
    private static final double GAUSSIAN_SIGMA2 = 2.7595e6;

    // ===================== Double-Debye model parameters at reference temperature (273.75 K) =====================
    private static final double EPSILON_0 = computeEpsilon0();
    private static final double EPSILON_1 = computeEpsilon1();
    private static final double EPSILON_2 = computeEpsilon2();
    private static final double PRINCIPAL_RELAXATION_FREQUENCY = computePrincipalRelaxationFrequency();
    private static final double SECONDARY_RELAXATION_FREQUENCY = computeSecondaryRelaxationFrequency();

    // ===================== Grid specifications (see README Table 1) =====================
    private static final int NUMBER_OF_LATITUDE_POINTS = 721;
    private static final int NUMBER_OF_LONGITUDE_POINTS = 1441;
    private static final double GRID_LATITUDE_START_DEGREE = -90.0;
    private static final double GRID_LATITUDE_STEP_DEGREE = 0.25;
    private static final double GRID_LONGITUDE_START_DEGREE = -180.0;
    private static final double GRID_LONGITUDE_STEP_DEGREE = 0.25;

    // ==================================================================================
    //                                  Physical Formulas
    // ==================================================================================

    /**
     * <p>Computes ε<sub>0</sub>, a permittivity parameter at the fixed reference temperature of 273.75 K.</p>
     * <p>Defined in ITU-R P.840-9, Equation (6).</p>
     *
     * <p>At T = 273.75 K, ε<sub>0</sub> ≈ 87.5655 (dimensionless)</p>
     *
     * @return ε<sub>0</sub> (dimensionless)
     */
    public static double computeEpsilon0() {
        return 77.66 + 103.3 * (300.0 / REFERENCE_TEMPERATURE_KELVIN - 1.0);
    }

    /**
     * <p>Computes ε<sub>1</sub>, a permittivity parameter.</p>
     * <p>Defined in ITU-R P.840-9, Equation (7).</p>
     *
     * <p>At T = 273.75 K, ε<sub>1</sub> ≈ 5.8756 (dimensionless)</p>
     *
     * @return ε<sub>1</sub> (dimensionless)
     */
    public static double computeEpsilon1() {
        return 0.0671 * EPSILON_0;
    }

    /**
     * <p>Computes ε<sub>2</sub>, a constant permittivity parameter.</p>
     * <p>Defined in ITU-R P.840-9, Equation (8).</p>
     *
     * <p>At T = 273.75 K, ε<sub>2</sub> = 3.52 (dimensionless)</p>
     *
     * @return ε<sub>2</sub> (dimensionless)
     */
    public static double computeEpsilon2() {
        return 3.52;
    }

    /**
     * <p>Computes f<sub>p</sub>, the principal relaxation frequency at the fixed reference temperature of 273.75 K.</p>
     * <p>Defined in ITU-R P.840-9, Equation (9).</p>
     *
     * <p>At T = 273.75 K, f<sub>p</sub> ≈ 9.1056 GHz</p>
     *
     * @return f<sub>p</sub> in gigahertz
     */
    public static double computePrincipalRelaxationFrequency() {
        double x = 300.0 / REFERENCE_TEMPERATURE_KELVIN - 1.0;

        return 20.20 - 146.0 * x + 316.0 * x * x;
    }

    /**
     * <p>Computes f<sub>s</sub>, the secondary relaxation frequency.</p>
     * <p>Defined in ITU-R P.840-9, Equation (10).</p>
     *
     * <p>At T = 273.75 K, f<sub>s</sub> ≈ 362.4033 GHz</p>
     *
     * @return f<sub>s</sub> in gigahertz
     */
    public static double computeSecondaryRelaxationFrequency() {
        return 39.8 * PRINCIPAL_RELAXATION_FREQUENCY;
    }

    /**
     * <p>Computes ε&Prime;(f), the imaginary part of the complex permittivity of liquid water.</p>
     * <p>Defined in ITU-R P.840-9, Equation (4).</p>
     *
     * @param frequency Frequency in gigahertz
     * @return ε&Prime;(f) (dimensionless)
     */
    public static double computeEpsilonImaginary(double frequency) {
        double term1 = frequency * (EPSILON_0 - EPSILON_1) /
                       (PRINCIPAL_RELAXATION_FREQUENCY * (1 + Math.pow(frequency / PRINCIPAL_RELAXATION_FREQUENCY, 2)));
        double term2 = frequency * (EPSILON_1 - EPSILON_2) /
                       (SECONDARY_RELAXATION_FREQUENCY * (1 + Math.pow(frequency / SECONDARY_RELAXATION_FREQUENCY, 2)));

        return term1 + term2;
    }

    /**
     * <p>Computes ε&prime;(f), the real part of the complex permittivity of liquid water.</p>
     * <p>Defined in ITU-R P.840-9, Equation (5).</p>
     *
     * @param frequency Frequency in gigahertz
     * @return ε&prime;(f) (dimensionless)
     */
    public static double computeEpsilonReal(double frequency) {
        double term1 = (EPSILON_0 - EPSILON_1) / (1 + Math.pow(frequency / PRINCIPAL_RELAXATION_FREQUENCY, 2));
        double term2 = (EPSILON_1 - EPSILON_2) / (1 + Math.pow(frequency / SECONDARY_RELAXATION_FREQUENCY, 2));

        return term1 + term2 + EPSILON_2;
    }

    /**
     * <p>Computes &eta;(f), a parameter derived from the real and imaginary parts of the permittivity.</p>
     * <p>Defined in ITU-R P.840-9, Equation (3).</p>
     *
     * @param frequency Frequency in gigahertz
     * @return &eta;(f) (dimensionless)
     */
    public static double computeEta(double frequency) {
        double epsilonReal = computeEpsilonReal(frequency);
        double epsilonImaginary = computeEpsilonImaginary(frequency);

        return (2 + epsilonReal) / epsilonImaginary;
    }

    /**
     * <p>Computes K<sub>l</sub>(f), the cloud liquid water specific attenuation coefficient
     *    at the fixed reference temperature of 273.75 K.</p>
     * <p>Defined in ITU-R P.840-9, Equation (2).</p>
     *
     * @param frequency Frequency in gigahertz
     * @return K<sub>l</sub>(f) in (dB/km)/(g/m<sup>3</sup>)
     */
    public static double computeCloudLiquidWaterSpecificAttenuationCoefficient(double frequency) {
        double epsilonImaginary = computeEpsilonImaginary(frequency);
        double eta = computeEta(frequency);

        return 0.819 * frequency / (epsilonImaginary * (1.0 + eta * eta));
    }

    /**
     * <p>Computes K<sub>L</sub>(f), the cloud liquid mass absorption coefficient
     *    at the fixed reference temperature of 273.75 K.</p>
     * <p>Defined in ITU-R P.840-9, Equations (12, 14, 16).</p>
     *
     * @param frequency Frequency in gigahertz
     * @return K<sub>L</sub>(f) in dB/(kg/m<sup>2</sup>) (equivalently dB/mm)
     */
    public static double computeCloudLiquidMassAbsorptionCoefficient(double frequency) {
        double cloudLiquidWaterSpecificAttenuationCoefficient = computeCloudLiquidWaterSpecificAttenuationCoefficient(frequency);
        double term1 = GAUSSIAN_A1 * Math.exp(-Math.pow(frequency - GAUSSIAN_F1, 2) / GAUSSIAN_SIGMA1);
        double term2 = GAUSSIAN_A2 * Math.exp(-Math.pow(frequency - GAUSSIAN_F2, 2) / GAUSSIAN_SIGMA2);

        return cloudLiquidWaterSpecificAttenuationCoefficient * (term1 + term2 + GAUSSIAN_A3);
    }

    // ==================================================================================
    //                                Prediction Methods
    // ==================================================================================

    /**
     * <p>Slant path instantaneous cloud attenuation prediction method</p>
     * <p>Defined in ITU-R P.840-9, Equation (11).</p>
     *
     * <p>A<sub>C</sub>(f) = K<sub>L</sub>(f) · L / sin(&theta;)</p>
     *
     * @param frequency Frequency in gigahertz
     * @param integratedCloudLiquidWaterContent L in kg/m<sup>2</sup> (equivalently mm)
     * @param elevationAngle Elevation angle (&theta;) in degrees
     * @return Attenuation in dB
     */
    public static double computeSlantPathInstantaneousCloudAttenuation(
            double frequency,
            double integratedCloudLiquidWaterContent,
            double elevationAngle) {

        double cloudLiquidMassAbsorptionCoefficient = computeCloudLiquidMassAbsorptionCoefficient(frequency);
        double sineOfElevationAngle = Math.sin(Math.toRadians(elevationAngle));

        return cloudLiquidMassAbsorptionCoefficient * integratedCloudLiquidWaterContent / sineOfElevationAngle;
    }

    /**
     * <p>Slant path statistical cloud attenuation prediction method</p>
     * <p>Defined in ITU-R P.840-9, Equation (13).</p>
     *
     * <p>A<sub>C</sub>(f, p) = K<sub>L</sub>(f) · L(p) / sin(&theta;)</p>
     *
     * <p>where L(p) is obtained by:</p>
     * <ol>
     *   <li>Bilinear spatial interpolation at the two bracketing probabilities
     *       (p<sub>below</sub> and p<sub>above</sub>) (see Section 4.2.1a, 4.2.1b, and 4.2.1c)</li>
     *   <li>Linear interpolation of L with respect to log<sub>10</sub>(p) (see Section 4.2.1d)</li>
     * </ol>
     *
     * @param frequency Frequency in gigahertz
     * @param latitude Latitude in degrees
     * @param longitude Longitude in degrees
     * @param exceedanceProbability p in percent
     * @param elevationAngle Elevation angle (&theta;) in degrees
     * @param gridsByProbability A map where each key is a p, and the corresponding value is the grid of L(p).
     * @return Attenuation in dB
     */
    public static double computeSlantPathStatisticalCloudAttenuation(
            double frequency,
            double latitude,
            double longitude,
            double exceedanceProbability,
            double elevationAngle,
            TreeMap<Double, double[][]> gridsByProbability) {

        double probabilityBelow = gridsByProbability.floorKey(exceedanceProbability);
        double probabilityAbove = gridsByProbability.ceilingKey(exceedanceProbability);

        if (probabilityBelow == probabilityAbove) {
            double integratedCloudLiquidWaterContent = bilinearInterpolation(latitude, longitude, gridsByProbability.get(probabilityBelow));

            return computeSlantPathInstantaneousCloudAttenuation(frequency, integratedCloudLiquidWaterContent, elevationAngle);
        }

        double integratedCloudLiquidWaterContentBelow = bilinearInterpolation(latitude, longitude, gridsByProbability.get(probabilityBelow));
        double integratedCloudLiquidWaterContentAbove = bilinearInterpolation(latitude, longitude, gridsByProbability.get(probabilityAbove));

        double logP = Math.log10(exceedanceProbability);
        double logPBelow = Math.log10(probabilityBelow);
        double logPAbove = Math.log10(probabilityAbove);

        double integratedCloudLiquidWaterContent =
               integratedCloudLiquidWaterContentBelow +
               (integratedCloudLiquidWaterContentAbove - integratedCloudLiquidWaterContentBelow) * (logP - logPBelow) / (logPAbove - logPBelow);

        return computeSlantPathInstantaneousCloudAttenuation(frequency, integratedCloudLiquidWaterContent, elevationAngle);
    }

    /**
     * <p>Log-normal approximation to the slant path statistical cloud attenuation</p>
     * <p>Defined in ITU-R P.840-9, Equation (15).</p>
     *
     * <p>For p &lt; P<sub>L</sub>: A<sub>C</sub>(f, p) =
     *    K<sub>L</sub>(f) · exp(m<sub>L</sub> + &sigma;<sub>L</sub> · Q<sup>-1</sup>(p / P<sub>L</sub>)) / sin(&theta;)</p>
     * <p>For p &ge; P<sub>L</sub>: A<sub>C</sub>(f, p) = 0</p>
     *
     * <p><b>Note:</b> If the location is at a grid point where P<sub>L</sub> &le; 0.02%, or if the location
     *    is between grid points where any of the four surrounding grid points has P<sub>L</sub> &le; 0.02%,
     *    then A<sub>C</sub>(f, p) = 0 dB.</p>
     *
     * @param frequency Frequency in gigahertz
     * @param latitude Latitude in degrees
     * @param longitude Longitude in degrees
     * @param exceedanceProbability p in percent
     * @param elevationAngle Elevation angle (&theta;) in degrees
     * @param logNormalMeanParameterGrid m<sub>L</sub> grid in natural log
     * @param logNormalStandardDeviationParameterGrid &sigma;<sub>L</sub> grid in natural log
     * @param cloudProbabilityGrid P<sub>L</sub> grid in percent
     * @return Attenuation in dB
     */
    public static double computeLogNormalApproximationToTheSlantPathStatisticalCloudAttenuation(
            double frequency,
            double latitude,
            double longitude,
            double exceedanceProbability,
            double elevationAngle,
            double[][] logNormalMeanParameterGrid,
            double[][] logNormalStandardDeviationParameterGrid,
            double[][] cloudProbabilityGrid) {

        if (isAnyCornerCloudProbabilityBelowThreshold(latitude, longitude, cloudProbabilityGrid)) {
            return 0.0;
        }

        double logNormalMeanParameter = bilinearInterpolation(latitude, longitude, logNormalMeanParameterGrid);
        double logNormalStandardDeviationParameter = bilinearInterpolation(latitude, longitude, logNormalStandardDeviationParameterGrid);
        double cloudProbability = bilinearInterpolation(latitude, longitude, cloudProbabilityGrid);

        if (cloudProbability <= CLOUD_PROBABILITY_THRESHOLD_PERCENT || exceedanceProbability >= cloudProbability) {
            return 0.0;
        }

        double cloudLiquidMassAbsorptionCoefficient = computeCloudLiquidMassAbsorptionCoefficient(frequency);
        double inverseStandardNormalCCDF = computeInverseStandardNormalCCDF(exceedanceProbability / cloudProbability);
        double sineOfElevationAngle = Math.sin(Math.toRadians(elevationAngle));

        return cloudLiquidMassAbsorptionCoefficient *
               Math.exp(logNormalMeanParameter + logNormalStandardDeviationParameter * inverseStandardNormalCCDF) / sineOfElevationAngle;
    }

    // ==================================================================================
    //                           File Reading and Grid Creation
    // ==================================================================================

    /**
     * <p>Reads a digital map of L, m<sub>L</sub>, &sigma;<sub>L</sub>, or P<sub>L</sub>
     *    and constructs the corresponding grid.</p>
     *
     * <p>Grid points with <code>NaN</code> are preserved as {@link Double#NaN} in the grid.</p>
     *
     * @param filePath Path to the digital map file
     * @return A grid
     * @throws IOException If:
     * <ul>
     *   <li>A row is missing or empty,</li>
     *   <li>A row contains fewer or more columns than expected,</li>
     *   <li>A value cannot be parsed as a number,</li>
     *   <li>The file contains more rows than expected, or</li>
     *   <li>The file cannot otherwise be read.</li>
     * </ul>
     */
    public static double[][] readGridFromFile(String filePath) throws IOException {
        double[][] grid = new double[NUMBER_OF_LATITUDE_POINTS][NUMBER_OF_LONGITUDE_POINTS];

        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            for (int row = 0; row < NUMBER_OF_LATITUDE_POINTS; row++) {
                String line = br.readLine();

                if (line == null || line.trim().isEmpty()) {
                    throw new IOException("Row " + (row + 1) + " is missing.");
                }

                String[] values = line.trim().split("\\s+");

                if (values.length != NUMBER_OF_LONGITUDE_POINTS) {
                    if (values.length < NUMBER_OF_LONGITUDE_POINTS) {
                        throw new IOException("Row " + (row + 1) + " has fewer columns than expected.");
                    }
                    else {
                        throw new IOException("Row " + (row + 1) + " has more columns than expected.");
                    }
                }

                for (int column = 0; column < NUMBER_OF_LONGITUDE_POINTS; column++) {
                    try {
                        grid[row][column] = Double.parseDouble(values[column]);
                    }
                    catch (NumberFormatException e) {
                        throw new IOException("Invalid number at row " + (row + 1) + ", column " + (column + 1) + ": " + values[column] + ".", e);
                    }
                }
            }

            if (br.readLine() != null) {
                throw new IOException("File has more rows than expected.");
            }
        }

        return grid;
    }

    /**
     * <p>Loads the ITU-R P.840-9 digital maps of L from the specified package
     *    (<code>L_*.TXT</code> files for annual or monthly statistics) and constructs the corresponding grids.</p>
     *
     * <p>Returned map:</p>
     * <ul>
     *   <li><b>Key</b>: p</li>
     *   <li><b>Value</b>: The grid of L(p)</li>
     * </ul>
     *
     * <p>These grids are later used for statistical attenuation calculations (Equation 13).</p>
     *
     * @param folderPath Path to the digital maps folder
     * @return A map where each key is a p, and the corresponding value is the grid of L(p).
     * @throws IOException If:
     * <ul>
     *   <li>The folder is missing or empty,</li>
     *   <li>A file cannot be parsed into a grid,</li>
     *   <li>No valid <code>L_*.TXT</code> files are found, or</li>
     *   <li>An I/O error occurs while reading files.</li>
     * </ul>
     */
    public static TreeMap<Double, double[][]> loadGridsByProbability(String folderPath) throws IOException {
        TreeMap<Double, double[][]> gridsByProbability = new TreeMap<>();

        File folder = new File(folderPath);
        File[] files = folder.listFiles();

        if (files == null || files.length == 0) {
            throw new IOException("Folder is empty or missing: " + folderPath + ".");
        }

        for (File file : files) {
            String fileName = file.getName();

            if (!fileName.startsWith("L_") || !fileName.toUpperCase().endsWith(".TXT")) {
                continue;
            }

            String numericPart = fileName.substring(fileName.indexOf('_') + 1, fileName.lastIndexOf('.'));

            if (!numericPart.matches("\\d+")) {
                continue; // Skips files like L_mean.TXT and L_std.TXT
            }

            double exceedanceProbability;

            // 001, 002, 003 ,005 => 0.01%, 0.02%, 0.03%, 0.05%
            if (numericPart.length() == 3 && numericPart.startsWith("00")) {
                int val = Integer.parseInt(numericPart);
                exceedanceProbability = val / 100.0;
            }
            // 01, 02, 03, 05 => 0.1%, 0.2%, 0.3%, 0.5%
            else if (numericPart.length() == 2 && numericPart.startsWith("0")) {
                int val = Integer.parseInt(numericPart);
                exceedanceProbability = val / 10.0;
            }
            // 1, 2, 3, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99, 100 => 1%, 2%, 3%, 5%, 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, 95%, 99%, 100%
            else {
                int val = Integer.parseInt(numericPart);
                exceedanceProbability = val;
            }

            try {
                double[][] grid = readGridFromFile(file.getPath());
                gridsByProbability.put(exceedanceProbability, grid);
            }
            catch (IOException e) {
                throw new IOException("Failed to read " + fileName + " in " + folderPath + ".", e);
            }
        }

        if (gridsByProbability.isEmpty()) {
            throw new IOException("No valid L_*.TXT probability maps found in " + folderPath + ".");
        }

        return gridsByProbability;
    }

    // ==================================================================================
    //                                  Utility Functions
    // ==================================================================================

    /**
     * <p>Performs bilinear interpolation on a grid (per ITU-R P.840-9, Sections 4.2.1 and 4.2.2; P.1144 Annex 1).</p>
     *
     * <p><code>NaN</code> corner values are ignored, and the interpolation weights are renormalized to the valid neighbors.
     *    If all four corners are <code>NaN</code>, the function returns 0.</p>
     *
     * @param latitude Latitude in degrees
     * @param longitude Longitude in degrees
     * @param grid A grid
     * @return Interpolated value at (latitude, longitude), or 0 if no valid neighbors exist
     */
    public static double bilinearInterpolation(double latitude, double longitude, double[][] grid) {
        int southernLatitudeIndex = (int) Math.floor((latitude - GRID_LATITUDE_START_DEGREE) / GRID_LATITUDE_STEP_DEGREE);
        southernLatitudeIndex = Math.max(0, Math.min(southernLatitudeIndex, NUMBER_OF_LATITUDE_POINTS - 2));
        int northernLatitudeIndex = southernLatitudeIndex + 1;

        int westernLongitudeIndex = (int) Math.floor((longitude - GRID_LONGITUDE_START_DEGREE) / GRID_LONGITUDE_STEP_DEGREE);
        int easternLongitudeIndex = (westernLongitudeIndex + 1) % NUMBER_OF_LONGITUDE_POINTS;

        if (easternLongitudeIndex == 0) {
            easternLongitudeIndex = 1;
        }

        double southernLatitude = GRID_LATITUDE_START_DEGREE + southernLatitudeIndex * GRID_LATITUDE_STEP_DEGREE;
        double northernLatitude = GRID_LATITUDE_START_DEGREE + northernLatitudeIndex * GRID_LATITUDE_STEP_DEGREE;
        double westernLongitude = GRID_LONGITUDE_START_DEGREE + westernLongitudeIndex * GRID_LONGITUDE_STEP_DEGREE;
        double easternLongitude = GRID_LONGITUDE_START_DEGREE + easternLongitudeIndex * GRID_LONGITUDE_STEP_DEGREE;

        double xFraction = (longitude - westernLongitude) / (easternLongitude - westernLongitude);
        double yFraction = (latitude - southernLatitude) / (northernLatitude - southernLatitude);

        double valueSouthWest = grid[southernLatitudeIndex][westernLongitudeIndex];
        double valueSouthEast = grid[southernLatitudeIndex][easternLongitudeIndex];
        double valueNorthWest = grid[northernLatitudeIndex][westernLongitudeIndex];
        double valueNorthEast = grid[northernLatitudeIndex][easternLongitudeIndex];

        double weightSouthWest = (1 - xFraction) * (1 - yFraction);
        double weightSouthEast = xFraction * (1 - yFraction);
        double weightNorthWest = (1 - xFraction) * yFraction;
        double weightNorthEast = xFraction * yFraction;

        double weightedSum = 0.0;
        double totalWeight = 0.0;

        if (!Double.isNaN(valueSouthWest)) {
            weightedSum += valueSouthWest * weightSouthWest;
            totalWeight += weightSouthWest;
        }
        if (!Double.isNaN(valueSouthEast)) {
            weightedSum += valueSouthEast * weightSouthEast;
            totalWeight += weightSouthEast;
        }
        if (!Double.isNaN(valueNorthWest)) {
            weightedSum += valueNorthWest * weightNorthWest;
            totalWeight += weightNorthWest;
        }
        if (!Double.isNaN(valueNorthEast)) {
            weightedSum += valueNorthEast * weightNorthEast;
            totalWeight += weightNorthEast;
        }

        return (totalWeight > 0.0) ? (weightedSum / totalWeight) : 0.0;
    }

    /**
     * <p>Checks the four surrounding grid points of the specified location in the P<sub>L</sub> grid
     *    to determine whether the special condition in ITU-R P.840-9 applies.</p>
     *
     * <p>According to Equation (15), NOTE: If the desired location is at or between
     *    grid points where any P<sub>L</sub> &le; 0.02%, then the predicted attenuation
     *    A<sub>C</sub>(f, p) must be set to 0 dB.</p>
     *
     * @param latitude Latitude in degrees
     * @param longitude Longitude in degrees
     * @param cloudProbabilityGrid P<sub>L</sub> grid in percent
     * @return <code>true</code> If any of the four surrounding grid points has P<sub>L</sub> &le; 0.02%, else <code>false</code>
     */
    public static boolean isAnyCornerCloudProbabilityBelowThreshold(double latitude, double longitude, double[][] cloudProbabilityGrid) {
        int southernLatitudeIndex = (int) Math.floor((latitude - GRID_LATITUDE_START_DEGREE) / GRID_LATITUDE_STEP_DEGREE);
        southernLatitudeIndex = Math.max(0, Math.min(southernLatitudeIndex, NUMBER_OF_LATITUDE_POINTS - 2));
        int northernLatitudeIndex = southernLatitudeIndex + 1;

        int westernLongitudeIndex = (int) Math.floor((longitude - GRID_LONGITUDE_START_DEGREE) / GRID_LONGITUDE_STEP_DEGREE);
        int easternLongitudeIndex = (westernLongitudeIndex + 1) % NUMBER_OF_LONGITUDE_POINTS;

        if (easternLongitudeIndex == 0) {
            easternLongitudeIndex = 1;
        }

        double cloudProbabilitySouthWest = cloudProbabilityGrid[southernLatitudeIndex][westernLongitudeIndex];
        double cloudProbabilitySouthEast = cloudProbabilityGrid[southernLatitudeIndex][easternLongitudeIndex];
        double cloudProbabilityNorthWest = cloudProbabilityGrid[northernLatitudeIndex][westernLongitudeIndex];
        double cloudProbabilityNorthEast = cloudProbabilityGrid[northernLatitudeIndex][easternLongitudeIndex];

        return (cloudProbabilitySouthWest <= CLOUD_PROBABILITY_THRESHOLD_PERCENT) || (cloudProbabilitySouthEast <= CLOUD_PROBABILITY_THRESHOLD_PERCENT) ||
               (cloudProbabilityNorthWest <= CLOUD_PROBABILITY_THRESHOLD_PERCENT) || (cloudProbabilityNorthEast <= CLOUD_PROBABILITY_THRESHOLD_PERCENT);
    }

    /**
     * <p>Computes Q<sup>-1</sup>(x), the inverse standard normal CCDF.</p>
     * <p>Defined in ITU-R P.1057-7, Equations (5c)–(5e).</p>
     *
     * @param x A dimensionless argument between 0 and 1, equal to <code>p / P<sub>L</sub></code> in ITU-R P.840-9
     * @return Value y such that Q(y) = x for the standard normal CCDF
     */
    public static double computeInverseStandardNormalCCDF(double x) {
        final double EPSILON = 1e-16;
        double p = Math.max(EPSILON, Math.min(1.0 - EPSILON, x));

        // Coefficients from ITU-R P.1057-7 (Equations 5d and 5e)
        final double[] c = {
                2.938163982698783,
                4.374664141464968,
                -2.549732539343734,
                -2.400758277161838,
                -0.3223964580411365,
                -0.007784894002430293
        };
        final double[] d = {
                3.754408661907416,
                2.445134137142996,
                0.3224671290700398,
                0.007784695709041462
        };

        final double[] a = {
                2.506628277459239,
                -30.66479806614716,
                138.3577518672690,
                -275.9285104469687,
                220.9460984245205,
                -39.69683028665376
        };
        final double[] b = {
                -13.28068155288572,
                66.80131188771972,
                -155.6989798598866,
                161.5858368580409,
                -54.47609879822406
        };

        double y;

        if (p <= 0.5) {
            if (p <= 0.02425) {
                // Tail expansion (Equation 5d)
                double t = Math.sqrt(-2.0 * Math.log(p));
                y = (((((c[5] * t + c[4]) * t + c[3]) * t + c[2]) * t + c[1]) * t + c[0]) /
                    ((((d[3] * t + d[2]) * t + d[1]) * t + d[0]) * t + 1.0);
            }
            else {
                // Central expansion (Equation 5e)
                double difference = p - 0.5;
                double t = difference * difference;
                y = (((((a[5] * t + a[4]) * t + a[3]) * t + a[2]) * t + a[1]) * t + a[0]) * difference /
                    (((((b[4] * t + b[3]) * t + b[2]) * t + b[1]) * t + b[0]) * t + 1.0);
            }

            return -y;
        }
        else {
            p = 1.0 - p;

            if (p <= 0.02425) {
                // Tail expansion (Equation 5d)
                double t = Math.sqrt(-2.0 * Math.log(p));
                y = (((((c[5] * t + c[4]) * t + c[3]) * t + c[2]) * t + c[1]) * t + c[0]) /
                    ((((d[3] * t + d[2]) * t + d[1]) * t + d[0]) * t + 1.0);
            }
            else {
                // Central expansion (Equation 5e)
                double difference = p - 0.5;
                double t = difference * difference;
                y = (((((a[5] * t + a[4]) * t + a[3]) * t + a[2]) * t + a[1]) * t + a[0]) * difference /
                    (((((b[4] * t + b[3]) * t + b[2]) * t + b[1]) * t + b[0]) * t + 1.0);
            }

            return y;
        }
    }

    // ==================================================================================
    //                           Utility Functions for Testing
    // ==================================================================================

    /**
     * <p>Computes the log-normal term in Equation (15) in the ITU-R P.840-9.</p>
     *
     * <p>For p &lt; P<sub>L</sub>: Log-normal term = exp(m<sub>L</sub> + &sigma;<sub>L</sub> · Q<sup>-1</sup>(p / P<sub>L</sub>))</p>
     * <p>For p &ge; P<sub>L</sub>: Log-normal term = exp(m<sub>L</sub>)</p>
     *
     * @param exceedanceProbability p in percent
     * @param logNormalMeanParameter m<sub>L</sub> in natural log
     * @param logNormalStandardDeviationParameter &sigma;<sub>L</sub> in natural log
     * @param cloudProbability P<sub>L</sub> in percent
     * @return Log-normal term
     */
    public static double computeLogNormalTerm(double exceedanceProbability, double logNormalMeanParameter, double logNormalStandardDeviationParameter, double cloudProbability) {
        if (exceedanceProbability < cloudProbability) {
            double inverseStandardNormalCCDF = computeInverseStandardNormalCCDF(exceedanceProbability / cloudProbability);

            return Math.exp(logNormalMeanParameter + logNormalStandardDeviationParameter * inverseStandardNormalCCDF);
        }
        else {
            return Math.exp(logNormalMeanParameter);
        }
    }

    // ==================================================================================
    //                          Command Line Interface (Demo)
    // ==================================================================================

    /**
     * <p>Command Line Interface demo for the three prediction methods:</p>
     * <ol>
     *   <li>Slant path instantaneous cloud attenuation prediction method (Equation 11)</li>
     *   <li>Slant path statistical cloud attenuation prediction method (Equation 13) with annual or monthly statistics</li>
     *   <li>Log-normal approximation to the slant path statistical cloud attenuation (Equation 15) with annual log-normal statistics</li>
     * </ol>
     *
     * <p>This demo is designed for simple calculations and testing.</p>
     *
     * @param args Command line arguments (unused)
     */
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in).useLocale(Locale.US);

        System.out.println("Select slant path cloud attenuation prediction method:");
        System.out.println("1 - Instantaneous (Equation 11)");
        System.out.println("2 - Statistical (Equation 13) using annual/monthly statistics");
        System.out.println("3 - Log-normal approximation (Equation 15) using annual log-normal statistics");
        System.out.print("Enter your choice (1, 2, or 3): ");

        String methodToken = scanner.next().trim();

        if (!methodToken.matches("[123]")) {
            System.err.println("Invalid method choice. Please enter exactly 1, 2, or 3.");

            return;
        }

        int methodChoice = Integer.parseInt(methodToken);

        System.out.printf("Enter latitude in degrees [%.0f, %.0f]: ", LATITUDE_MIN_DEGREES, LATITUDE_MAX_DEGREES);
        double latitude = scanner.nextDouble();

        if (latitude < LATITUDE_MIN_DEGREES || latitude > LATITUDE_MAX_DEGREES) {
            System.err.printf("Latitude must be within [%.0f, %.0f] degrees.%n", LATITUDE_MIN_DEGREES, LATITUDE_MAX_DEGREES);

            return;
        }

        System.out.printf("Enter longitude in degrees [%.0f, %.0f]: ", LONGITUDE_MIN_DEGREES, LONGITUDE_MAX_DEGREES);
        double longitude = scanner.nextDouble();

        if (longitude < LONGITUDE_MIN_DEGREES || longitude > LONGITUDE_MAX_DEGREES) {
            System.err.printf("Longitude must be within [%.0f, %.0f] degrees.%n", LONGITUDE_MIN_DEGREES, LONGITUDE_MAX_DEGREES);

            return;
        }

        double exceedanceProbability = 0.0;

        if (methodChoice == 2 || methodChoice == 3) {
            System.out.printf(Locale.US, "Enter exceedance probability in percent [0, %.0f]: ", P_MAX_PERCENT);
            exceedanceProbability = scanner.nextDouble();

            if (exceedanceProbability < 0 || exceedanceProbability > P_MAX_PERCENT) {
                System.err.printf(Locale.US, "Exceedance probability must be within [0, %.0f] percent.%n", P_MAX_PERCENT);

                return;
            }
        }

        System.out.printf("Enter frequency in GHz [%.0f, %.0f]: ", FREQUENCY_MIN_GHZ, FREQUENCY_MAX_GHZ);
        double frequency = scanner.nextDouble();

        if (frequency < FREQUENCY_MIN_GHZ || frequency > FREQUENCY_MAX_GHZ) {
            System.err.printf("Frequency must be within [%.0f, %.0f] GHz.%n", FREQUENCY_MIN_GHZ, FREQUENCY_MAX_GHZ);

            return;
        }

        System.out.printf("Enter elevation angle (θ) in degrees (%.0f, %.0f]: ", ELEVATION_ANGLE_MIN_DEGREES, ELEVATION_ANGLE_MAX_DEGREES);
        double elevationAngle = scanner.nextDouble();

        if (elevationAngle <= ELEVATION_ANGLE_MIN_DEGREES || elevationAngle > ELEVATION_ANGLE_MAX_DEGREES) {
            System.err.printf("Elevation angle (θ) must be in (%.0f, %.0f] degrees.%n", ELEVATION_ANGLE_MIN_DEGREES, ELEVATION_ANGLE_MAX_DEGREES);

            return;
        }

        switch (methodChoice) {
            case 1: {
                System.out.print("Enter integrated cloud liquid water content in kg/m^2 or mm: ");
                double integratedCloudLiquidWaterContent = scanner.nextDouble();

                double attenuation = computeSlantPathInstantaneousCloudAttenuation(frequency, integratedCloudLiquidWaterContent, elevationAngle);
                System.out.printf("Instantaneous attenuation: %.10f dB%n", attenuation);

                break;
            }

            case 2: {
                System.out.print("Select statistics type: Annual (A) or Monthly (M): ");
                String statisticsType = scanner.next().trim().toUpperCase(Locale.ROOT);

                String folderPath;

                if ("A".equals(statisticsType)) {
                    folderPath = "data/annual/";
                }
                else if ("M".equals(statisticsType)) {
                    System.out.print("Enter month number (1-12): ");
                    int month = scanner.nextInt();

                    if (month < 1 || month > 12) {
                        System.err.println("Month must be between 1 and 12.");

                        return;
                    }

                    folderPath = String.format("data/month%02d/", month);
                }
                else {
                    System.err.println("Please enter 'A' or 'M'.");

                    return;
                }

                double P_MIN_PERCENT = "A".equals(statisticsType) ? P_MIN_PERCENT_ANNUAL : P_MIN_PERCENT_MONTHLY;

                if (exceedanceProbability < P_MIN_PERCENT) {
                    System.out.printf(Locale.US,
                                      "\u001B[33mWarning: The entered exceedance probability = %.3f%% is below the minimum allowed (%.2f%%). It has been adjusted to %.2f%%.\u001B[0m%n",
                                      exceedanceProbability, P_MIN_PERCENT, P_MIN_PERCENT);

                    exceedanceProbability = P_MIN_PERCENT;
                }

                TreeMap<Double, double[][]> gridsByProbability;

                try {
                    gridsByProbability = loadGridsByProbability(folderPath);
                }
                catch (IOException io) {
                    System.err.println("Failed to load grids by probability from " + folderPath + ".");
                    System.err.println(io.getMessage());

                    return;
                }

                double attenuation = computeSlantPathStatisticalCloudAttenuation(
                        frequency,
                        latitude,
                        longitude,
                        exceedanceProbability,
                        elevationAngle,
                        gridsByProbability);
                System.out.printf("Statistical attenuation: %.10f dB%n", attenuation);

                break;
            }

            case 3: {
                String folderPath = "data/logNormalAnnual/";

                double[][] logNormalMeanParameterGrid;
                double[][] logNormalStandardDeviationParameterGrid;
                double[][] cloudProbabilityGrid;

                try {
                    logNormalMeanParameterGrid = readGridFromFile(folderPath + "mL.TXT");
                    logNormalStandardDeviationParameterGrid = readGridFromFile(folderPath + "sL.TXT");
                    cloudProbabilityGrid = readGridFromFile(folderPath + "PL.TXT");
                }
                catch (IOException io) {
                    System.err.println("Failed to load log-normal parameter grids from " + folderPath + ".");
                    System.err.println(io.getMessage());

                    return;
                }

                double attenuation = computeLogNormalApproximationToTheSlantPathStatisticalCloudAttenuation(
                        frequency,
                        latitude,
                        longitude,
                        exceedanceProbability,
                        elevationAngle,
                        logNormalMeanParameterGrid,
                        logNormalStandardDeviationParameterGrid,
                        cloudProbabilityGrid);
                System.out.printf("Log-normal attenuation: %.10f dB%n", attenuation);

                break;
            }
        }
    }
}