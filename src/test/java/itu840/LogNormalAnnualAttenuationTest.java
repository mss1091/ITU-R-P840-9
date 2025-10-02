package itu840;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.*;
import java.util.*;
import java.nio.charset.*;
import java.nio.file.*;

import static org.junit.jupiter.api.Assertions.*;

/**
 * <p>Validator for the ITU-R P.840-9 log-normal approximation to the slant path statistical cloud attenuation,
 *    using the digital maps in <code>data/logNormalAnnual/</code>.</p>
 *
 * <p>To validate the implementation, this test uses <code>src/test/resources/logNormalAnnualAttenuationTest.csv</code>,
 *    which contains the following columns:</p>
 *
 * <ol>
 *   <li><b>Inputs:</b>
 *     <ul>
 *       <li>Latitude (°N)</li>
 *       <li>Longitude (°E)</li>
 *       <li>p (%)</li>
 *       <li>Frequency (GHz)</li>
 *       <li>Elevation angle (&theta;) (°)</li>
 *     </ul>
 *   </li>
 *   <li><b>Expected intermediate values:</b>
 *     <ul>
 *       <li>K<sub>L</sub>, The cloud liquid mass absorption coefficient (dB/(kg/m<sup>2</sup>))</li>
 *       <li>m<sub>L</sub>, Log-normal mean parameter (natural log scale)</li>
 *       <li>&sigma;<sub>L</sub>, Log-normal standard deviation parameter (natural log scale)</li>
 *       <li>P<sub>L</sub>, Probability of cloud (%)</li>
 *       <li>Log-normal term, The exponential expression from Equation (15) (dimensionless)</li>
 *       <li>A<sub>C</sub> at zenith, Attenuation at θ = 90° (dB)</li>
 *     </ul>
 *   </li>
 *   <li><b>Expected final result:</b>
 *     <ul>
 *       <li>A<sub>C</sub>, Attenuation (dB)</li>
 *     </ul>
 *   </li>
 * </ol>
 */
public class LogNormalAnnualAttenuationTest {

    private static final String CSV_RESOURCE = "/logNormalAnnualAttenuationTest.csv";
    private static final Path DIGITAL_MAPS_DIRECTORY = Path.of("data", "logNormalAnnual");

    private static final double ABSOLUTE_TOLERANCE = 1e-6;
    private static final double RELATIVE_TOLERANCE = 5e-4;

    private static double[][] logNormalMeanParameterGrid, logNormalStandardDeviationParameterGrid, cloudProbabilityGrid;

    @BeforeAll
    static void loadAnnualDigitalMaps() throws Exception {
        logNormalMeanParameterGrid = Itu840.readGridFromFile(DIGITAL_MAPS_DIRECTORY.resolve("mL.TXT").toString());
        assertNotNull(logNormalMeanParameterGrid, "Failed to load annual log-normal mean parameter grid (mL.TXT) from " + DIGITAL_MAPS_DIRECTORY + ".");

        logNormalStandardDeviationParameterGrid = Itu840.readGridFromFile(DIGITAL_MAPS_DIRECTORY.resolve("sL.TXT").toString());
        assertNotNull(logNormalStandardDeviationParameterGrid, "Failed to annual load log-normal standard deviation parameter grid (sL.TXT) from " + DIGITAL_MAPS_DIRECTORY + ".");

        cloudProbabilityGrid = Itu840.readGridFromFile(DIGITAL_MAPS_DIRECTORY.resolve("PL.TXT").toString());
        assertNotNull(cloudProbabilityGrid, "Failed to load annual cloud probability grid (PL.TXT) from " + DIGITAL_MAPS_DIRECTORY + ".");
    }

    @Test
    void validateTestData() throws Exception {
        InputStream CSVFile = LogNormalAnnualAttenuationTest.class.getResourceAsStream(CSV_RESOURCE);
        assertNotNull(CSVFile, "CSV file could not be found at " + CSV_RESOURCE + ".");

        try (BufferedReader reader = new BufferedReader(new InputStreamReader(CSVFile, StandardCharsets.UTF_8))) {
            // Skip two header rows
            String header1 = reader.readLine();
            String header2 = reader.readLine();
            assertNotNull(header1, "First header row is missing.");
            assertNotNull(header2, "Second header row is missing.");

            String line;
            int testIndex = 1;

            while ((line = reader.readLine()) != null) {
                if (line.isBlank()) {
                    continue;
                }

                String[] columns = line.split(",", -1);
                assertEquals(12, columns.length, "Test case " + testIndex + " does not have the expected number of columns.");

                // ===================== Inputs =====================
                double latitude = Double.parseDouble(columns[0]);
                double longitude = Double.parseDouble(columns[1]);
                double exceedanceProbability = Double.parseDouble(columns[2]);
                double frequency = Double.parseDouble(columns[3]);
                double elevationAngle = Double.parseDouble(columns[4]);

                // ===================== Expected values =====================
                double cloudLiquidMassAbsorptionCoefficientExpected = Double.parseDouble(columns[5]);
                double logNormalMeanParameterExpected = Double.parseDouble(columns[6]);
                double logNormalStandardDeviationParameterExpected = Double.parseDouble(columns[7]);
                double cloudProbabilityExpected = Double.parseDouble(columns[8]);
                double logNormalTermExpected = Double.parseDouble(columns[9]);
                double attenuationZenithExpected = Double.parseDouble(columns[10]);
                double attenuationExpected = Double.parseDouble(columns[11]);

                // ===================== Computed values from Itu840 =====================
                double cloudLiquidMassAbsorptionCoefficient = Itu840.computeCloudLiquidMassAbsorptionCoefficient(frequency);
                double logNormalMeanParameter = Itu840.bilinearInterpolation(latitude, longitude, logNormalMeanParameterGrid);
                double logNormalStandardDeviationParameter = Itu840.bilinearInterpolation(latitude, longitude, logNormalStandardDeviationParameterGrid);
                double cloudProbability = Itu840.bilinearInterpolation(latitude, longitude, cloudProbabilityGrid);
                double logNormalTerm = Itu840.computeLogNormalTerm(exceedanceProbability, logNormalMeanParameter, logNormalStandardDeviationParameter, cloudProbability);
                double attenuationZenith = Itu840.computeLogNormalApproximationToTheSlantPathStatisticalCloudAttenuation(
                        frequency,
                        latitude,
                        longitude,
                        exceedanceProbability,
                        90,
                        logNormalMeanParameterGrid,
                        logNormalStandardDeviationParameterGrid,
                        cloudProbabilityGrid);
                double attenuation = Itu840.computeLogNormalApproximationToTheSlantPathStatisticalCloudAttenuation(
                        frequency,
                        latitude,
                        longitude,
                        exceedanceProbability,
                        elevationAngle,
                        logNormalMeanParameterGrid,
                        logNormalStandardDeviationParameterGrid,
                        cloudProbabilityGrid);

                System.out.printf(Locale.US, "Test %d%n", testIndex);
                System.out.printf(Locale.US, "Latitude = %.3f°, Longitude = %.3f°, p = %.3f%%, Frequency = %.3f GHz, Elevation Angle (θ) = %.3f°%n",
                                  latitude, longitude, exceedanceProbability, frequency, elevationAngle);
                System.out.printf(Locale.US, "K_L             = %.10f (expected = %.10f)%n", cloudLiquidMassAbsorptionCoefficient, cloudLiquidMassAbsorptionCoefficientExpected);
                System.out.printf(Locale.US, "m_L             = %.10f (expected = %.10f)%n", logNormalMeanParameter, logNormalMeanParameterExpected);
                System.out.printf(Locale.US, "s_L             = %.10f (expected = %.10f)%n", logNormalStandardDeviationParameter, logNormalStandardDeviationParameterExpected);
                System.out.printf(Locale.US, "P_L             = %.10f (expected = %.10f)%n", cloudProbability, cloudProbabilityExpected);
                System.out.printf(Locale.US, "Log-normal term = %.10f (expected = %.10f)%n", logNormalTerm, logNormalTermExpected);
                System.out.printf(Locale.US, "A_C at zenith   = %.10f (expected = %.10f)%n", attenuationZenith, attenuationZenithExpected);
                System.out.printf(Locale.US, "A_C             = %.10f (expected = %.10f)%n", attenuation, attenuationExpected);
                System.out.println("-----------------------------------------------------------------------");

                assertAll("Test " + testIndex,
                        () -> assertEquals(cloudLiquidMassAbsorptionCoefficientExpected, cloudLiquidMassAbsorptionCoefficient,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(cloudLiquidMassAbsorptionCoefficientExpected))),

                        () -> assertEquals(logNormalMeanParameterExpected, logNormalMeanParameter,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(logNormalMeanParameterExpected))),

                        () -> assertEquals(logNormalStandardDeviationParameterExpected, logNormalStandardDeviationParameter,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(logNormalStandardDeviationParameterExpected))),

                        () -> assertEquals(cloudProbabilityExpected, cloudProbability,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(cloudProbabilityExpected))),

                        () -> assertEquals(logNormalTermExpected, logNormalTerm,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(logNormalTermExpected))),

                        () -> assertEquals(attenuationZenithExpected, attenuationZenith,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(attenuationZenithExpected))),

                        () -> assertEquals(attenuationExpected, attenuation,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(attenuationExpected)))
                );

                testIndex++;
            }
        }
    }
}