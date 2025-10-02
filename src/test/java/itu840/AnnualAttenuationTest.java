package itu840;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.*;
import java.util.*;
import java.nio.charset.*;
import java.nio.file.*;

import static org.junit.jupiter.api.Assertions.*;

/**
 * <p>Validator for the ITU-R P.840-9 slant path statistical cloud attenuation prediction method (Equation 13),
 *    using the digital maps in <code>data/annual/</code>.</p>
 *
 * <p>To validate the implementation, this test uses <code>src/test/resources/annualAttenuationTest.csv</code>,
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
 *       <li>&epsilon;&prime;, The real part of the complex permittivity of liquid water (dimensionless)</li>
 *       <li>&epsilon;&Prime;, The imaginary part of the complex permittivity of liquid water (dimensionless)</li>
 *       <li>&eta;, A parameter derived from the real and imaginary parts of the permittivity (dimensionless)</li>
 *       <li>K<sub>L</sub>, The cloud liquid mass absorption coefficient (dB/(kg/m<sup>2</sup>))</li>
 *       <li>L(p) (kg/m<sup>2</sup>)</li>
 *     </ul>
 *   </li>
 *   <li><b>Expected final result:</b>
 *     <ul>
 *       <li>A<sub>C</sub>, Attenuation (dB)</li>
 *     </ul>
 *   </li>
 * </ol>
 */
public class AnnualAttenuationTest {

    private static final String CSV_RESOURCE = "/annualAttenuationTest.csv";
    private static final Path DIGITAL_MAPS_DIRECTORY = Path.of("data", "annual");

    private static final double ABSOLUTE_TOLERANCE = 1e-6;
    private static final double RELATIVE_TOLERANCE = 5e-4;

    private static TreeMap<Double, double[][]> gridsByProbability;

    @BeforeAll
    static void loadAnnualDigitalMaps() throws Exception {
        gridsByProbability = Itu840.loadGridsByProbability(DIGITAL_MAPS_DIRECTORY.toString());
        assertFalse(gridsByProbability.isEmpty(), "Failed to load annual digital maps from " + DIGITAL_MAPS_DIRECTORY + ".");
    }

    @Test
    void validateTestData() throws Exception {
        InputStream CSVFile = AnnualAttenuationTest.class.getResourceAsStream(CSV_RESOURCE);
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
                assertEquals(11, columns.length, "Test case " + testIndex + " does not have the expected number of columns.");

                // ===================== Inputs =====================
                double latitude = Double.parseDouble(columns[0]);
                double longitude = Double.parseDouble(columns[1]);
                double exceedanceProbability = Double.parseDouble(columns[2]);
                double frequency = Double.parseDouble(columns[3]);
                double elevationAngle = Double.parseDouble(columns[4]);

                // ===================== Expected values =====================
                double epsilonRealExpected = Double.parseDouble(columns[5]);
                double epsilonImaginaryExpected = Double.parseDouble(columns[6]);
                double etaExpected = Double.parseDouble(columns[7]);
                double cloudLiquidMassAbsorptionCoefficientExpected = Double.parseDouble(columns[8]);
                double integratedCloudLiquidWaterContentExpected = Double.parseDouble(columns[9]);
                double attenuationExpected = Double.parseDouble(columns[10]);

                // ===================== Computed values from Itu840 =====================
                double epsilonReal = Itu840.computeEpsilonReal(frequency);
                double epsilonImaginary = Itu840.computeEpsilonImaginary(frequency);
                double eta = Itu840.computeEta(frequency);
                double cloudLiquidMassAbsorptionCoefficient = Itu840.computeCloudLiquidMassAbsorptionCoefficient(frequency);
                double integratedCloudLiquidWaterContent = Itu840.computeSlantPathStatisticalCloudAttenuation(
                        frequency,
                        latitude,
                        longitude,
                        exceedanceProbability,
                        90.0,
                        gridsByProbability) / cloudLiquidMassAbsorptionCoefficient;
                double attenuation = Itu840.computeSlantPathStatisticalCloudAttenuation(
                        frequency,
                        latitude,
                        longitude,
                        exceedanceProbability,
                        elevationAngle,
                        gridsByProbability);

                System.out.printf(Locale.US, "Test %d%n", testIndex);
                System.out.printf(Locale.US, "Latitude = %.3f°, Longitude = %.3f°, p = %.3f%%, Frequency = %.3f GHz, Elevation Angle (θ) = %.3f°%n",
                                  latitude, longitude, exceedanceProbability, frequency, elevationAngle);
                System.out.printf(Locale.US, "ε′   = %.10f (expected = %.10f)%n", epsilonReal, epsilonRealExpected);
                System.out.printf(Locale.US, "ε″   = %.10f (expected = %.10f)%n", epsilonImaginary, epsilonImaginaryExpected);
                System.out.printf(Locale.US, "Eta  = %.10f (expected = %.10f)%n", eta, etaExpected);
                System.out.printf(Locale.US, "K_L  = %.10f (expected = %.10f)%n", cloudLiquidMassAbsorptionCoefficient, cloudLiquidMassAbsorptionCoefficientExpected);
                System.out.printf(Locale.US, "L(p) = %.10f (expected = %.10f)%n", integratedCloudLiquidWaterContent, integratedCloudLiquidWaterContentExpected);
                System.out.printf(Locale.US, "A_C  = %.10f (expected = %.10f)%n", attenuation, attenuationExpected);
                System.out.println("-----------------------------------------------------------------------");

                assertAll("Test " + testIndex,
                        () -> assertEquals(epsilonRealExpected, epsilonReal,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(epsilonRealExpected))),

                        () -> assertEquals(epsilonImaginaryExpected, epsilonImaginary,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(epsilonImaginaryExpected))),

                        () -> assertEquals(etaExpected, eta,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(etaExpected))),

                        () -> assertEquals(cloudLiquidMassAbsorptionCoefficientExpected, cloudLiquidMassAbsorptionCoefficient,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(cloudLiquidMassAbsorptionCoefficientExpected))),

                        () -> assertEquals(integratedCloudLiquidWaterContentExpected, integratedCloudLiquidWaterContent,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(integratedCloudLiquidWaterContentExpected))),

                        () -> assertEquals(attenuationExpected, attenuation,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(attenuationExpected)))
                );

                testIndex++;
            }
        }
    }
}