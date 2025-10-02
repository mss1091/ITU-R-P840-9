package itu840;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.*;
import java.util.*;
import java.nio.charset.*;
import java.nio.file.*;

import static org.junit.jupiter.api.Assertions.*;

/**
 * <p>Validator for the ITU-R P.840-9 annual L(p) values computed with bilinear interpolation,
 *    using the digital maps in <code>data/annual/</code>.</p>
 *
 * <p>To validate the implementation, this test uses <code>src/test/resources/annualLpTest.csv</code>,
 *    which contains the following columns:</p>
 *
 * <ol>
 *   <li><b>Inputs:</b>
 *     <ul>
 *       <li>Latitude (°N)</li>
 *       <li>Longitude (°E)</li>
 *       <li>p (%)</li>
 *     </ul>
 *   </li>
 *   <li><b>Expected result:</b>
 *     <ul>
 *       <li>L(p) (kg/m<sup>2</sup>)</li>
 *     </ul>
 *   </li>
 * </ol>
 */
public class AnnualLpTest {

    private static final String CSV_RESOURCE = "/annualLpTest.csv";
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
        InputStream CSVFile = AnnualLpTest.class.getResourceAsStream(CSV_RESOURCE);
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
                assertEquals(4, columns.length, "Test case " + testIndex + " does not have the expected number of columns.");

                // ===================== Inputs =====================
                double latitude = Double.parseDouble(columns[0]);
                double longitude = Double.parseDouble(columns[1]);
                double exceedanceProbability = Double.parseDouble(columns[2]);

                // ===================== Expected value =====================
                double integratedCloudLiquidWaterContentExpected = Double.parseDouble(columns[3]);

                // ===================== Computed value from Itu840 =====================
                double integratedCloudLiquidWaterContent = Itu840.computeSlantPathStatisticalCloudAttenuation(
                        1,
                        latitude,
                        longitude,
                        exceedanceProbability,
                        90.0,
                        gridsByProbability) / Itu840.computeCloudLiquidMassAbsorptionCoefficient(1);

                System.out.printf(Locale.US, "Test %d%n", testIndex);
                System.out.printf(Locale.US, "Latitude = %.3f°, Longitude = %.3f°, p = %.3f%%%n",
                                  latitude, longitude, exceedanceProbability);
                System.out.printf(Locale.US, "L(p) = %.10f (expected = %.10f)%n", integratedCloudLiquidWaterContent, integratedCloudLiquidWaterContentExpected);
                System.out.println("-----------------------------------------------------------------------");

                assertAll("Test " + testIndex,
                        () -> assertEquals(integratedCloudLiquidWaterContentExpected, integratedCloudLiquidWaterContent,
                                Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(integratedCloudLiquidWaterContentExpected)))
                );

                testIndex++;
            }
        }
    }
}