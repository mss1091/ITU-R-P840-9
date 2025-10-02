package itu840;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.*;
import java.util.*;
import java.nio.charset.*;
import java.nio.file.*;

import static org.junit.jupiter.api.Assertions.*;

import java.time.Month;
import java.time.format.TextStyle;
import java.util.Locale;

/**
 * <p>Validator for the ITU-R P.840-9 monthly L(p) values computed with bilinear interpolation,
 *    using the digital maps in <code>data/monthXX/</code> for February, May, August, and November.</p>
 *
 * <p>To validate the implementation, this test uses <code>src/test/resources/monthlyLpTest.csv</code>,
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
 *   <li><b>Expected results:</b>
 *     <ul>
 *       <li>L(p) for February (kg/m<sup>2</sup>)</li>
 *       <li>L(p) for May (kg/m<sup>2</sup>)</li>
 *       <li>L(p) for August (kg/m<sup>2</sup>)</li>
 *       <li>L(p) for November (kg/m<sup>2</sup>)</li>
 *     </ul>
 *   </li>
 * </ol>
 */
public class MonthlyLpTest {

    private static final String CSV_RESOURCE = "/monthlyLpTest.csv";
    private static final Path DATA_DIRECTORY = Path.of("data");

    private static final int[] MONTHS = {2, 5, 8, 11};

    private static final double ABSOLUTE_TOLERANCE = 1e-6;
    private static final double RELATIVE_TOLERANCE = 5e-4;

    private static final HashMap<Integer, TreeMap<Double, double[][]>> monthlyGridsByProbability = new HashMap<>();

    @BeforeAll
    static void loadMonthlyDigitalMaps() throws Exception {
        for (int month : MONTHS) {
            Path DIGITAL_MAPS_DIRECTORY = DATA_DIRECTORY.resolve(String.format("month%02d", month));

            TreeMap<Double, double[][]> gridsByProbability = Itu840.loadGridsByProbability(DIGITAL_MAPS_DIRECTORY.toString());
            assertFalse(gridsByProbability.isEmpty(), "Failed to load monthly digital maps from " + DIGITAL_MAPS_DIRECTORY + ".");

            monthlyGridsByProbability.put(month, gridsByProbability);
        }
    }

    @Test
    void validateTestData() throws Exception {
        InputStream CSVFile = MonthlyLpTest.class.getResourceAsStream(CSV_RESOURCE);
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
                assertEquals(7, columns.length, "Test case " + testIndex + " does not have the expected number of columns.");

                // ===================== Inputs =====================
                double latitude = Double.parseDouble(columns[0]);
                double longitude = Double.parseDouble(columns[1]);
                double exceedanceProbability = Double.parseDouble(columns[2]);

                for (int i = 0; i < MONTHS.length; i++) {
                    int month = MONTHS[i];

                    // ===================== Expected value per month =====================
                    double integratedCloudLiquidWaterContentExpected = Double.parseDouble(columns[i + 3]);

                    // ===================== Computed value from Itu840 =====================
                    double integratedCloudLiquidWaterContent = Itu840.computeSlantPathStatisticalCloudAttenuation(
                            1,
                            latitude,
                            longitude,
                            exceedanceProbability,
                            90.0,
                            monthlyGridsByProbability.get(month)) / Itu840.computeCloudLiquidMassAbsorptionCoefficient(1);

                    System.out.printf(Locale.US, "Test %d (%s)%n", testIndex, Month.of(month).getDisplayName(TextStyle.FULL, Locale.ENGLISH));
                    System.out.printf(Locale.US, "Latitude = %.3f°, Longitude = %.3f°, p = %.3f%%%n",
                                      latitude, longitude, exceedanceProbability);
                    System.out.printf(Locale.US, "L(p) = %.10f (expected = %.10f)%n", integratedCloudLiquidWaterContent, integratedCloudLiquidWaterContentExpected);
                    System.out.println("-----------------------------------------------------------------------");

                    assertAll("Test " + testIndex,
                            () -> assertEquals(integratedCloudLiquidWaterContentExpected, integratedCloudLiquidWaterContent,
                                    Math.max(ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE * Math.abs(integratedCloudLiquidWaterContentExpected)))
                    );
                }

                testIndex++;
            }
        }
    }
}