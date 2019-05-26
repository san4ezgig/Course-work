
import java.awt.*;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.ListIterator;

import core.GeneratorCreator;
import kurs.BigDecimalMatrix;
import core.*;
import org.knowm.xchart.*;

import javax.swing.*;

import static java.math.BigDecimal.*;


public class Main {

    private static double cKor(BigDecimal lambda, BigDecimalMatrix d0, BigDecimalMatrix d1) {
        BigDecimal lambdaT = new BigDecimal(1 / lambda.doubleValue());
        BigDecimalMatrix tetta = getTetta(d0, d1);
        BigDecimalMatrix d0TNeg = (d0.multiply(new BigDecimal(-1).setScale(12, RoundingMode.HALF_UP))).inverse();
        BigDecimalMatrix e = BigDecimalMatrix.eCol(2, ONE);
        BigDecimal v = new BigDecimal(tetta.multiply(lambda.multiply(new BigDecimal(2))).multiply(d0TNeg).multiply(e).getElement(0, 0)
                .subtract(new BigDecimal(-1)).doubleValue() / lambda.pow(2).doubleValue());
        System.out.println("Cvar");
        System.out.println((tetta.multiply(lambda.multiply(new BigDecimal(2))).multiply(d0TNeg).multiply(e).getElement(0, 0)
                .subtract(new BigDecimal(-1))));
        return (tetta.multiply(lambdaT).multiply(d0TNeg).multiply(d1).multiply(d0TNeg).multiply(e).getElement(0, 0)
                .subtract(lambdaT.pow(2))).doubleValue() / v.doubleValue();
    }

    private static BigDecimalMatrix getTetta(BigDecimalMatrix d0, BigDecimalMatrix d1) {
        BigDecimalMatrix D = d0.add(d1);
        D.setElement(0, 0, new BigDecimal(1).setScale(12, RoundingMode.HALF_UP));
        D.setElement(1, 0, new BigDecimal(1).setScale(12, RoundingMode.HALF_UP));
        return new BigDecimalMatrix(new BigDecimal[][]{
                {new BigDecimal(1), new BigDecimal(0)}
        }, 6).multiply(D.inverse());
    }

    private static BigDecimal sqrt(BigDecimal value) {
        BigDecimal x = new BigDecimal(Math.sqrt(value.doubleValue()));
        return x.add(new BigDecimal(value.subtract(x.multiply(x)).doubleValue() / (x.doubleValue() * 2.0)));
    }

    private static String getNameOfGraphics(int number) {
        switch (number) {
            case 11: {
                return "L";
            }
            case 12: {
                return "N";
            }
            case 13: {
                return "P0 custom";
            }
            case 14: {
                return "P0 energy";
            }
            case 15: {
                return "P idle";
            }
            case 16: {
                return "~L";
            }
            default: {
                return "";
            }
        }
    }

    private static void DrawGraphics(ArrayList<Double>[] data, ArrayList<String> nameOfLines, final int numberOfPerfomanceFunction, final int systemSize) {
        String nameOfGraphics = getNameOfGraphics(numberOfPerfomanceFunction);
        final XYChart chart = new XYChartBuilder().width(800).height(800).title(nameOfGraphics).xAxisTitle("K").yAxisTitle(nameOfGraphics).build();

        ListIterator<String> iterator = nameOfLines.listIterator();
        for (ArrayList<Double> list : data) {
            double[] arrayOfPoints = new double[list.size()];
            Double[] array = list.toArray(new Double[list.size()]);
            for (int i = 0; i < array.length; i++) {
                arrayOfPoints[i] = array[i];
            }

            chart.addSeries(iterator.next(), arrayOfPoints);
        }

        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {

                // Create and set up the window.
                JFrame frame = new JFrame("Advanced Example");
                frame.setLayout(new BorderLayout());
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

                // chart
                JPanel chartPanel = new XChartPanel<XYChart>(chart);
                frame.add(chartPanel, BorderLayout.CENTER);

                // label
                JLabel label = new JLabel("Blah blah blah.", SwingConstants.CENTER);
                frame.add(label, BorderLayout.SOUTH);

                // Display the window.
                frame.pack();
                frame.setVisible(true);
                try {
                    //BitmapEncoder.saveBitmap(chart, "./" + (numberOfPerfomanceFunction + (systemSize == 1 ? 10 : 0)), BitmapEncoder.BitmapFormat.JPG);

                } catch (Exception e) {
                    System.out.println(1);

                }
            }
        });
    }

    public static void main(String[] args) {
        //Variable sector
        BigDecimal gamma = new BigDecimal("0.4");
        BigDecimal lambda = (new BigDecimal("0.3"));
        BigDecimal HALF = new BigDecimal("0.5");
        BigDecimal accuracy = new BigDecimal("0.000001");
        int systemSize = 2;
        BigDecimalMatrix d0, d1;
        switch (systemSize) {
            case 1: {
                d0 = new BigDecimalMatrix(1, lambda.negate(), 8);
                d1 = new BigDecimalMatrix(1, lambda, 8);
                break;
            }
            case 2: {
                d0 = new BigDecimalMatrix(new BigDecimal[][]{
                        {new BigDecimal("-0.405780"), ZERO},
                        {ZERO, new BigDecimal("-0.013173")}
                }, 8);
                d1 = new BigDecimalMatrix(new BigDecimal[][]{
                        {new BigDecimal("0.403080"), new BigDecimal("0.002700")},
                        {new BigDecimal("0.007338"), new BigDecimal("0.005835")}
                }, 8);
                break;
            }
            default: {
                d0 = new BigDecimalMatrix(1, lambda.negate(), 8);
                d1 = new BigDecimalMatrix(1, lambda, 8);
            }
        }
        BigDecimalMatrix matrixExponent = new BigDecimalMatrix(new BigDecimal[][]{
                {new BigDecimal("0.998653"), new BigDecimal("0.00134662")},
                {new BigDecimal("0.00365981"), new BigDecimal("0.99634")},
        }, 8);
        System.out.println(d0);
        System.out.println(d1);
        System.out.println(d0.add(d1).multiply(HALF));
        // Initital value sector
        System.out.println("tetta:");
        System.out.println(getTetta(d0, d1));
        // System.out.println(cKor(lambda, d0, d1));
        System.out.println();
        Test.runAllTests(gamma, lambda, d0, d1, accuracy, matrixExponent, systemSize);
        System.out.println(ErgodicityCondition.getLambdaRestriction(ONE, gamma, 10));

        MatrixContainer.reInit();
        GMatrixCreator gMatrixCreator;
        PSlashMatrixCreator pSlashMatrixCreator;
        PhiMatrixCreator phiMatrixCreator;
        StationaryDistributionCreator sdCreator;
        PerformanceParameters pParameters;
        GeneratorCreator generatorCreator;
        BigDecimalMatrix g0;
        ArbitararyTimeGenerator arbitararyTimeGenerator;

        ArrayList<Double>[] paramsValueList = new ArrayList[4];
        ArrayList<String> nameOfLines = new ArrayList<>();
        int numberOfPerfomanceFunction = 11;
        //PP sector
        /*for (int n = 0; n < 6; n++) {
            for (BigDecimal t = HALF, i = ZERO; t.compareTo(valueOf(2)) <= 0; t = t.add(HALF), i = i.add(ONE)) {
                paramsValueList[i.intValue()] = new ArrayList<>();
                nameOfLines.add(t.toString());
                System.out.println();
                for (int k = 1; k <= 20; k++) {
                    BigDecimal K = new BigDecimal(k);
                    g0 = BigDecimalMatrix.identity(systemSize * (K.intValue() + 1));
                    generatorCreator = new GeneratorCreator(gamma, lambda, K.intValue(), d0, d1, t, accuracy, systemSize);
                    gMatrixCreator = new GMatrixCreator(generatorCreator);
                    pSlashMatrixCreator = new PSlashMatrixCreator(generatorCreator, gMatrixCreator.create(g0));
                    phiMatrixCreator = new PhiMatrixCreator(pSlashMatrixCreator, K.intValue(), systemSize);
                    sdCreator = new StationaryDistributionCreator(pSlashMatrixCreator, phiMatrixCreator.getPhiMatrices(), K.intValue(), systemSize);
                    ArrayList<BigDecimalMatrix> piVector = sdCreator.getPiVectors();
                    arbitararyTimeGenerator = new ArbitararyTimeGenerator(generatorCreator, piVector);
                    pParameters = new PerformanceParameters(piVector, lambda, d1, systemSize, arbitararyTimeGenerator);
                    BigDecimal param = pParameters.requestCorrectPerfomanceFunction(numberOfPerfomanceFunction, k);
                    System.out.println(param.toString().substring(0, 11));
                    paramsValueList[i.intValue()].add(param.doubleValue());
                    MatrixContainer.reInit();
                }
                System.out.println(t + " Completed");
            }
            DrawGraphics(paramsValueList, nameOfLines, numberOfPerfomanceFunction, systemSize);
            numberOfPerfomanceFunction++;
        }*/

        for (int n = 0; n < 6; n++) {
            for (BigDecimal t = HALF, i = ZERO; t.compareTo(valueOf(2)) <= 0; t = t.add(HALF), i = i.add(ONE)) {
                int k = 10;
                BigDecimal lambdaRestriction = valueOf(ErgodicityCondition.getLambdaRestriction(t, gamma, k)).subtract(valueOf(0.0001));
                paramsValueList[i.intValue()] = new ArrayList<>();
                nameOfLines.add(t.toString());
                System.out.println();
                for (BigDecimal newLambda = lambdaRestriction.subtract(valueOf(0.0020)); newLambda.compareTo(lambdaRestriction) <= 0; newLambda = newLambda.add(valueOf(0.0001))) {
                    BigDecimal coeff = valueOf(newLambda.doubleValue() / lambda.doubleValue());
                    // System.out.println(coeff);
                    BigDecimalMatrix d0New = d0.multiply(coeff);
                    BigDecimalMatrix d1New = d1.multiply(coeff);
                    //System.out.println(newLambda);
                    BigDecimal K = new BigDecimal(k);
                    g0 = BigDecimalMatrix.identity(systemSize * (K.intValue() + 1));
                    generatorCreator = new GeneratorCreator(gamma, newLambda, K.intValue(), d0New, d1New, t, accuracy, systemSize);
                    gMatrixCreator = new GMatrixCreator(generatorCreator);
                    pSlashMatrixCreator = new PSlashMatrixCreator(generatorCreator, gMatrixCreator.create(g0));
                    phiMatrixCreator = new PhiMatrixCreator(pSlashMatrixCreator, K.intValue(), systemSize);
                    sdCreator = new StationaryDistributionCreator(pSlashMatrixCreator, phiMatrixCreator.getPhiMatrices(), K.intValue(), systemSize);
                    ArrayList<BigDecimalMatrix> piVector = sdCreator.getPiVectors();
                    arbitararyTimeGenerator = new ArbitararyTimeGenerator(generatorCreator, piVector);
                    pParameters = new PerformanceParameters(piVector, newLambda, d1New, systemSize, arbitararyTimeGenerator);
                    BigDecimal param = pParameters.requestCorrectPerfomanceFunction(numberOfPerfomanceFunction, k);
                    System.out.println(param.toString().substring(0, 11));
                    paramsValueList[i.intValue()].add(param.doubleValue());
                    MatrixContainer.reInit();
                }
                System.out.println(t + " Completed");
            }
            DrawGraphics(paramsValueList, nameOfLines, numberOfPerfomanceFunction, systemSize);
            numberOfPerfomanceFunction++;
        }
        //
        /*for (int n = 0; n < 5; n++) {
            for (BigDecimal coeff = valueOf(1), i = ZERO; coeff.compareTo(valueOf(2)) <= 0; coeff = coeff.add(valueOf(0.25)), i = i.add(ONE)) {
                BigDecimal newGamma = gamma.multiply(coeff);
                BigDecimal newLambda = lambda.multiply(coeff);
                BigDecimalMatrix newD0 = d0.multiply(coeff);
                BigDecimalMatrix newD1 = d1.multiply(coeff);
                paramsValueList[i.intValue()] = new ArrayList<>();
                BigDecimal t = valueOf(1);
                nameOfLines.add(lambda.multiply(coeff).toString());
                System.out.println();
                for (int k = 1; k <= 20; k++) {
                    BigDecimal K = new BigDecimal(k);
                    g0 = BigDecimalMatrix.identity(systemSize * (K.intValue() + 1));
                    generatorCreator = new GeneratorCreator(newGamma, newLambda, K.intValue(), newD0, newD1, t, accuracy, systemSize);
                    gMatrixCreator = new GMatrixCreator(generatorCreator);
                    pSlashMatrixCreator = new PSlashMatrixCreator(generatorCreator, gMatrixCreator.create(g0));
                    phiMatrixCreator = new PhiMatrixCreator(pSlashMatrixCreator, K.intValue(), systemSize);
                    sdCreator = new StationaryDistributionCreator(pSlashMatrixCreator, phiMatrixCreator.getPhiMatrices(), K.intValue(), systemSize);
                    ArrayList<BigDecimalMatrix> piVector = sdCreator.getPiVectors();
                    arbitararyTimeGenerator = new ArbitararyTimeGenerator(generatorCreator, piVector);
                    pParameters = new PerformanceParameters(piVector, newLambda, newD1, systemSize, arbitararyTimeGenerator);
                    BigDecimal param = pParameters.requestCorrectPerfomanceFunction(numberOfPerfomanceFunction);
                    System.out.println(param.toString().substring(0, 11));
                    paramsValueList[i.intValue()].add(param.doubleValue());
                    MatrixContainer.reInit();
                }
                System.out.println(coeff + " Completed");
            }
            DrawGraphics(paramsValueList, nameOfLines, numberOfPerfomanceFunction);
            numberOfPerfomanceFunction++;
        }*/

        /*for (BigDecimal t = new BigDecimal("2.3610800"); t.compareTo(new BigDecimal("2.3610840")) <= 0; t = t.add(new BigDecimal("0.0000002"))) {
            System.out.println();
            for (int k = 1; k < 21; k++) {
                System.out.println(t + " " + k);
                if (!ErgodicityCondition.check(gamma, lambda, k, t)) {
                    throw new IllegalArgumentException("Ergodicity condition does not perform. " + t + " " + k);
                }
            }
            System.out.println(t + " Completed");
        }*/
        //Extremal condition test
        /*for (BigDecimal t = new BigDecimal("2.3610800"); t.compareTo(new BigDecimal("2.3610840")) <= 0; t = t.add(new BigDecimal("0.000001"))) {
            System.out.println();
            for (int k = 1; k < 21; k++) {
                BigDecimal K = new BigDecimal(k);
                g0 = BigDecimalMatrix.identity(2 * K.intValue() + 2);
                generatorCreator = new GeneratorCreator(gamma, lambda, K.intValue(), d0, d1, t, accuracy, systemSize);
                gMatrixCreator = new GMatrixCreator(generatorCreator);
                pSlashMatrixCreator = new PSlashMatrixCreator(generatorCreator, gMatrixCreator.create(g0));
                phiMatrixCreator = new PhiMatrixCreator(pSlashMatrixCreator, K.intValue(), systemSize);
                sdCreator = new StationaryDistributionCreator(pSlashMatrixCreator, phiMatrixCreator.getPhiMatrices(), K.intValue(), systemSize);
                ArrayList<BigDecimalMatrix> piVector = sdCreator.getPiVectors();
                arbitararyTimeGenerator = new ArbitararyTimeGenerator(generatorCreator, piVector);
                pParameters = new PerformanceParameters(piVector, lambda, d1, systemSize, arbitararyTimeGenerator);
                // list[k].add((pParameters.getAverageNumberOfRequests().toString().substring(0, 11)));
                System.out.println(pParameters.getAverageNumberOfRequests().toString().substring(0, 11));
                MatrixContainer.reInit();
            }
            System.out.println(t + " Completed");
        }*/
        //Eps grapher
        //GNU Plot
    }
}
