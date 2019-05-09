
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;

import core.GeneratorCreator;
import kurs.BigDecimalMatrix;
import core.*;

import static java.math.BigDecimal.*;


public class Main {

    static double cKor(BigDecimal lambda, BigDecimalMatrix d0, BigDecimalMatrix d1) {
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

    static BigDecimalMatrix getTetta(BigDecimalMatrix d0, BigDecimalMatrix d1) {
        BigDecimalMatrix D = d0.add(d1);
        D.setElement(0, 0, new BigDecimal(1).setScale(12, RoundingMode.HALF_UP));
        D.setElement(1, 0, new BigDecimal(1).setScale(12, RoundingMode.HALF_UP));
        return new BigDecimalMatrix(new BigDecimal[][]{
                {new BigDecimal(1), new BigDecimal(0)}
        }, 6).multiply(D.inverse());
    }

    static BigDecimal sqrt(BigDecimal value) {
        BigDecimal x = new BigDecimal(Math.sqrt(value.doubleValue()));
        return x.add(new BigDecimal(value.subtract(x.multiply(x)).doubleValue() / (x.doubleValue() * 2.0)));
    }

    private static void showList(ArrayList<String>[] list, Boolean isExcel) {
        if (isExcel) {
            for (int i = 0; i < 4; i++) {
                for (ArrayList item : list) {
                    System.out.println(item.get(i));
                }
                System.out.println();
            }
        } else {
            for (ArrayList item : list) {
                for (Object k : item) {
                    System.out.print(k + " ");
                }
                System.out.println();
            }
        }
    }

    public static void main(String[] args) {
        //Variable sector
        BigDecimal gamma = new BigDecimal("0.4");
        BigDecimal lambda = (new BigDecimal("0.3"));
        BigDecimal HALF = new BigDecimal("0.5");
        int scale = 20;
        BigDecimal accuracy = new BigDecimal("0.00000001");

        BigDecimalMatrix d0 = new BigDecimalMatrix(new BigDecimal[][]{
                {new BigDecimal("-0.405780"), ZERO},
                {ZERO, new BigDecimal("-0.013173")}
        }, 8);
        BigDecimalMatrix d1 = new BigDecimalMatrix(new BigDecimal[][]{
                {new BigDecimal("0.403080"), new BigDecimal("0.002700")},
                {new BigDecimal("0.007338"), new BigDecimal("0.005835")}
        }, 8);
        BigDecimalMatrix matrixExponent = new BigDecimalMatrix(new BigDecimal[][]{
                {new BigDecimal("0.998653"), new BigDecimal("0.00134662")},
                {new BigDecimal("0.00365981"), new BigDecimal("0.99634")},
        }, 8);
        // Initital value sector
        BigDecimalMatrix tetta = getTetta(d0, d1);
        System.out.println(d0);
        System.out.println(d1);
        System.out.println(d0.add(d1).multiply(HALF));
        System.out.println(tetta);
        System.out.println(cKor(lambda, d0, d1));
        System.out.println();
        Test.runAllTests(gamma, lambda, d0, d1, accuracy, matrixExponent);

        MatrixContainer.reInit();
        GMatrixCreator gMatrixCreator;
        PSlashMatrixCreator pSlashMatrixCreator;
        PhiMatrixCreator phiMatrixCreator;
        StationaryDistributionCreator sdCreator;
        PerformanceParameters pParameters;
        GeneratorCreator generatorCreator;
        BigDecimalMatrix g0;
        //PP sector
        for (BigDecimal t = HALF; t.compareTo(valueOf(2)) <= 0; t = t.add(HALF)) {
            System.out.println();
            for (int k = 1; k <= 20; k++) {
                BigDecimal K = new BigDecimal(k);
                g0 = BigDecimalMatrix.identity(2 * K.intValue() + 2);
                generatorCreator = new GeneratorCreator(gamma, lambda, K.intValue(), d0, d1, t, accuracy);
                gMatrixCreator = new GMatrixCreator(generatorCreator);
                pSlashMatrixCreator = new PSlashMatrixCreator(generatorCreator, gMatrixCreator.create(g0));
                phiMatrixCreator = new PhiMatrixCreator(pSlashMatrixCreator, K.intValue());
                sdCreator = new StationaryDistributionCreator(pSlashMatrixCreator, phiMatrixCreator.getPhiMatrices(), K.intValue());
                pParameters = new PerformanceParameters(sdCreator.getPiVectors(), lambda, d1);
                System.out.println(pParameters.getNoEnergySystemIdleProbability().toString().substring(0, 11));
               // pParameters.check(tetta);
                MatrixContainer.reInit();
            }
            System.out.println(t + " Completed");
        }
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
                generatorCreator = new GeneratorCreator(gamma, lambda, K.intValue(), d0, d1, t, accuracy);
                gMatrixCreator = new GMatrixCreator(generatorCreator);
                pSlashMatrixCreator = new PSlashMatrixCreator(generatorCreator, gMatrixCreator.create(g0));
                phiMatrixCreator = new PhiMatrixCreator(pSlashMatrixCreator, K.intValue());
                sdCreator = new StationaryDistributionCreator(pSlashMatrixCreator, phiMatrixCreator.getPhiMatrices(), K.intValue());
                pParameters = new PerformanceParameters(sdCreator.getPiVectors(), lambda, d1);
                // list[k].add((pParameters.getAverageNumberOfRequests().toString().substring(0, 11)));
                System.out.println(pParameters.getAverageNumberOfRequests().toString().substring(0, 11));
                MatrixContainer.reInit();
            }
            System.out.println(t + " Completed");
        }*/
        //Eps grapher
        //GNU Plot


        /*System.out.println();
        System.out.println(d0.add(d1).multiply(T));
        System.out.println(generatorCreator.checkPhi());*//*
       // System.out.println(generatorCreator.checkP());
        *//*BigDecimalMatrix[][] checkMatrix = generatorCreator.create(2, 2);*/

        //System.out.println(pParameters.getAverageNumberOfEnergyUnits());

        /*for (BigDecimal t = HALF; t.compareTo(T) <= 0; t = t.add(HALF)) {
            System.out.println("T = " + t);
            System.out.println();
            for (double l = 0.1; l < 0.4; l += 0.02) {
                BigDecimal l1 = new BigDecimal(l);
                d0 = d0.multiply(l1.divide(lambda, RoundingMode.HALF_UP))-;
                d1 = d1.multiply(l1.divide(lambda, RoundingMode.HALF_UP));
                K = new BigDecimal(5);
                g0 = BigDecimalMatrix.identity(2 * K.intValue() + 2);
                generatorCreator = new GeneratorCreator(gamma, l1, K.intValue(), d0, d1, t, accuracy);
                gMatrixCreator = new GMatrixCreator(generatorCreator);
                pSlashMatrixCreator = new PSlashMatrixCreator(generatorCreator, gMatrixCreator.create(g0));
                phiMatrixCreator = new PhiMatrixCreator(pSlashMatrixCreator, K.intValue());
                sdCreator = new StationaryDistributionCreator(pSlashMatrixCreator, phiMatrixCreator.getPhiMatrices(), K.intValue());
                pParameters = new PerformanceParameters(sdCreator.getPiVectors(), lambda);
                System.out.println(pParameters.getAverageNumberOfRequests().toString().substring(0, 11));
                MatrixContainer.reInit();
            }
            System.out.println();
        }*/
    }
}
