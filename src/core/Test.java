package core;

import kurs.BigDecimalMatrix;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Objects;

import static java.math.BigDecimal.*;

public class Test {
    public static void runAllTests(BigDecimal gamma,
                                   BigDecimal lambda,
                                   BigDecimalMatrix d0,
                                   BigDecimalMatrix d1,
                                   BigDecimal accuracy,
                                   BigDecimalMatrix matrixExponent,
                                   int systemSize
    ) {
        BigDecimal HALF = new BigDecimal("0.5");
        GMatrixCreator gMatrixCreator;
        PSlashMatrixCreator pSlashMatrixCreator;
        PhiMatrixCreator phiMatrixCreator;
        StationaryDistributionCreator sdCreator;
        PerformanceParameters pParameters;
        GeneratorCreator generatorCreator;
        BigDecimalMatrix g0;

        MatrixContainer.reInit();
        BigDecimal tForCheck = HALF;
        System.out.println();
        int j = 0;
        BigDecimal K = new BigDecimal("5");
        BigDecimalMatrix vMatrix;
        generatorCreator = new GeneratorCreator(gamma, lambda, K.intValue(), d0, d1, tForCheck, accuracy, systemSize);

        BigDecimalMatrix result = BigDecimalMatrix.eCol(systemSize * (K.intValue() + 1), ZERO);
        do {
            vMatrix = generatorCreator.create(0, j).multiply(BigDecimalMatrix.eCol(systemSize * (K.intValue() + 1), ONE));
            result = result.add(vMatrix);
            j++;
        } while (vMatrix.squaredEuclidianNorm().doubleValue() > accuracy.doubleValue());
        System.out.println("Vmatrix:");
        System.out.println(result);
        System.out.println(BigDecimalMatrix.eRow(systemSize * (K.intValue() + 1), ONE).multiply(BigDecimalMatrix.eCol(systemSize * (K.intValue() + 1), ONE).subtract(result)));
        System.out.println();

        int l = 0;
        BigDecimalMatrix yMatrix;
        result = BigDecimalMatrix.eCol(systemSize * (K.intValue() + 1), ZERO);
        do {
            yMatrix = generatorCreator.create(1, l).multiply(BigDecimalMatrix.eCol(systemSize * (K.intValue() + 1), ONE));
            result = result.add(yMatrix);
            l++;
        } while (yMatrix.squaredEuclidianNorm().doubleValue() > accuracy.doubleValue());
        System.out.println("Ymatrix:");
        System.out.println(result);
        System.out.println(BigDecimalMatrix.eRow(systemSize * (K.intValue() + 1), ONE).multiply(BigDecimalMatrix.eCol(systemSize * (K.intValue() + 1), ONE).subtract(result)));

        result = BigDecimalMatrix.eCol(systemSize, ZERO);
        BigDecimalMatrix mMatrix;
        j = 0;
        do {
            mMatrix = generatorCreator.funcM(j).multiply(BigDecimalMatrix.eCol(systemSize, ONE));
            result = result.add(mMatrix);
            j++;
        } while (mMatrix.squaredEuclidianNorm().doubleValue() > accuracy.doubleValue());
        System.out.println("M:");
        System.out.println(result);

        result = BigDecimalMatrix.eCol(systemSize, ZERO);
        BigDecimalMatrix nMatrix;
        j = 0;
        do {
            nMatrix = generatorCreator.funcN(j).multiply(BigDecimalMatrix.eCol(systemSize, ONE));
            result = result.add(nMatrix);
            j++;
        } while (nMatrix.squaredEuclidianNorm().doubleValue() > accuracy.doubleValue());
        System.out.println("N:");
        System.out.println(result);

        BigDecimal notMatrixResult = ZERO;
        BigDecimal smalPhiK;
        j = 0;
        do {
            smalPhiK = GeneratorCreator.funcSmallPhiK(gamma, HALF, j);
            notMatrixResult = notMatrixResult.add(smalPhiK);
            j++;
        } while (smalPhiK.doubleValue() > accuracy.doubleValue());
        System.out.println("phi");
        System.out.println(notMatrixResult);

        System.out.println("P:");
        System.out.println(checkP(accuracy, tForCheck, generatorCreator, systemSize));
        System.out.println("P should be");
        switch (systemSize) {
            case 1: {
                System.out.println(CheckPForStationaryDist(lambda, tForCheck, accuracy));
                System.out.println(checkP(accuracy, tForCheck, generatorCreator, systemSize).getElement(0, 0).subtract(CheckPForStationaryDist(lambda, tForCheck, accuracy)));
            }
            case 2: {
                System.out.println(matrixExponent);
                System.out.println(checkP(accuracy, tForCheck, generatorCreator, systemSize).subtract(matrixExponent));
            }
        }
        System.out.println("Phi");
        System.out.println(checkPhi(generatorCreator, K.intValue(), accuracy, systemSize));
        System.out.println("aaaaaaa");
//        System.out.println(PForStationaryDist(lambda, tForCheck, 5));
//        System.out.println(generatorCreator.funcP(5, tForCheck));
//        System.out.println(PForStationaryDist(lambda, tForCheck, 5).subtract(generatorCreator.funcP(5, tForCheck).getElement(0, 0)));
        MatrixContainer.reInit();
        //Check system
        g0 = BigDecimalMatrix.identity(systemSize * K.intValue() + systemSize);
        gMatrixCreator = new GMatrixCreator(generatorCreator);
        pSlashMatrixCreator = new PSlashMatrixCreator(generatorCreator, gMatrixCreator.create(g0));
        phiMatrixCreator = new PhiMatrixCreator(pSlashMatrixCreator, K.intValue(), systemSize);
        sdCreator = new StationaryDistributionCreator(pSlashMatrixCreator, phiMatrixCreator.getPhiMatrices(), K.intValue(), systemSize);
        sdCreator.checkPIMatrix(tForCheck, gamma, d0, lambda);

        //MatrixContainer.reInit();
        System.out.println("test arbitry time");
        ArrayList<BigDecimalMatrix> piVector = sdCreator.getPiVectors();
        ArbitararyTimeGenerator arbitararyTimeGenerator = new ArbitararyTimeGenerator(generatorCreator, piVector);
        result = BigDecimalMatrix.eRow(systemSize, ZERO);
        for (j = 0; j < piVector.size(); j++) {
            result = result.add(arbitararyTimeGenerator.calculateP(j).multiply(calculateKronekerMatrix(K.intValue(), systemSize)));
        }
        System.out.println(result);

        System.out.println("test arbitry time vector");
        result = BigDecimalMatrix.eRow(1, ZERO);
        for (j = 0; j < piVector.size(); j++) {
            result = result.add(arbitararyTimeGenerator.calculateP(j).multiply(BigDecimalMatrix.eCol(systemSize * (K.intValue() + 1), ONE)));
        }
        System.out.println(result);
//        if (systemSize == 2) {
//            System.out.println(result.multiply(d0.add(d1)));
//        }

        MatrixContainer.reInit();
        System.out.println("test Y matrix in arbitry time");
        arbitararyTimeGenerator = new ArbitararyTimeGenerator(generatorCreator, piVector);
        result = BigDecimalMatrix.eCol(systemSize * (K.intValue() + 1), ZERO);
        j = 0;
        BigDecimalMatrix eCol = BigDecimalMatrix.eCol(systemSize * (K.intValue() + 1), ONE);
        BigDecimalMatrix value;
        do {
            value = arbitararyTimeGenerator.create(1, j + 1).multiply(eCol);
            result = result.add(value);
            j++;
        } while (value.squaredEuclidianNorm().doubleValue() > accuracy.doubleValue());
        System.out.println(result);

        result = BigDecimalMatrix.zeroMatrix(systemSize);
        BigDecimalMatrix sum;
        int n = 0;
        do {
            sum = BigDecimalMatrix.zeroMatrix(systemSize);
            for (int k = 0; k <= K.intValue(); k++) {
                sum = sum.add(arbitararyTimeGenerator.funcPhiWithWave(n, k));
            }
            result = result.add(sum);
            n++;
        } while (sum.squaredEuclidianNorm().doubleValue() > accuracy.doubleValue());

        System.out.println("PhiWithWave:");
        System.out.println(result.multiply(BigDecimalMatrix.eCol(systemSize, ONE)));
        switch (systemSize) {
            case 1: {
                System.out.println("PhiWithWaveCheck:");
                BigDecimal phiWithWaveCheck = tForCheck;
                System.out.println(phiWithWaveCheck);
                break;
            }
            case 2: {
                BigDecimal phiWithWaveCheck = tForCheck;
                System.out.println("PhiWithWaveCheck:");
                System.out.println(phiWithWaveCheck);
            }
        }
        System.out.println(arbitararyTimeGenerator.calculateP(0).multiply(BigDecimalMatrix.eCol(systemSize * (K.intValue() + 1), ONE)));

        MatrixContainer.reInit();
    }

    private static BigDecimalMatrix calculateKronekerMatrix(int K, int systemSize) {
        BigDecimalMatrix result = new BigDecimalMatrix((K + 1) * systemSize, systemSize, 12);
        BigDecimalMatrix I = BigDecimalMatrix.identity(systemSize);
        for (int i = 0; i < K + 1; i++) {
            switch (systemSize) {
                case 1: {
                    result.setElement(i, 0, I.getElement(0, 0));
                    break;
                }
                case 2: {
                    result.setElement(i * systemSize, 0, I.getElement(0, 0));
                    result.setElement(i * systemSize, 1, I.getElement(0, 1));
                    result.setElement(1 + i * systemSize, 0, I.getElement(1, 0));
                    result.setElement(1 + i * systemSize, 1, I.getElement(1, 1));
                }
            }
        }
        return result;
    }

    private static BigDecimal PForStationaryDist(BigDecimal lambda, BigDecimal T, int k) {
        return new BigDecimal(lambda.multiply(T).pow(k).doubleValue() / GeneratorCreator.factor(k)
                * Math.exp(-lambda.multiply(T).doubleValue()));
    }

    private static BigDecimal CheckPForStationaryDist(BigDecimal lambda, BigDecimal T, BigDecimal accuracy) {
        BigDecimal sum = ZERO;
        BigDecimal val;
        int i = 0;
        do {
            val = PForStationaryDist(lambda, T, i);
            sum = sum.add(val);
            i++;
        } while (val.doubleValue() >= accuracy.doubleValue());
        return sum;
    }

    public static BigDecimalMatrix checkP(BigDecimal accuracy, BigDecimal T, GeneratorCreator generatorCreator, int systemSize) {
        BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, systemSize, 12);
        BigDecimalMatrix val;
        int i = 0;
        do {
            val = generatorCreator.funcP(i, T);
            sum = sum.add(val);
            i++;
        } while (val.squaredEuclidianNorm().doubleValue() >= accuracy.doubleValue());
        return sum;
    }


    public static BigDecimalMatrix checkPhi(GeneratorCreator generatorCreator, int K, BigDecimal accuracy, int systemSize) {
        BigDecimalMatrix sum = BigDecimalMatrix.eCol(systemSize, ZERO);
        BigDecimalMatrix e = BigDecimalMatrix.eCol(systemSize, ONE);
        BigDecimalMatrix val;
        int i = 0;
        do {
            val = BigDecimalMatrix.eCol(systemSize, ZERO);
            for (int k = 0; k <= K; k++) {
                val = val.add(generatorCreator.funcPhi(i, k).multiply(e));
            }
            sum = sum.add(val);
            i++;
        } while (val.squaredEuclidianNorm().compareTo(accuracy) > 0);
        return sum;
    }
}
