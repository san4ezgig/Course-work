package core;

import kurs.BigDecimalMatrix;

import java.math.BigDecimal;

import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;

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
        BigDecimal K = new BigDecimal("10");
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
        System.out.println(result);

        result = BigDecimalMatrix.eCol(systemSize, ZERO);
        BigDecimalMatrix nMatrix;
        j = 0;
        do {
            nMatrix = generatorCreator.funcN(j).multiply(BigDecimalMatrix.eCol(systemSize, ONE));
            result = result.add(nMatrix);
            j++;
        } while (nMatrix.squaredEuclidianNorm().doubleValue() > accuracy.doubleValue());
        System.out.println(result);

        BigDecimal notMatrixResult = ZERO;
        BigDecimal smalPhiK;
        j = 0;
        do {
            smalPhiK = GeneratorCreator.funcSmallPhiK(gamma, HALF, j);
            notMatrixResult = notMatrixResult.add(smalPhiK);
            j++;
        } while (smalPhiK.doubleValue() > accuracy.doubleValue());
        System.out.println(notMatrixResult);

        System.out.println(Test.checkP(accuracy, tForCheck, generatorCreator, systemSize));
        System.out.println(matrixExponent);
        if (systemSize == 2) {
            System.out.println(Test.checkP(accuracy, tForCheck, generatorCreator, systemSize).subtract(matrixExponent));
        }
        System.out.println(Test.checkPhi(generatorCreator, K.intValue(), accuracy, systemSize));
        MatrixContainer.reInit();
        //Check system
        g0 = BigDecimalMatrix.identity(systemSize * K.intValue() + systemSize);
        gMatrixCreator = new GMatrixCreator(generatorCreator);
        pSlashMatrixCreator = new PSlashMatrixCreator(generatorCreator, gMatrixCreator.create(g0));
        phiMatrixCreator = new PhiMatrixCreator(pSlashMatrixCreator, K.intValue(), systemSize);
        sdCreator = new StationaryDistributionCreator(pSlashMatrixCreator, phiMatrixCreator.getPhiMatrices(), K.intValue(), systemSize);
        if (systemSize == 2) {
            sdCreator.checkPIMatrix(tForCheck, gamma, d0, lambda);
        }
        MatrixContainer.reInit();
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
