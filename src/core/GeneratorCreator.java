package core;

import kurs.BigDecimalMatrix;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Objects;

import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;
import static java.math.BigDecimal.valueOf;


public class GeneratorCreator {
    public BigDecimal gamma;
    public BigDecimal lambda;
    public Integer K;
    public BigDecimal accuracy;
    public BigDecimalMatrix d0;
    public BigDecimalMatrix d1;
    public BigDecimal T;
    private MatrixElementCreator elementCreator;
    public BigDecimal tetta;
    public int systemSize;

    public GeneratorCreator(
            BigDecimal gamma,
            BigDecimal lambda,
            int K,
            BigDecimalMatrix d0,
            BigDecimalMatrix d1,
            BigDecimal T,
            BigDecimal accuracy,
            int systemSize
    ) {
        boolean validness = gamma.compareTo(ONE) < 0 &&
                gamma.compareTo(ZERO) > 0 &&
                lambda.compareTo(ONE) < 0 &&
                lambda.compareTo(ZERO) > 0 &&
                K > 0 &&
                gamma.compareTo(lambda) > 0;
        /*boolean validness = K > 0 &&
                gamma.compareTo(lambda) > 0;*/

        if (!ErgodicityCondition.check(gamma, lambda, K, T)) {
            throw new IllegalArgumentException("Ergodicity condition does not perform. " + T + " " + K);
        }
        if (!validness) {
            throw new IllegalArgumentException("Input parameters are not valid.");
        }

        this.gamma = gamma;
        this.lambda = lambda;
        this.K = K;
        this.accuracy = accuracy;
        this.d0 = d0;
        this.tetta = tetta();
        this.T = T;
        this.d1 = d1;
        this.elementCreator = new MatrixElementCreator();
        this.systemSize = systemSize;
    }

    public BigDecimalMatrix create(int i, int j) {
        BigDecimalMatrix tempMatrix, matrix = new BigDecimalMatrix(systemSize * (K + 1), systemSize * (K + 1), 12);

        if (i - j > 1) {
            return matrix;
        }

        int hash = i >= 1 ? Objects.hash(i, j - i + 1) : Objects.hash(0, j);

        // int hash = Objects.hash(i, j);

        try {
            if (MatrixContainer.getGenerators().containsKey(hash)) {
                matrix = MatrixContainer.getGenerators().get(hash).clone();
            } else {
                for (int k = 0; k < K + 1; k++) {
                    for (int bK = 0; bK < K + 1; bK++) {
                        tempMatrix = elementCreator.create(i, k, j, bK);
                        switch (systemSize) {
                            case 1: {
                                matrix.setElement(k, bK, tempMatrix.getElement(0, 0));
                                break;
                            }
                            case 2: {
                                matrix.setElement(k * systemSize, bK * systemSize, tempMatrix.getElement(0, 0));
                                matrix.setElement(k * systemSize, 1 + bK * systemSize, tempMatrix.getElement(0, 1));
                                matrix.setElement(1 + k * systemSize, bK * systemSize, tempMatrix.getElement(1, 0));
                                matrix.setElement(1 + k * systemSize, 1 + bK * systemSize, tempMatrix.getElement(1, 1));
                            }
                        }
                    }
                }
                MatrixContainer.getGenerators().put(hash, matrix.clone());
            }
        } catch (CloneNotSupportedException e) {
            System.out.println(e.getMessage());
        }
        return matrix;
    }

    private class MatrixElementCreator {

        private BigDecimalMatrix create(int i, int k, int j, int bK) {
            boolean validness = i >= 0 &&
                    j >= 0 &&
                    k >= 0 &&
                    k <= K &&
                    bK >= 0 &&
                    bK <= K;
            if (!validness) {
                throw new IllegalArgumentException("Indexes out of bounds.");
            }

            BigDecimalMatrix element;

            boolean isVMatrix = (i == 0);

            if (k - bK > 1) {
                element = BigDecimalMatrix.zeroMatrix(systemSize);
                return element;
            }

            if (isVMatrix) {
                if (k == 0) {
                    if (bK == K) {
                        return createFirstRowLastColumnElementV(j);
                    }

                    if (bK == K - 1) {
                        return createFirstRowPenultimateElementV(j);
                    }

                    return createFirstRowNotLastColumnElementV(j, bK);
                } else {
                    if (bK == K) {
                        return createNotFirstRowLastColumnElementV(k, j);
                    }

                    if (bK == K - 1) {
                        return createNotFirstRowPenultimateElementV(j, k);
                    }

                    return createNotFirstRowNotLastColumnElementV(k, j, bK);
                }
            }

            if (bK == K) {
                element = k == 0 ? createFirstRowLastColumnElement(j, i) : createLastColumnElement(k, j, i);
                return element;
            }

            if (k == 0) {
                element = createFirstRowNotLastColumnElement(j, bK, i);
                return element;
            }

            element = createNotFirstRowNotLastColumnElement(k, j, bK, i);
            return element;
        }

        // YMatrix region
        private BigDecimalMatrix createFirstRowNotLastColumnElement(int j, int bK, int i) {
            BigDecimalMatrix result = new BigDecimalMatrix(systemSize, systemSize, 10);
            for (int n = 0; n <= j - i + 1; n++) {
                result = result.add(funcN(n).multiply(funcPhi(j - i + 1 - n, bK)));
            }
            return result;
        }

        private BigDecimalMatrix createFirstRowLastColumnElement(int j, int i) {
            BigDecimalMatrix result = new BigDecimalMatrix(systemSize, systemSize, 10);
            for (int n = 0; n <= j - i + 1; n++) {
                result = result.add(funcN(n).multiply(funcPhiWithHat(j - i + 1 - n, K)));
            }
            return result;
        }

        private BigDecimalMatrix createLastColumnElement(int k, int j, int i) {
            return funcPhiWithHat(j - i + 1, K - k + 1);
        }

        private BigDecimalMatrix createNotFirstRowNotLastColumnElement(int k, int j, int bK, int i) {
            return funcPhi(j - i + 1, bK - k + 1);
        }
        // end Region

        //VMatrix region

        private BigDecimalMatrix createFirstRowLastColumnElementV(int j) {
            BigDecimalMatrix result = new BigDecimalMatrix(systemSize, systemSize, 12);
            for (int n = 0; n <= j; n++) {
                result = result.add(funcN(n).multiply(funcPhiWithHat(j - n, K)));
            }
            result = funcM(0).multiply(result);

            BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, systemSize, 12);
            for (int m = 0; m <= K - 1; m++) {
                sum = sum.add(funcM(m).multiply(funcPhiWithHat(j, K - m)));
            }
            sum = sum.add(funcMWithHat(K).multiply(funcPhiWithHat(j, 1)));
            sum = funcN(0).multiply(sum);

            result = result.add(sum);
            return result;
        }

        private BigDecimalMatrix createFirstRowPenultimateElementV(int j) {
            BigDecimalMatrix result = new BigDecimalMatrix(systemSize, systemSize, 12);
            for (int n = 0; n <= j; n++) {
                result = result.add(funcN(n).multiply(funcPhi(j - n, K - 1)));
            }
            result = funcM(0).multiply(result);

            BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, systemSize, 12);
            for (int m = 0; m <= K - 1; m++) {
                sum = sum.add(funcM(m).multiply(funcPhi(j, K - 1 - m)));
            }

            sum = sum.add(funcMWithHat(K).multiply(funcPhi(j, 0)));
            sum = funcN(0).multiply(sum);
            result = result.add(sum);
            return result;
        }

        private BigDecimalMatrix createFirstRowNotLastColumnElementV(int j, int bK) {
            BigDecimalMatrix result = new BigDecimalMatrix(systemSize, systemSize, 12);

            for (int n = 0; n <= j; n++) {
                result = result.add(funcN(n).multiply(funcPhi(j - n, bK)));
            }

            result = funcM(0).multiply(result);

            BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, systemSize, 12);

            for (int m = 0; m <= bK; m++) {
                sum = sum.add(funcM(m).multiply(funcPhi(j, bK - m)));
            }

            sum = funcN(0).multiply(sum);

            result = result.add(sum);

            return result;
        }

        private BigDecimalMatrix createNotFirstRowPenultimateElementV(int j, int k) {
            BigDecimalMatrix result = new BigDecimalMatrix(systemSize, systemSize, 12);
            for (int m = 0; m <= K - k; m++) {
                result = result.add(funcM(m).multiply(funcPhi(j, K - k - m)));
            }
            result = result.add(funcMWithHat(K - k + 1).multiply(funcPhi(j, 0)));
            return result;
        }

        private BigDecimalMatrix createNotFirstRowLastColumnElementV(int k, int j) {
            BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, systemSize, 12);
            for (int m = 0; m <= K - k; m++) {
                sum = sum.add(funcM(m).multiply(funcPhiWithHat(j, K - k + 1 - m)));
            }
            sum = sum.add(funcMWithHat(K - k + 1).multiply(funcPhiWithHat(j, 1)));
            return sum;
        }

        private BigDecimalMatrix createNotFirstRowNotLastColumnElementV(int k, int j, int bK) {
            BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, systemSize, 12);
            for (int m = 0; m <= bK - k + 1; m++) {
                sum = sum.add(funcM(m).multiply(funcPhi(j, bK - k + 1 - m)));
            }
            return sum;
        }
        // End region
    }

    public BigDecimalMatrix funcMWithHat(int n) {
        BigDecimalMatrix val;
        BigDecimalMatrix negativeD0 = d0.multiply(new BigDecimal(-1));
        BigDecimal gammaN = gamma.pow(n);
        val = (negativeD0.add(BigDecimalMatrix.identity(systemSize).multiply(gamma))).inverse();
        val = BigDecimalMatrix.powMatrix(val, n);
        val = val.multiply(gammaN).multiply(negativeD0.inverse()).multiply(d1);
        return val;
    }

    public BigDecimalMatrix funcM(int n) {
        BigDecimalMatrix val;
        BigDecimalMatrix negativeD0 = d0.multiply(ONE.negate());
        BigDecimal gammaN = gamma.pow(n);
        val = (negativeD0.add(BigDecimalMatrix.identity(systemSize).multiply(gamma))).inverse();
        val = BigDecimalMatrix.powMatrix(val, n + 1);
        val = val.multiply(gammaN).multiply(d1);
        return val;
    }

    public BigDecimalMatrix funcK(int n, int j) {
        try {
            int hash = Objects.hash(n + 10, j * 20);

            if (MatrixContainer.getKMatrix().containsKey(hash)) {
                return MatrixContainer.getKMatrix().get(hash).clone();
            }

            BigDecimal negTetta = ONE.divide(tetta, 12, RoundingMode.HALF_DOWN);
            if (j == 0) {
                BigDecimalMatrix value;
                if (n == 0) {
                    value = new BigDecimalMatrix(systemSize, ONE, 12);
                    MatrixContainer.getKMatrix().put(hash, value.clone());
                    return value;
                }

                value = BigDecimalMatrix.zeroMatrix(systemSize);
                MatrixContainer.getKMatrix().put(hash, value.clone());
                return value;
            }

            if (n == 0) {
                BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, ONE, 12).add(d0.multiply(negTetta));
                sum = funcK(0, j - 1).multiply(sum);
                MatrixContainer.getKMatrix().put(hash, sum.clone());
                return sum;
            }
            BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, systemSize, 12);

            if (n >= 1) {
                sum = sum.add(funcK(n - 1, j - 1).multiply(d1).multiply(negTetta));
            }

            sum = sum.add(funcK(n, j - 1)
                    .multiply(new BigDecimalMatrix(systemSize, ONE, 12).add(d0.multiply(negTetta))));
            MatrixContainer.getKMatrix().put(hash, sum.clone());
            return sum;

        } catch (CloneNotSupportedException e) {
            System.out.println(e.getMessage());
            return new BigDecimalMatrix(systemSize, systemSize, 12);
        }
    }

    public BigDecimalMatrix funcP(int i, BigDecimal t) {
        BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, systemSize, 12);
        BigDecimalMatrix val;
        int j = i;
        int hash = Objects.hash(i, t);
        boolean vallidness;
        try {
            if (MatrixContainer.getPMatrix().containsKey(hash)) {
                return MatrixContainer.getPMatrix().get(hash).clone();
            }
            do {
                val = funcK(i, j);

                double valueUnderSum = Math.exp(tetta.multiply(t).doubleValue() * -1)
                        * Math.pow(tetta.multiply(t).doubleValue(), j)
                        / GeneratorCreator.factor(j);
                val = val.multiply(valueOf(valueUnderSum));
                sum = sum.add(val);
                j++;

                vallidness = val.squaredEuclidianNorm().doubleValue() >= accuracy.doubleValue();
            } while (vallidness);
            MatrixContainer.getPMatrix().put(hash, sum.clone());
        } catch (CloneNotSupportedException e) {
            System.out.println(e.getMessage());
        }

        return sum;
    }

    public static BigDecimal funcSmallPhiK(BigDecimal gamma, BigDecimal t, int k) {
        return new BigDecimal(gamma.multiply(t).pow(k).doubleValue() / factor(k)
                * Math.exp(-gamma.multiply(t).doubleValue()));
    }

    private BigDecimal funcSmallPhiKWithHat(int k) {
        BigDecimal result = ZERO;
        BigDecimal sum;
        int i = k;
        do {
            sum = GeneratorCreator.funcSmallPhiK(gamma, T, i);
            result = result.add(sum);
            i++;
        } while (sum.doubleValue() >= accuracy.doubleValue());
        return result;
    }

    public BigDecimalMatrix funcPhi(int i, int k) {
        return funcP(i, T).multiply(GeneratorCreator.funcSmallPhiK(gamma, T, k));
    }

    private BigDecimalMatrix funcPhiWithHat(int i, int k) {
        return funcP(i, T).multiply(funcSmallPhiKWithHat(k));
    }

    public BigDecimalMatrix funcN(int n) {
        BigDecimalMatrix result = new BigDecimalMatrix(systemSize, systemSize, 12);
        BigDecimalMatrix sum;
        int j = n;
        do {
            sum = funcK(n, j).multiply(tetta.pow(j).divide((gamma.add(tetta)).pow(j + 1), 12, RoundingMode.CEILING));
            result = result.add(sum);
            j++;
        } while (sum.squaredEuclidianNorm().compareTo(accuracy) > 0);
        return result.multiply(gamma);
    }

    private BigDecimal factor(BigDecimal n) {
        if (n.signum() < 0) {
            throw new IllegalArgumentException("Negative value in factorial function.");
        }
        if (n.scale() < 0) {
            throw new IllegalArgumentException("Real value in factorial function.");
        }
        if (n.equals(ZERO)) {
            return ONE;
        }
        return n.multiply(factor(n.subtract(ONE)));
    }

    public static int factor(int n) {
        int result = 1;
        for (int i = 1; i <= n; i++) {
            result = result * i;
        }
        return result;
    }

    private BigDecimal tetta() {
        BigDecimal tetta = ZERO;
        for (int i = 0; i < d0.getHeight(); i++) {
            BigDecimal val = d0.getElement(i, i).multiply(new BigDecimal(-1));
            tetta = val.compareTo(tetta) == 1 ? val : tetta;
        }
        return tetta;
    }

    public BigDecimal getAccuracy() {
        return new BigDecimal(accuracy.toString());
    }

    public int getK() {
        return K;
    }
}
