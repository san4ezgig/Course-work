package core;

import kurs.BigDecimalMatrix;

import java.math.BigDecimal;
import java.util.Objects;

import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;


public class GeneratorCreator {
    private BigDecimal gamma;
    private BigDecimal lambda;
    private Integer K;
    private BigDecimal accuracy;
    private BigDecimalMatrix d0;
    private BigDecimalMatrix d1;
    private BigDecimal T;
    private MatrixElementCreator elementCreator;

    public GeneratorCreator(
            BigDecimal gamma,
            BigDecimal lambda,
            int K,
            BigDecimalMatrix d0,
            BigDecimalMatrix d1,
            BigDecimal T,
            BigDecimal accuracy
    ) {
        boolean validness = gamma.compareTo(ONE) < 0 &&
                gamma.compareTo(ZERO) > 0 &&
                lambda.compareTo(ONE) < 0 &&
                lambda.compareTo(ZERO) > 0 &&
                K > 0 &&
                gamma.compareTo(lambda) > 0;

        if (!ErgodicityCondition.check(gamma.doubleValue(), lambda, K, T)) {
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
        this.T = T;
        this.d1 = d1;
        this.elementCreator = new MatrixElementCreator();
    }

    public BigDecimalMatrix[][] create(int i, int j) {
        BigDecimalMatrix[][] matrix = new BigDecimalMatrix[K + 1][K + 1];

        if (i - j > 1) {
            return matrix;
        }

        for (int k = 0; k < K + 1; k++) {
            for (int bK = 0; bK < K + 1; bK++) {
                matrix[k][bK] = elementCreator.create(i, k, j, bK);
                //System.out.println(matrix[r][c]);
            }
        }

        return matrix;
    }

    private class MatrixElementCreator {

        //Approved
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

            if (bK == K) {
                element = createLastColumnElement(k, j);
                return element;
            }

            if (k - bK > 1) {
                element = new BigDecimalMatrix(2, 2, 10);
                return element;
            }

            if (k == 0) {
                element = createFirstRowNotLastColumnElement(j, bK);
                return element;
            }

            element = createNotFirstRowNotLastColumnElement(k, j, bK);
            return element;
        }

        private BigDecimalMatrix createFirstRowNotLastColumnElement(int j, int bK) {
            BigDecimalMatrix result = new BigDecimalMatrix(2, 2, 10);
            for (int n = 0; n <= j - 2; n++) {
                result = result.add(
                        funcN(n).multiply(funcPhi(j - 2 - n, bK))
                );
            }
            return result;
        }

        //Approved
        private BigDecimalMatrix createLastColumnElement(int k, int j) {
            return funcPhiWithHat(new BigDecimal(j), new BigDecimal(K - k + 1));
        }

        //Approved
        private BigDecimalMatrix createNotFirstRowNotLastColumnElement(int k, int j, int bK) {
            return funcPhi(j, bK - k + 1);
        }
    }

    public BigDecimalMatrix funcMWithHat(int m) {
        int n = m;
        BigDecimalMatrix mWithHat = new BigDecimalMatrix(2, 2, 10);
        BigDecimalMatrix mn = funcM(n);
        while (Math.sqrt(mn.squaredEuclidianNorm().doubleValue()) >= accuracy.doubleValue()) {
            mWithHat = mWithHat.add(mn);
            n++;
            mn = funcM(n);
        }
        return mWithHat;
    }

    public BigDecimalMatrix funcM(int n) {
        BigDecimalMatrix minusD0 = d0.multiply(new BigDecimal(-1));
        BigDecimalMatrix gammaI = new BigDecimalMatrix(2, new BigDecimal(1), 10).multiply(gamma);
        BigDecimalMatrix sumInverse = minusD0.add(gammaI).inverse();
        BigDecimalMatrix val = sumInverse;
        for (int i = 0; i < n; i++) {
            val = val.multiply(sumInverse);
        }
        return val.multiply(gamma.pow(n)).multiply(d1);
    }

    public BigDecimalMatrix funcK(int n, int j) {
        try {
            int hash = Objects.hash(n, j);

            if (MatrixContainer.getKMatrix().containsKey(hash)) {
                return MatrixContainer.getKMatrix().get(hash).clone();
            }

            BigDecimal negTetta = new BigDecimal(1 / tetta().doubleValue());
            if (j == 0) {
                return n == 0 ? new BigDecimalMatrix(2, new BigDecimal(1), 6) : new BigDecimalMatrix(2, 2, 6);
            } else {
                if (n == 0) {
                    BigDecimalMatrix sum = new BigDecimalMatrix(2, new BigDecimal(1), 6).add(d0.multiply(negTetta));
                    sum = funcK(0, j - 1).multiply(sum);
                    MatrixContainer.getKMatrix().put(hash, sum.clone());
                    return sum;
                } else {
                    BigDecimalMatrix sum = new BigDecimalMatrix(2, 2, 6);

                    if (n >= 1) {
                        sum = sum.add(funcK(n - 1, j - 1).multiply(d1).multiply(negTetta));
                    }

                    sum = sum.add(funcK(n, j - 1)
                            .multiply(new BigDecimalMatrix(2, new BigDecimal(1), 6).add(d0.multiply(negTetta))));
                    MatrixContainer.getKMatrix().put(hash, sum.clone());
                    return sum;
                }
            }
        } catch (CloneNotSupportedException e) {
            System.out.println(e.getMessage());
            return new BigDecimalMatrix(2, 2, 6);
        }
    }

    public BigDecimalMatrix funcP(int i, BigDecimal t) {
        BigDecimalMatrix sum = new BigDecimalMatrix(2, 2, 6);
        BigDecimalMatrix val;
        int j = 0;
        do {
            if (i == 0) {
                val = new BigDecimalMatrix(2, new BigDecimal(1), 6)
                        .add(d0.multiply(new BigDecimal(1 / tetta().doubleValue())));
                BigDecimalMatrix newVal = val;
                if (j == 0) {
                    val = new BigDecimalMatrix(2, new BigDecimal(1), 6);
                }

                for (int k = 0; k < j; k++) {
                    val = val.multiply(newVal);
                }
            } else {
                val = funcK(i, j);
            }

            double valueUnderSum = Math.exp(tetta().multiply(t).doubleValue() * -1)
                    * Math.pow(tetta().multiply(t).doubleValue(), j)
                    / this.factor(j);
            val = val.multiply(new BigDecimal(valueUnderSum));
            sum = sum.add(val);
            j++;
        } while (j <= 10/*val.cubicNorm().doubleValue() >= accuracy*/);
        return sum;
    }

    public BigDecimal funcSmallPhiK(BigDecimal t, int k) {
        return new BigDecimal(gamma.multiply(t).pow(k).doubleValue() / factor(k)
                * Math.exp(gamma.multiply(t).doubleValue() * -1));
    }

    public BigDecimal funcSmalPhiKWithHat(BigDecimal t, int k) {
        BigDecimal val = funcSmallPhiK(t, k);
        BigDecimal result = ZERO;
        int i = k;
        while (val.compareTo(accuracy) == 1) {
            result = result.add(val);
            i++;
            val = funcSmallPhiK(t, i);
        }
        return result;
    }

    public BigDecimalMatrix funcPhi(int i, int k) {
        return funcP(i, T).multiply(funcSmallPhiK(T, k));
    }

    public BigDecimalMatrix funcPhiWithHat(BigDecimal i, BigDecimal k) {
        return funcP(i.intValue(), T).multiply(funcSmalPhiKWithHat(T, k.intValue()));
    }

    public BigDecimalMatrix funcN(int n) {
        return funcP(n, T).multiply(gamma
                .multiply(new BigDecimal(Math.exp(gamma.multiply(new BigDecimal(-1)).multiply(T).doubleValue()))));
    }

    public BigDecimalMatrix checkP() {
        BigDecimalMatrix sum = new BigDecimalMatrix(2, 2, 6);
        BigDecimalMatrix val;
        int i = 0;
        do {
            val = funcP(i, T);
            System.out.println(val);
            sum = sum.add(val);
            i++;
        } while (i <= 5);
        sum.setScale(6);
        return sum;
    }

    public BigDecimalMatrix checkPhi() {
        BigDecimalMatrix sum = BigDecimalMatrix.eCol(2, ZERO);
        BigDecimalMatrix e = BigDecimalMatrix.eCol(2, ONE);
        BigDecimalMatrix val;
        int i = 0;
        do {
            val = BigDecimalMatrix.eCol(2, ZERO);
            for (int k = 0; k <= K; k++) {
                val = val.add(funcPhi(i, k).multiply(e));
                System.out.println(val);
            }
            sum = sum.add(val);
            i++;
        } while (i <= 5);
        sum.setScale(6);
        return sum;
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

    private int factor(int n) {
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
