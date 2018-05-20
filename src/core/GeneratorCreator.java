package core;

import kurs.BigDecimalMatrix;

import java.math.BigDecimal;
import java.math.RoundingMode;

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

        if(!ErgodicityCondition.check(gamma, lambda, K, T)) {
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

        for (int r = 0; r < K + 1; r++) {
            for (int c = 0; c < K + 1; c++) {
                matrix[r][c] =  elementCreator.create(i, r, j, c);
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
            if (k - bK > 1) {
                element = new BigDecimalMatrix(2, 2, 10);
            } else if (k == 0) {
                if (bK == K) {
                    element = createLastColumnElement(i, k, j, bK);
                } else {
                    element = createFirstRowNotLastColumnElement(i, k, j, bK);
                }
            } else {
                if (bK == K) {
                    element = createLastColumnElement(i, k, j, bK);
                } else {
                    element = createNotFirstRowNotLastColumnElement(i, k, j, bK);
                }
            }
            //System.out.println(element);
            return element;
        }
        private BigDecimalMatrix createFirstRowNotLastColumnElement(int i, int k, int j, int bK) {
            BigDecimalMatrix result = new BigDecimalMatrix(2, 2, 10);
            for(int n = 0; n <= j - 2; n++) {
                result = result.add(
                        funcN(n).multiply(funcPhi(new BigDecimal(j - 2 - n), new BigDecimal(bK)))
                );
            }
            return result;
        }

        //Approved
        private BigDecimalMatrix createLastColumnElement(int i, int k, int j, int bK) {
            return funcPhiWithHat(new BigDecimal(j), new BigDecimal(bK - k + 1));
        }

        //Approved
        private BigDecimalMatrix createNotFirstRowNotLastColumnElement(int i, int k, int j, int bK) {
            return funcPhi(new BigDecimal(j), new BigDecimal(bK - k + 1));
        }
    }

    public BigDecimalMatrix funcMWithHat(int m) {
        int n = m;
        BigDecimalMatrix mWithHat = new BigDecimalMatrix(2, 2, 10);
        BigDecimalMatrix mn = funcM(n);
        while(Math.sqrt(mn.squaredEuclidianNorm().doubleValue()) >= accuracy.doubleValue()) {
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
        for(int i = 0; i < n; i++) {
            val = val.multiply(sumInverse);
        }
        return val.multiply(gamma.pow(n)).multiply(d1);
    }

    public BigDecimalMatrix funcK(int n, int j) {
        if (j == 0) {
            return n == 0 ? new BigDecimalMatrix(2, new BigDecimal(1), 10) : new BigDecimalMatrix(2, 2, 10);
        } else {
            if (n == 0) {
                BigDecimalMatrix sum = new BigDecimalMatrix(2, new BigDecimal(1), 10)
                        .add(d0.multiply(new BigDecimal(1)
                                .divide(tetta(), RoundingMode.HALF_UP)));
                sum = funcK(0, j - 1).multiply(sum);
                return sum;
            } else {
                BigDecimalMatrix sum = new BigDecimalMatrix(2, 2, 10);
                for (int i = 0; i <= n - 1; i++) {
                    if( n - i == 1) {
                        sum = sum.add(funcK(i, j - 1).multiply(d1));
                    }
                    sum = sum.add(funcK(n, j - 1)
                                    .multiply(new BigDecimalMatrix(2, new BigDecimal(1), 10)
                                            .add(d0.multiply(new BigDecimal(1).divide(tetta(), RoundingMode.HALF_UP)))));
                }
                return sum.multiply(new BigDecimal(1).divide(tetta(), RoundingMode.HALF_UP));
            }
        }
    }

    public BigDecimalMatrix funcP(int n, BigDecimal t) {
        BigDecimalMatrix sum = new BigDecimalMatrix(2, 2, 10);
        BigDecimalMatrix val;
        int j = 1;
        do {
            //val = funcK(j, n);
            val = funcK(n, j);
            val = val.multiply(new BigDecimal(Math.exp(tetta().multiply(new BigDecimal(-1)).multiply(t).doubleValue())))
                    .multiply((tetta().multiply(t)).pow(j)
                            .divide(factor(new BigDecimal(j)), RoundingMode.HALF_UP));
            sum = sum.add(val);
            j++;
        } while (Math.sqrt(val.squaredEuclidianNorm().doubleValue()) >= accuracy.doubleValue());
        return sum;
    }

    public BigDecimal funcSmalPhiK(BigDecimal t, int k) {
        return (gamma.multiply(t).pow(k))
                .divide(factor(new BigDecimal(k)), RoundingMode.HALF_UP)
                .multiply(new BigDecimal(Math.exp(gamma.multiply(t).doubleValue())));
    }

    public BigDecimal funcSmalPhiKWithHat(BigDecimal t, int k) {
        BigDecimal val = funcSmalPhiK(t, k);
        BigDecimal result = ZERO;
        int i = k;
        while(val.compareTo(accuracy) == 1 ) {
            result = result.add(val);
            i++;
            val = funcSmalPhiK(t, i);
        }
        return result;
    }

    public BigDecimalMatrix funcPhi(BigDecimal i, BigDecimal k) {
        return funcP(i.intValue(), T).multiply(funcSmalPhiK(T, k.intValue()));
    }

    public BigDecimalMatrix funcPhiWithHat(BigDecimal i, BigDecimal k) {
        return funcP(i.intValue(), T).multiply(funcSmalPhiKWithHat(T, k.intValue()));
    }

    public BigDecimalMatrix funcN(int n) {
        return funcP(n, T).multiply(gamma
                .multiply(new BigDecimal(Math.exp(gamma.multiply(new BigDecimal(-1)).multiply(T).doubleValue()))));
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

    private BigDecimal tetta() {
        BigDecimal tetta = ZERO;
        for (int i = 0; i < d0.getHeight(); i++) {
            tetta = (d0.getElement(i, i).multiply(new BigDecimal(-1)).compareTo(tetta) == 1) ? d0.getElement(i, i).multiply(new BigDecimal(-1)) : tetta;
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
