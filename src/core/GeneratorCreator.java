package core;

import kurs.BigDecimalMatrix;

import java.math.BigDecimal;
import java.util.Objects;

import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;
import static java.math.BigDecimal.valueOf;


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

    public BigDecimalMatrix create(int i, int j) {
        BigDecimalMatrix tempMatrix, matrix = new BigDecimalMatrix(2 * K + 2, 2 * K + 2, 6);

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
                        matrix.setElement(k * 2, bK * 2, tempMatrix.getElement(0, 0));
                        matrix.setElement(k * 2, 1 + bK * 2, tempMatrix.getElement(0, 1));
                        matrix.setElement(1 + k * 2, bK * 2, tempMatrix.getElement(1, 0));
                        matrix.setElement(1 + k * 2, 1 + bK * 2, tempMatrix.getElement(1, 1));
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

            if (bK == K - 1) {
                if (isVMatrix) {
                    element = createFirstRowPenultimateElementV(j, bK);
                    return element;
                }
            }

            if (bK == K) {
                if (isVMatrix) {
                    element = k == 0 ? createFirstRowLastColumnElementV(j) : createLastColumnElementV(k, j);
                    return element;
                }
                element = createLastColumnElement(k, j);
                return element;
            }

            if (k - bK > 1) {
                element = BigDecimalMatrix.zeroMatrix(2);
                return element;
            }

            if (k == 0) {
                element = isVMatrix ? createFirstRowNotLastColumnElementV(j, bK) : createFirstRowNotLastColumnElement(j, bK);
                return element;
            }

            element = isVMatrix ? createNotFirstRowNotLastColumnElementV(k, j, bK) : createNotFirstRowNotLastColumnElement(k, j, bK);
            return element;
        }

        // YMatrix region
        private BigDecimalMatrix createFirstRowNotLastColumnElement(int j, int bK) {
            BigDecimalMatrix result = new BigDecimalMatrix(2, 2, 10);
            for (int n = 0; n <= j; n++) {
                result = result.add(funcN(n).multiply(funcPhi(j - n, bK)));
            }
            return result;
        }

        private BigDecimalMatrix createLastColumnElement(int k, int j) {
            return funcPhiWithHat(j, K - k + 1);
        }

        private BigDecimalMatrix createNotFirstRowNotLastColumnElement(int k, int j, int bK) {
            return funcPhi(j, bK - k + 1);
        }
        // end Region

        //VMatrix region

        private BigDecimalMatrix createFirstRowLastColumnElementV(int j) {
            BigDecimalMatrix result = new BigDecimalMatrix(2, 2, 6);
            for (int n = 0; n <= j; n++) {
                result = result.add(funcN(n).multiply(funcPhiWithHat(j - n, K)));
            }
            result = result.multiply(funcM(0));

            BigDecimalMatrix sum = new BigDecimalMatrix(2, 2, 6);
            for (int m = 0; m <= K - 1; m++) {
                sum = sum.add(funcM(m).multiply(funcPhiWithHat(j, K - m)));
            }
            sum = sum.add(funcMWithHat(K).multiply(funcPhiWithHat(j, 1)));
            sum = sum.multiply(funcN(0));

            result = result.add(sum);
            return result;
        }


        private BigDecimalMatrix createFirstRowPenultimateElementV(int j, int bK) {
            BigDecimalMatrix result = new BigDecimalMatrix(2, 2, 6);
            for (int m = 0; m <= K - bK; m++) {
                result = result.add(funcM(m).multiply(funcPhi(j, K - bK - m)));
            }
            result = result.add(funcMWithHat(K - bK + 1).multiply(funcPhi(j, 0)));
            return result;
        }

        private BigDecimalMatrix createFirstRowNotLastColumnElementV(int j, int bK) {
            BigDecimalMatrix result = new BigDecimalMatrix(2, 2, 6);

            for (int n = 0; n <= j; n++) {
                result = result.add(funcN(n).multiply(funcPhi(j - n, bK)));
            }

            result = funcM(0).multiply(result);

            BigDecimalMatrix sum = new BigDecimalMatrix(2, 2, 6);

            for (int m = 0; m <= bK; m++) {
                sum = sum.add(funcM(m).multiply(funcPhi(j, bK - m)));
            }

            sum = funcN(0).multiply(sum);

            result = result.add(sum);

            return result;
        }

        private BigDecimalMatrix createLastColumnElementV(int k, int j) {
            BigDecimalMatrix sum = new BigDecimalMatrix(2, 2, 6);
            for (int m = 0; m <= K - k + 1; m++) {
                sum = sum.add(funcM(m).multiply(funcPhiWithHat(j, K - k + 1 - m)));
            }
            sum = sum.add(funcMWithHat(K - k + 1).multiply(funcPhiWithHat(j, 1)));
            return sum;
        }

        private BigDecimalMatrix createNotFirstRowNotLastColumnElementV(int k, int j, int bK) {
            BigDecimalMatrix sum = new BigDecimalMatrix(2, 2, 6);
            for (int m = 0; m <= bK - k + 1; m++) {
                sum = sum.add(funcM(m).multiply(funcPhi(j, bK - k + 1 - m)));
            }
            return sum;
        }
        // End region
    }

    private BigDecimalMatrix funcMWithHat(int m) {
        int n = m;
        BigDecimalMatrix mWithHat = new BigDecimalMatrix(2, 2, 10);
        BigDecimalMatrix mn;
        do {
            mn = funcM(n);
            mWithHat = mWithHat.add(mn);
            n++;
        } while (mn.squaredEuclidianNorm().doubleValue() >= accuracy.doubleValue());
        return mWithHat;
    }

    public BigDecimalMatrix funcM(int n) {
        BigDecimalMatrix val;
        BigDecimalMatrix negativeD0 = d0.multiply(new BigDecimal(-1));
        BigDecimal gammaN = gamma.pow(n);
        val = (negativeD0.add(BigDecimalMatrix.identity(2).multiply(gamma))).inverse();
        val = BigDecimalMatrix.powMatrix(val, n + 1);
        val = val.multiply(gammaN).multiply(d1);
        return val;
    }

    private BigDecimalMatrix funcK(int n, int j) {
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

    private BigDecimalMatrix funcP(int i, BigDecimal t) {
        BigDecimalMatrix sum = new BigDecimalMatrix(2, 2, 6);
        BigDecimalMatrix val;
        int j = 0;
        int hash = Objects.hash(i, t);
        boolean vallidness;
        try {
            if (MatrixContainer.getPMatrix().containsKey(hash)) {
                return MatrixContainer.getPMatrix().get(hash).clone();
            }
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

               // System.out.println(val);
                vallidness = (j == 0) || (j == 1) || (j == 2) || (j == 3) || val.cubicNorm().doubleValue() >= accuracy.doubleValue();
                /*vallidness = (val.cubicNorm().doubleValue() == 0 && valueUnderSum >= accuracy.doubleValue())
                        || val.cubicNorm().doubleValue() >= accuracy.doubleValue();*/
                //    } while (j <= 15/*val.cubicNorm().doubleValue() >= accuracy.doubleValue()*/);
            } while (vallidness);
            MatrixContainer.getPMatrix().put(hash, sum.clone());
        } catch (CloneNotSupportedException e) {
            System.out.println(e.getMessage());
        }

        return sum;
    }

    public BigDecimal funcSmallPhiK(BigDecimal t, int k) {
        return new BigDecimal(gamma.multiply(t).pow(k).doubleValue() / factor(k)
                * Math.exp(-gamma.multiply(t).doubleValue()));
    }

    private BigDecimal funcSmallPhiKWithHat(BigDecimal t, int k) {
        BigDecimal result = ZERO;
        BigDecimal sum;
        int i = k;
        do {
            sum = funcSmallPhiK(t, i);
            result = result.add(sum);
            i++;
        } while (sum.doubleValue() >= accuracy.doubleValue());
        return result;
    }

    private BigDecimalMatrix funcPhi(int i, int k) {
        return funcP(i, T).multiply(funcSmallPhiK(T, k));
    }

    private BigDecimalMatrix funcPhiWithHat(int i, int k) {
        return funcP(i, T).multiply(funcSmallPhiKWithHat(T, k));
    }

    public BigDecimalMatrix funcN(int n) {
        BigDecimalMatrix result = new BigDecimalMatrix(2, 2, 6);
        BigDecimalMatrix sum;
        boolean vallidness;
        int j = 0;
        do {
            sum = funcK(n, j).multiply(gamma)
                    .multiply(new BigDecimal(tetta().pow(j).doubleValue() / (gamma.add(tetta()).pow(j + 1).doubleValue())));
            result = result.add(sum);
            vallidness = (j == 0) || (j == 1) || (j == 2) || (j == 3) || sum.cubicNorm().doubleValue() >= accuracy.doubleValue();
            j++;
        } while (vallidness);
        return result;
    }

    public BigDecimalMatrix checkP() {
        BigDecimalMatrix sum = new BigDecimalMatrix(2, 2, 6);
        BigDecimalMatrix val;
        int i = 0;
        do {
            val = funcP(i, T);
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
              //  System.out.println(val);
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
