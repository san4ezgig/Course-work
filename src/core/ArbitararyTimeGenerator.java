package core;

import kurs.BigDecimalMatrix;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Objects;

import static java.math.BigDecimal.*;

public class ArbitararyTimeGenerator {
    private GeneratorCreator creator;
    private BigDecimal gamma;
    private BigDecimal lambda;
    private Integer K;
    private BigDecimal accuracy;
    private BigDecimalMatrix d0;
    private BigDecimalMatrix d1;
    private BigDecimal T;
    private BigDecimal tetta;
    private int systemSize;
    private ArbitararyTimeGenerator.MatrixElementCreator elementCreator;
    private ArrayList<BigDecimalMatrix> piVector;

    public ArbitararyTimeGenerator(GeneratorCreator creator, ArrayList<BigDecimalMatrix> piVector) {
        this.creator = creator;
        this.gamma = creator.gamma;
        this.lambda = creator.lambda;
        this.K = creator.K;
        this.accuracy = creator.accuracy;
        this.d0 = creator.d0;
        this.tetta = creator.tetta;
        this.T = creator.T;
        this.d1 = creator.d1;
        this.systemSize = creator.systemSize;
        this.piVector = piVector;
        this.elementCreator = new MatrixElementCreator();
    }

    public BigDecimalMatrix create(int i, int j) {
        BigDecimalMatrix tempMatrix, resultMatrix = new BigDecimalMatrix(systemSize * (K + 1), systemSize * (K + 1), 12);

        int hash = i > 0 ? Objects.hash(i, j - i) : Objects.hash(0, j);

//        int hash = Objects.hash(i, j);

        try {
            if (MatrixContainer.getArbitryTimeGenerators().containsKey(hash)) {
                resultMatrix = MatrixContainer.getArbitryTimeGenerators().get(hash).clone();
            } else {
                for (int k = 0; k < K + 1; k++) {
                    for (int bK = 0; bK < K + 1; bK++) {
                        tempMatrix = elementCreator.create(i, k, j, bK);
                        switch (systemSize) {
                            case 1: {
                                resultMatrix.setElement(k, bK, tempMatrix.getElement(0, 0));
                                break;
                            }
                            case 2: {
                                resultMatrix.setElement(k * systemSize, bK * systemSize, tempMatrix.getElement(0, 0));
                                resultMatrix.setElement(k * systemSize, 1 + bK * systemSize, tempMatrix.getElement(0, 1));
                                resultMatrix.setElement(1 + k * systemSize, bK * systemSize, tempMatrix.getElement(1, 0));
                                resultMatrix.setElement(1 + k * systemSize, 1 + bK * systemSize, tempMatrix.getElement(1, 1));
                            }
                        }
                    }
                }
                MatrixContainer.getArbitryTimeGenerators().put(hash, resultMatrix.clone());
            }
        } catch (CloneNotSupportedException e) {
            System.out.println(e.getMessage());
        }
        return resultMatrix;
    }

    private class MatrixElementCreator {

        private BigDecimalMatrix create(int i, int k, int j, int bK) {
            BigDecimalMatrix zeroMatrix = new BigDecimalMatrix(systemSize, ZERO, 15);

            boolean isVMatrix = (i == 0);

            if (isVMatrix) {
                if (j == 0) {

                    if (bK == K) {
                        return createLastColumnElementZeroV(k);
                    }
                    if (k > bK) {
                        return zeroMatrix;
                    }
                    return createNotLastColumnElementZeroV(k, bK);
                }
                else {
                    if (k - bK > 1) {
                        return zeroMatrix;
                    }
                    if (k == 0) {
                        if (bK == K - 1) {
                            return createFirstRowPenultimateElementV(j);
                        }

                        if (bK == K) {
                            return createFirstRowLastColumnElementV(j);
                        }

                        return createFirstRowNotLastColumnElementV(j, bK);
                    }
                    else {
                        if (bK == K - 1) {
                            return createNotFirstRowPenultimateElementV(j , k);
                        }
                        if (bK == K) {
                            return createNotFirstRowLastColumnElementV(k, j);
                        }

                        return createNotFirstRowNotLastColumnElementV(k, j, bK);
                    }
                }
            }

            if (k - bK > 1) {
                return zeroMatrix;
            }

            if (bK == K) {
                return k == 0 ? createFirstRowLastColumnElement(j, i) : createNotFirstRowLastColumnElement(k, j, i);
            }

            if (k == 0) {
                return createFirstRowNotLastColumnElement(j, bK, i);
            }

            return createNotFirstRowNotLastColumnElement(k, j, bK, i);
        }

        // VMatrix region
        // VMatrix j == 0
        private BigDecimalMatrix createNotLastColumnElementZeroV(int k, int bK) {
            BigDecimalMatrix result = BigDecimalMatrix.identity(systemSize).multiply(gamma).subtract(d0).inverse();
            result = BigDecimalMatrix.powMatrix(result, bK - k + 1);
            result = result.multiply(valueOf(Math.pow(gamma.doubleValue(), bK - k)));
            return result;
        }

        private BigDecimalMatrix createLastColumnElementZeroV(int k) {
            BigDecimalMatrix result = d0.multiply(ONE.negate()).inverse();
            BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, ZERO, 10);
            BigDecimalMatrix val;
            for (int l = 0; l <= K - k -1; l++) {
                val = BigDecimalMatrix.powMatrix(BigDecimalMatrix.identity(systemSize).multiply(gamma).subtract(d0).inverse(), l + 1);
                val = val.multiply(valueOf(Math.pow(gamma.doubleValue(), l)));
                sum = sum.add(val);
            }
            result = result.subtract(sum);
            return result;
        }
        //end

        private BigDecimalMatrix createFirstRowLastColumnElementV(int j) {
            BigDecimalMatrix result;
            BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, ZERO, 10);

            for (int l = 0; l <= K - 1; l++) {
                sum = sum.add(creator.funcM(l).multiply(funcPhiWithWaveWithHat(j - 1, K - l)));
            }
            sum = sum.add(creator.funcMWithHat(K).multiply(funcPhiWithWaveWithHat(j - 1, 1)));

            result = creator.funcN(0).multiply(sum);


            sum = BigDecimalMatrix.zeroMatrix(systemSize);
            for (int l = 0; l <= j - 1; l++) {
                sum = sum.add(
                        creator.funcN(l).multiply(
                                funcPhiWithWaveWithHat(j - l - 1, K)
                        )
                );
            }
            result = result.add(creator.funcM(0).multiply(sum));
            return result;
        }

        private BigDecimalMatrix createFirstRowPenultimateElementV(int j) {
            BigDecimalMatrix result = new BigDecimalMatrix(systemSize, ZERO, 10);
            BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, ZERO, 10);
            int cronecer = (K == 1) ? 1 : 0;
            result = creator.funcN(j).multiply(valueOf(1 / gamma.doubleValue())).multiply(valueOf(cronecer));

            for (int l = 0; l <= j - 1; l++) {
                sum = sum.add(
                        creator.funcN(l).multiply(funcPhiWithWave( j - l - 1, K - 1))
                );
            }
            sum = creator.funcM(0).multiply(sum);
            result = result.add(sum);

            sum = new BigDecimalMatrix(systemSize, systemSize, 12);

            for (int l = 0; l <= K - 1; l++) {
                sum = sum.add(
                        creator.funcM(l).multiply(funcPhiWithWave(j - 1, K - 1 - l))
                );
            }
            sum = sum.add(creator.funcMWithHat(K).multiply(funcPhiWithWave(j - 1, 0)));
            sum = creator.funcN(0).multiply(sum);
            result = result.add(sum);

            return result;
        }

        private BigDecimalMatrix createFirstRowNotLastColumnElementV(int j, int bK) {
            BigDecimalMatrix result;
            BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, ZERO, 10);
            int cronecer = (bK == 0) ? 1 : 0;
            result = creator.funcN(j).multiply(valueOf(1 / gamma.doubleValue())).multiply(valueOf(cronecer));

            for (int l = 0; l <= bK; l++) {
                sum = sum.add(
                        creator.funcM(l).multiply(funcPhiWithWave(j - 1, bK - l))
                );
            }
            sum = creator.funcN(0).multiply(sum);

            result = result.add(sum);

            sum = new BigDecimalMatrix(systemSize, systemSize, 12);
            for (int l = 0; l <= j - 1; l++) {
                sum = sum.add(
                        creator.funcN(l).multiply(funcPhiWithWave(j - l - 1, bK))
                );
            }
            sum = creator.funcM(0).multiply(sum);

            result = result.add(sum);

            return result;
        }

        private BigDecimalMatrix createNotFirstRowPenultimateElementV(int j, int k) {
            BigDecimalMatrix result = new BigDecimalMatrix(systemSize, systemSize, 12);

            for (int l = 0; l <= K - k; l++) {
                result = result.add(
                        creator.funcM(l).multiply(funcPhiWithWave(j - 1, K - k - l))
                );
            }

            result = result.add(
                    creator.funcMWithHat(K - k + 1).multiply(funcPhiWithWave(j - 1, 0))
            );

            return result;
        }

        private BigDecimalMatrix createNotFirstRowLastColumnElementV(int k, int j) {
            BigDecimalMatrix result = new BigDecimalMatrix(systemSize, systemSize, 12);

            for (int l = 0; l <= K - k; l++) {
                result = result.add(creator.funcM(l).multiply(funcPhiWithWaveWithHat(j - 1, K - k - l + 1)));
            }

            result = result.add(
                    creator.funcMWithHat(K - k + 1).multiply(funcPhiWithWaveWithHat(j - 1, 1))
            );

            return result;
        }

        private BigDecimalMatrix createNotFirstRowNotLastColumnElementV(int k, int j, int bK) {
            BigDecimalMatrix sum = new BigDecimalMatrix(systemSize, systemSize, 12);
            for (int l = 0; l <= bK - k + 1; l++) {
                sum = sum.add(creator.funcM(l).multiply(funcPhiWithWave(j - 1, bK - k - l + 1)));
            }
            return sum;
        }

        // end region

        // YMatrixRegion

        private BigDecimalMatrix createFirstRowNotLastColumnElement(int j, int bK, int i) {
            int r = j - i;
            BigDecimalMatrix result = new BigDecimalMatrix(systemSize, systemSize, 10);
            for (int l = 0; l <= r; l++) {
                result = result.add(creator.funcN(l).multiply(funcPhiWithWave(r - l, bK)));
            }
            return result;
        }

        private BigDecimalMatrix createFirstRowLastColumnElement(int j, int i) {
            BigDecimalMatrix result = new BigDecimalMatrix(systemSize, systemSize, 10);
            int r = j - i;
            for (int l = 0; l <= r; l++) {
                result = result.add(creator.funcN(l).multiply(funcPhiWithWaveWithHat(r - l, K)));
            }
            return result;
        }

        private BigDecimalMatrix createNotFirstRowLastColumnElement(int k, int j, int i) {
            return funcPhiWithWaveWithHat(j - i, K - k + 1);
        }

        private BigDecimalMatrix createNotFirstRowNotLastColumnElement(int k, int j, int bK, int i) {
            return funcPhiWithWave(j - i, bK - k + 1);
        }

        // end region
    }

    public BigDecimalMatrix calculateP(int j) {
        BigDecimalMatrix result;
        BigDecimalMatrix sum = BigDecimalMatrix.eRow(systemSize * (K + 1), ZERO);
        result = piVector.get(0).multiply(create(0, j));

        for (int r = 0; r <= j - 1; r++) {
            sum = sum.add(piVector.get(j - r).multiply(create(1, r + 1)));
        }
        result = result.add(sum).multiply(lambda);

        return result;
    }

    public BigDecimalMatrix funcPhiWithWave(int n, int k) {
        BigDecimalMatrix result = new BigDecimalMatrix(systemSize, ZERO, 12);
        BigDecimal sum;
        BigDecimal sumTettaGamma = tetta.add(gamma);
        BigDecimalMatrix val;
        int j = n;
        do {
            sum = ZERO;
            val = creator.funcK(n, j);
            double doubleValue = GeneratorCreator.factor(valueOf(j + k)).doubleValue() / GeneratorCreator.factor(valueOf(j)).doubleValue();
            doubleValue = doubleValue == Double.POSITIVE_INFINITY  ? 0 : doubleValue;
            val = val.multiply(
                    valueOf(
                            doubleValue
                    )
            );
            val = val.multiply(
                    valueOf(
                            Math.pow(tetta.doubleValue() / sumTettaGamma.doubleValue(), j)
                    )
            );

            for (int l = 0; l <= j + k; l++) {
                doubleValue = Math.pow(T.multiply(sumTettaGamma).doubleValue(), l)
                        / GeneratorCreator.factor(valueOf(l)).doubleValue();
                doubleValue = doubleValue == Double.POSITIVE_INFINITY  ? 0 : doubleValue;
                sum = sum.add(
                        valueOf(
                                doubleValue
                        )
                );
            }
            sum = sum.multiply(
                    valueOf(
                            Math.exp(T.multiply(sumTettaGamma).negate().doubleValue())
                    )
            );
            val = val.multiply(ONE.subtract(sum));

            result = result.add(val);
            j++;
        } while (val.squaredEuclidianNorm().doubleValue() >= accuracy.doubleValue());
        result = result.multiply(
                valueOf(
                        1 / (sumTettaGamma.multiply(valueOf(GeneratorCreator.factor(valueOf(k)).doubleValue()))).doubleValue()
                )
        );
        result = result.multiply(
                valueOf(
                        Math.pow(gamma.doubleValue() / sumTettaGamma.doubleValue(), k)
                )
        );
        return result;
    }

    private BigDecimalMatrix funcPhiWithWaveWithHat(int n, int k) {
        BigDecimalMatrix result = BigDecimalMatrix.zeroMatrix(systemSize);
        BigDecimalMatrix val;
        int l = k;
        do {
            val = funcPhiWithWave(n, l);
            result = result.add(val);
            l++;
        } while (val.squaredEuclidianNorm().doubleValue() >= accuracy.doubleValue());
        return result;
    }

}
