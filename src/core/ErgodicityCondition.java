package core;


import java.math.BigDecimal;

import static java.math.BigDecimal.*;


public class ErgodicityCondition {

    //Approved
    public static boolean check(BigDecimal gamma, BigDecimal lambda, int K, BigDecimal T) {
        BigDecimal minuend = lambda.multiply(T);
        BigDecimal subtrahend = getY0(gamma, T, K).multiply(valueOf(lambda.doubleValue() / gamma.doubleValue()));
        return minuend.add(subtrahend).doubleValue() < 1;
    }

    public static double getLambdaRestriction(BigDecimal t, BigDecimal gamma, int K) {
        return 1 / (t.doubleValue() + getY0(gamma, t, K).doubleValue() / gamma.doubleValue());
    }
    private static BigDecimal getY0(BigDecimal gamma, BigDecimal T, Integer K) {
        BigDecimal[] psiVector = getPsiVector(gamma, T, K);
        BigDecimal sum = ZERO;
        for (BigDecimal element : psiVector) {
            sum = sum.add(element);
        }
        return new BigDecimal(1 / sum.doubleValue());
    }

    //Approved
    private static BigDecimal[] getPsiVector(BigDecimal gamma, BigDecimal T, Integer K) {
        BigDecimal[] psiVector = new BigDecimal[K + 1];
        psiVector[0] = ONE;

        for (int k = 0; k <= K - 1; k++) {
            BigDecimal sum = ZERO;
            for (int s = 1; s <= k; s++) {
                sum = sum.add(psiVector[s].multiply(GeneratorCreator.funcSmallPhiK(gamma, T, k + 1 - s)));
            }
            psiVector[k + 1] = new BigDecimal(psiVector[k]
                    .subtract(GeneratorCreator.funcSmallPhiK(gamma, T, k))
                    .subtract(sum).doubleValue()
                    / GeneratorCreator.funcSmallPhiK(gamma, T, 0).doubleValue());
        }
        return psiVector;
    }
}