package core;


import java.math.BigDecimal;
import java.math.RoundingMode;

import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;
import static java.math.BigDecimal.valueOf;


public class ErgodicityCondition {

    //Approved
    public static boolean check(double gamma, BigDecimal lambda, Integer K, BigDecimal T) {
        BigDecimal minuend = lambda.multiply(T);
        if (K == 1) {
            K++;
        }
        BigDecimal subtrahend = getY0(new BigDecimal(gamma), T, K).multiply(new BigDecimal(lambda.doubleValue() / gamma));
        return minuend.add(subtrahend).compareTo(new BigDecimal(1)) == -1;
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
        BigDecimal[] psiVector = new BigDecimal[K];
        psiVector[0] = ONE;

        for (int k = 0; k < K - 1; k++) {
            BigDecimal sum = ZERO;
            for (int s = 1; s <= k; s++) {
                sum = sum.add(psiVector[s].multiply(smallPhi(gamma, T, valueOf(k + 1 - s).intValue())));
            }
            psiVector[k + 1] = new BigDecimal(psiVector[k]
                    .subtract(smallPhi(gamma, T, valueOf(k).intValue()))
                    .subtract(sum).doubleValue()
                    / smallPhi(gamma, T, valueOf(0).intValue()).doubleValue());
        }
        return psiVector;
    }

    private static BigDecimal smallPhi(BigDecimal gamma, BigDecimal t, int k) {
        return new BigDecimal(gamma.multiply(t).pow(k).doubleValue() / factor(k)
                * Math.exp(gamma.multiply(t).doubleValue() * -1));
    }

    //Approved
    public static BigDecimal factor(BigDecimal n) {
        if(n.signum() < 0) {
            throw new IllegalArgumentException("Negative value in factorial function.");
        }
        if(n.scale() < 0) {
            throw new IllegalArgumentException("Real value in factorial function.");
        }
        if(n.equals(ZERO)) {
            return ONE;
        }
        return n.multiply(factor(n.subtract(ONE)));
    }

    private static int factor(int n) {
        int result = 1;
        for (int i = 1; i <= n; i++) {
            result = result * i;
        }
        return result;
    }
}