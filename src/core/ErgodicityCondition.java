package core;

import kurs.BigDecimalMatrix;

import java.math.BigDecimal;
import java.math.RoundingMode;

import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;
import static java.math.BigDecimal.valueOf;

/**
 * Created by Lenovo on 09.05.2018.
 */
public class ErgodicityCondition {

    //Approved
    public static boolean check(BigDecimal gamma, BigDecimal lambda, Integer K, BigDecimal T) {
        BigDecimal minuend = lambda.multiply(T);
        BigDecimal subtrahend = getY0(gamma, T, K).multiply(lambda.divide(gamma, RoundingMode.HALF_UP));
        return minuend.add(subtrahend).compareTo(ONE) < 0;
    }

    private static BigDecimal getY0(BigDecimal gamma, BigDecimal T, Integer K) {
        BigDecimal[] psiVector = getPsiVector(gamma, T, K);
        BigDecimal sum = ZERO;
        for (BigDecimal element : psiVector) {
            sum = sum.add(element);
        }
        return ONE.divide(sum, RoundingMode.HALF_UP);
    }

    //Approved
    private static BigDecimal[] getPsiVector(BigDecimal gamma, BigDecimal T, Integer K) {
        BigDecimal[] psiVector = new BigDecimal[K];
        psiVector[0] = ONE;

        for (int k = 0; k < K - 1; k++) {
            BigDecimal sum = ZERO;
            for (int s = 1; s <= k; s++) {
                sum = sum.add(psiVector[s].multiply(smallPhi(gamma, T, valueOf(k + 1 - s))));
            }
            psiVector[k + 1] = psiVector[k].subtract(smallPhi(gamma, T,valueOf(k))).subtract(sum);
        }
        return psiVector;
    }

    /*private static BigDecimalMatrix getY0(BigDecimalMatrix d0, BigDecimalMatrix d1, BigDecimal z, int k) {
        BigDecimalMatrix minusDz = (d1.multiply(z).add(d0)).multiply(new BigDecimal(-1));
        BigDecimalMatrix gammaI = new BigDecimalMatrix(2, new BigDecimal(1), 10).multiply(z);
        BigDecimalMatrix sumInverse = minusDz.add(gammaI).inverse();
        BigDecimalMatrix val = sumInverse;
        for(int i = 0; i < k; i++) {
            val = val.multiply(sumInverse);
        }
        return val.multiply(z.pow(k));
    }*/

    //Approved
    private static BigDecimal smallPhi(BigDecimal gamma, BigDecimal T, BigDecimal k) {
        return (gamma.multiply(T).pow(k.intValue()))
                .divide(factor(k), RoundingMode.HALF_UP)
                .multiply(new BigDecimal(Math.exp(gamma.multiply(T).doubleValue())));
    }

    //Approved
    private static BigDecimal factor(BigDecimal n) {
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
}