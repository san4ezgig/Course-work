package core;

import kurs.BigDecimalMatrix;

import java.math.BigDecimal;

import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;

public class PoisonStream {
    private BigDecimal gamma;
    private BigDecimal lambda;
    private Integer K;
    private BigDecimal accuracy;
    private BigDecimalMatrix d0;
    private BigDecimalMatrix d1;
    private BigDecimal T;
    private BigDecimal tetta;

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
        this.elementCreator = new GeneratorCreator.MatrixElementCreator();
    }
}
