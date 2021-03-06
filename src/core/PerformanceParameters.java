package core;


import java.math.BigDecimal;
import java.util.ArrayList;

import kurs.BigDecimalMatrix;

import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;
import static java.math.BigDecimal.valueOf;

/**
 * Created by user on 04.06.2017.
 */
public class PerformanceParameters {
    private ArrayList<BigDecimalMatrix> piVector;
    private int piVectorSize;
    private int piSize;
    private BigDecimal lambda;
    private BigDecimalMatrix d1;
    private BigDecimalMatrix e;
    private BigDecimalMatrix eCol;
    private int systemSize;
    private ArbitararyTimeGenerator generator;

    public PerformanceParameters(
            ArrayList<BigDecimalMatrix> piVector,
            BigDecimal lambda,
            BigDecimalMatrix d1,
            int systemSize,
            ArbitararyTimeGenerator generator
    ) {
        this.piVector = piVector;
        this.piVectorSize = piVector.size();
        this.d1 = d1;
        this.piSize = piVector.get(0).getWidth();
        this.lambda = lambda;
        this.e = BigDecimalMatrix.eCol(this.piSize, ONE);
        this.eCol = BigDecimalMatrix.eCol(systemSize, ONE);
        this.systemSize = systemSize;
        this.generator = generator;
    }

    public BigDecimal requestCorrectPerfomanceFunction(int number, int K) {
        switch (number) {
            case 11: {
                return getAverageNumberOfRequests();
            }
            case 12: {
                return getAverageNumberOfEnergyUnits();
            }
            case 13: {
                return getNoRequestsProbability();
            }
            case 14: {
                return getNoEnergyUnitsProbability();
            }
            case 15: {
                return getNoEnergySystemIdleProbability();
            }
            case 16: {
                return getAverageNumberOfRequestsInArbitryTime(K);
            }
            default: {
                return ZERO;
            }
        }
    }

    public BigDecimal getAverageNumberOfRequestsInArbitryTime(int K) {
        BigDecimal sum = ZERO;
        // System.out.println(this.piVectorSize);
        BigDecimalMatrix eCol = BigDecimalMatrix.eCol(systemSize * (K + 1), ONE);
        for (int j = 0; j < this.piVectorSize; j++) {
            sum = sum.add(generator.calculateP(j).multiply(eCol).multiply(valueOf(j)).getElement(0, 0));
        }
        return sum;
    }

    public BigDecimal getAverageNumberOfRequests() {
        BigDecimal sum = ZERO;
        // System.out.println(this.piVectorSize);
        for (int i = 0; i < this.piVectorSize; i++) {
            sum = sum.add(piVector.get(i).multiply(e).multiply(valueOf(i)).getElement(0, 0));
        }
        return sum;
    }

    public BigDecimal getAverageNumberOfEnergyUnits() {
        BigDecimal sum = ZERO;
        BigDecimalMatrix piMatrix = BigDecimalMatrix.eRow(systemSize, ZERO);
        for (int i = 0; i < this.piVectorSize; i++) {
            for (int k = 1; k <= this.piSize / systemSize - 1; k++) {
                if (systemSize == 2) {
                    piMatrix.setElement(0, 0, piVector.get(i).getElement(0, 2 * k));
                    piMatrix.setElement(0, 1, piVector.get(i).getElement(0, 2 * k + 1));
                }
                else {
                    piMatrix.setElement(0, 0, piVector.get(i).getElement(0, k));
                }
                sum = sum.add(valueOf(k).multiply(piMatrix.multiply(eCol).getElement(0, 0)));
            }
        }
        return sum;
    }

    // It doesnt work
    public BigDecimal getEnergyUnitLossProbability() {
        BigDecimal sum = ZERO;
        for (int i = 0; i < this.piVectorSize; i++) {
            sum = sum.add(piVector.get(i).getElement(0, piSize - 1));
        }
        return sum;
    }

    public BigDecimal getNoRequestsProbability() {
        return piVector.get(0).multiply(e).getElement(0, 0);
    }

    public BigDecimal getNoEnergyUnitsProbability() {
        BigDecimal sum = ZERO;
        BigDecimalMatrix piMatrix = BigDecimalMatrix.eRow(systemSize, ZERO);
        for (int i = 0; i < this.piVectorSize; i++) {
            if (systemSize == 2) {
                piMatrix.setElement(0, 0, piVector.get(i).getElement(0, 0));
                piMatrix.setElement(0, 1, piVector.get(i).getElement(0, 1));
            }
            else {
                piMatrix.setElement(0, 0, piVector.get(i).getElement(0, 0));
            }
            sum = sum.add(piMatrix.multiply(eCol).getElement(0, 0));
        }
        return sum;
    }

    public BigDecimal getNoEnergySystemIdleProbability() {
        BigDecimal sum = ZERO;
        BigDecimalMatrix piMatrix = BigDecimalMatrix.eRow(systemSize, ZERO);
        for (int i = 1; i < this.piVectorSize; i++) {
            if (systemSize == 2) {
                piMatrix.setElement(0, 1, piVector.get(i).getElement(0, 1));
                piMatrix.setElement(0, 0, piVector.get(i).getElement(0, 0));
            }
            else {
                piMatrix.setElement(0, 0, piVector.get(i).getElement(0, 0));
            }
            sum = sum.add(piMatrix.multiply(eCol).getElement(0, 0));
        }
        return sum;
    }

    //Запрос начнет обслуживаться в момент прихода
    public BigDecimal getRequestStartWorkedWhenHisWent() {
        BigDecimal sum = ZERO;
        BigDecimalMatrix piMatrix = BigDecimalMatrix.eRow(systemSize, ZERO);

        for (int k = 0; k <= this.piSize / systemSize - 1; k++) {
            piMatrix.setElement(0, 0, piVector.get(0).getElement(0, 2 * k));
            piMatrix.setElement(0, 1, piVector.get(0).getElement(0, 2 * k + 1));
            sum = sum.add(piMatrix.multiply(d1).multiply(eCol).getElement(0, 0));
        }

        return new BigDecimal(1 / lambda.doubleValue() * sum.doubleValue());
    }
}
