package flowpro.user.equation;

import flowpro.api.ElementData;
import flowpro.api.FlowProProperties;
import static flowpro.user.equation.Aerodynamics.RHO_TOL;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

/**
 *
 * @author ales
 */
public class KOmega extends Aerodynamics {

    protected static final double P_TOL = 1e-1;

    double kIn;
    double omIn;

    double LIM_TOL;

    // parameters of turbulence model
    double bk;
    double bom;
    double aom;
    double Prt;
    double sk;
    double som;

    @Override
    public void init(FlowProProperties props) throws IOException {
        int dimension = props.getInt("dimension");
        super.init(props, dimension, dimension + 4, true);

        // inlet
        kIn = props.getDouble("kIn");
        omIn = props.getDouble("omIn");
        if (isInletSupersonic) {
            WIn[dimension + 2] = kIn;
            WIn[dimension + 3] = omIn;
        }

        // parameters of turbulence model
        bk = props.getDouble("bk");
        bom = props.getDouble("bom");
        aom = props.getDouble("aom");
        Prt = props.getDouble("Prt");
        sk = props.getDouble("sk");
        som = props.getDouble("som");
    }

    @Override
    public double pressure(double[] W) {
        W[0] = limiteRho(W[0]);

        double momentum2 = .0;
        for (int d = 0; d < dim; ++d) {
            momentum2 += W[d + 1] * W[d + 1];
        }

        double p = (kapa - 1) * (W[dim + 1] - momentum2 / (2 * W[0]));
        if (p < P_TOL) {
            p = P_TOL * Math.exp((p - P_TOL) / P_TOL);
        }
        return p;
    }

    @Override
    public void limitUnphysicalValues(double[] Ws, double[] W, int nBasis) { // limituje zaporne hodnoty
        if (Ws[0] < RHO_TOL) {
            for (int j = 0; j < nBasis; j++) {
                W[j] = RHO_TOL;
            }
        }

//        double momentum2 = .0;
//        for (int d = 0; d < dim; ++d) {
//            momentum2 += Ws[d+1] * Ws[d+1];
//        }
//        double Ek = momentum2 / (2 * Ws[0]);
//        if (Ws[dim+1] < Ek) {
//            for (int j = 0; j < nBasis; j++) {
//                W[j][dim+1] = Ek;
//            }
//        }
    }

    @Override
    public void saveReferenceValues(String filePath) throws IOException {
        FlowProProperties output = new FlowProProperties();

        output.setProperty("l", Double.toString(lRef));
        output.setProperty("p", Double.toString(pRef));
        output.setProperty("rho", Double.toString(rhoRef));
        output.setProperty("v", Double.toString(velocityRef));
        output.setProperty("t", Double.toString(tRef));

        output.store(new FileOutputStream(filePath), null);
    }

    @Override
    public double[] getReferenceValues() {
        return new double[]{lRef, pRef, rhoRef, velocityRef, tRef};
    }

    @Override
    public double[] constInitCondition() {
        if (isInletSupersonic) {
            return WIn;
        } else {
            double machIn = Math.sqrt(2 / (kapa - 1) * (Math.pow((1 / pOut), (kapa - 1) / kapa) - 1));
            double rhoIn = Math.pow(1 + ((kapa - 1) / 2) * machIn * machIn, 1 / (1 - kapa));
            double VIn = machIn * Math.sqrt((kapa * pOut) / rhoIn);

            double[] vIn = new double[dim];
            double alpha = 0;
            double beta = 0;
            if (attackAngle != null) {
                alpha = attackAngle[0];
                if (dim == 3) {
                    beta = attackAngle[1];
                }
            }
            switch (dim) {
                case 1:
                    vIn[0] = VIn;
                    break;
                case 2:
                    vIn[0] = VIn * Math.cos(alpha);
                    vIn[1] = VIn * Math.sin(alpha);
                    break;
                case 3:
                    vIn[0] = VIn * Math.cos(alpha) * Math.cos(beta);
                    vIn[1] = VIn * Math.sin(alpha) * Math.cos(beta);
                    vIn[2] = VIn * Math.sin(beta);
            }

            double[] W = new double[nEqs];
            W[0] = rhoIn;
            for (int d = 0; d < dim; ++d) {
                W[d + 1] = rhoIn * vIn[d];
            }
            W[dim + 1] = pOut / (kapa - 1) + 0.5 * rhoIn * VIn * VIn;
            W[dim + 2] = rhoIn * kIn;
            W[dim + 3] = rhoIn * omIn;

            return W;
        }
    }

    @Override
    public double[] numericalConvectiveFlux(double WL[], double WR[], double Vs, double[] n, int TT, ElementData elem) {
        WL[0] = limiteRho(WL[0]);
        WR[0] = limiteRho(WR[0]);

        double[] f = new double[nEqs];

        switch (TT) {
            case (Aerodynamics.BoundaryType.WALL):
            case (Aerodynamics.BoundaryType.INVISCID_WALL):
                double p = pressure(WL);
                f[0] = 0;
                for (int d = 0; d < dim; ++d) {
                    f[d + 1] = p * n[d];
                }
                f[dim + 1] = p * Vs;
                f[dim + 2] = 0;
                f[dim + 3] = 0;
                break;

            case (Aerodynamics.BoundaryType.INLET):
            case (Aerodynamics.BoundaryType.OUTLET):
                f = convectiveFlux(WR, Vs, n, elem);
                break;

            default: // interior edge
                double[] fL = convectiveFlux(WL, Vs, n, elem);
                double[] fR = convectiveFlux(WR, Vs, n, elem);
                double maxEigenValue = Math.max(maxEigenvalue(WL), maxEigenvalue(WR));
                for (int j = 0; j < nEqs; j++) {
                    f[j] = (fL[j] + fR[j] - maxEigenValue * (WR[j] - WL[j])) / 2;
                }
                break;
        }
        return f;
    }

    @Override
    public double[] convectiveFlux(double[] W, double Vs, double[] n, ElementData elem) {
        W[0] = limiteRho(W[0]);

        double[] f = new double[nEqs];

        double V = .0;
        for (int d = 0; d < dim; ++d) {
            V += W[d + 1] * n[d];
        }
        V /= W[0];

        double p = pressure(W);
        f[0] = W[0] * (V - Vs);
        for (int d = 0; d < dim; ++d) {
            f[d + 1] = W[d + 1] * (V - Vs) + p * n[d];
        }
        f[dim + 1] = (W[dim + 1] + p) * V - W[dim + 1] * Vs;
        f[dim + 2] = W[dim + 2] * (V - Vs);
        f[dim + 3] = W[dim + 3] * (V - Vs);
        return f;
    }

    @Override
    public double[] boundaryValue(double[] WL, double[] u, double[] n, int TT, ElementData elem) {
        WL[0] = limiteRho(WL[0]);
        double[] WR = new double[nEqs];
        double p = pressure(WL);

        switch (TT) {
            case (Aerodynamics.BoundaryType.WALL):
                if (isDiffusive) {
                    double absVelocity2 = .0;
                    for (int d = 0; d < dim; ++d) {
                        absVelocity2 += u[d] * u[d];
                    }

                    WR[0] = WL[0];
                    for (int d = 0; d < dim; ++d) {
                        WR[d + 1] = WR[0] * u[d];
                    }
                    WR[dim + 1] = p / (kapa - 1) + WR[0] * absVelocity2 / 2;
                    WR[dim + 2] = 0;
                    WR[dim + 3] = WL[dim + 3];
                    break;
                }
            case (Aerodynamics.BoundaryType.INVISCID_WALL):
                WR = Arrays.copyOf(WL, nEqs);
                break;

            case (Aerodynamics.BoundaryType.INLET):
                if (!isInletSupersonic) { // subsonic inlet
                    double mach;
                    if (p > pIn0) {
                        mach = 0;
                    } else {
                        mach = Math.sqrt((2 / (kapa - 1)) * (-1 + Math.pow(pIn0 / p, (kapa - 1) / kapa)));
                    }
                    double rhoIn = rhoIn0 * Math.pow((1 + ((kapa - 1) / 2) * mach * mach), 1 / (1 - kapa));
                    double normalVelocity = mach * Math.sqrt((kapa * p) / rhoIn);
                    double E = p / (kapa - 1) + 0.5 * rhoIn * normalVelocity * normalVelocity;

                    WR[0] = rhoIn;
                    for (int d = 0; d < dim; ++d) {
                        WR[d + 1] = -rhoIn * normalVelocity * n[d];
                    }
                    WR[dim + 1] = E;
                    WR[dim + 2] = rhoIn * kIn;
                    WR[dim + 3] = rhoIn * omIn;
                } else { // supersonic inlet
                    WR = Arrays.copyOf(WIn, nEqs);
                }
                break;

            case (Aerodynamics.BoundaryType.OUTLET):
                double rhoOut = WL[0];
                double absVelocity2 = .0;
                for (int d = 0; d < dim; ++d) {
                    absVelocity2 += WL[d + 1] * WL[d + 1];
                }
                absVelocity2 /= rhoOut * rhoOut;

                double a = Math.sqrt(kapa * pressure(WL) / WL[0]);
                double mach = Math.sqrt(absVelocity2) / a;

                if (mach < 1) {  // subsonic outlet
                    WR[0] = rhoOut;
                    for (int d = 0; d < dim; ++d) {
                        WR[d + 1] = rhoOut * WL[d + 1] / WL[0];
                    }
                    WR[dim + 1] = pOut / (kapa - 1) + rhoOut * absVelocity2 / 2;
                    WR[dim + 2] = WL[dim + 2];
                    WR[dim + 3] = WL[dim + 3];
                } else {  // supersonic outlet
                    WR = Arrays.copyOf(WL, nEqs);
                }
        }
        return WR;
    }

    @Override
    public double[] diffusiveFlux(double[] W, double[] dW, double[] n, ElementData elem) {
        W[0] = limiteRho(W[0]);
        double rho = W[0];

        double lam = -2. / 3; // Stokesuv vztah

        double[] velocity = new double[dim];
        double velocity2 = .0;
        for (int d = 0; d < dim; ++d) {
            velocity[d] = W[d + 1] / rho;
            velocity2 += velocity[d] * velocity[d];
        }

        double[] velocityJac = new double[dim * dim];
        for (int d = 0; d < dim; ++d) {
            for (int f = 0; f < dim; ++f) {
                velocityJac[dim * d + f] = (dW[f * nEqs + d + 1] - dW[f * nEqs] * velocity[d]) / rho;
            }
        }

        double p = pressure(W);
        double[] pOverRhoDer = new double[dim];
        for (int d = 0; d < dim; ++d) {
            double temp = .0;
            for (int f = 0; f < dim; ++f) {
                temp += velocity[f] * velocityJac[dim * f + d];
            }
            double pDer = (kapa - 1) * (dW[d * nEqs + dim + 1] - dW[d * nEqs] * velocity2 / 2 - rho * temp);
            pOverRhoDer[d] = (rho * pDer - p * dW[d * nEqs]) / (rho * rho); // rho * rho * p * p !!!!!!!!!!!!!!!!
        }

        // turbulentStress tensor calculation
        double[] stress = new double[dim * dim];
        double trace = .0;
        for (int d = 0; d < dim; ++d) {
            trace += velocityJac[dim * d + d];
            for (int f = 0; f < dim; ++f) {
                stress[dim * d + f] = velocityJac[dim * d + f] + velocityJac[dim * f + d];
            }
        }
        for (int d = 0; d < dim; ++d) {
            stress[dim * d + d] += lam * trace;
        }

        double k = W[dim + 2] / rho;
        double om = W[dim + 3] / rho;
        double mut = max(0, W[dim + 2] / Math.exp(om));

        double[] turbulentStress = new double[dim * dim];
        for (int i = 0; i < stress.length; i++) {
            turbulentStress[i] = stress[i] * mut;
        }
        double kMax = 2.0 / 3 * max(0, W[dim + 2]);
        for (int d = 0; d < dim; ++d) {
            turbulentStress[dim * d + d] -= kMax;
        }

        double[] kDer = new double[dim];
        double[] omDer = new double[dim];
        for (int d = 0; d < dim; ++d) {
            kDer[d] = (dW[d * nEqs + dim + 2] - dW[d * nEqs] * k) / rho;
            omDer[d] = (dW[d * nEqs + dim + 3] - dW[d * nEqs] * om) / rho;
        }

        double constant = kapa / (kapa - 1) * (1 / Pr + mut / Prt);
        double[] flux = new double[nEqs];
        for (int d = 0; d < dim; ++d) {
            double tmp = .0;
            for (int f = 0; f < dim; ++f) {
                flux[f + 1] += (stress[dim * d + f] + turbulentStress[dim * d + f]) * n[d] / Re;
                tmp += velocity[f] * (stress[dim * d + f] + turbulentStress[dim * d + f]);
            }
            flux[dim + 1] += (tmp + constant * pOverRhoDer[d]) * n[d] / Re;
            flux[dim + 2] += (1 + sk * mut) / Re * kDer[d] * n[d];
            flux[dim + 3] += (1 + som * mut) / Re * omDer[d] * n[d];
        }

        return flux;
    }

    @Override
    public double[] numericalDiffusiveFlux(double WL[], double WR[], double dWL[], double dWR[],
            double[] n, int TT, ElementData elem) {
        WL[0] = limiteRho(WL[0]);
        WR[0] = limiteRho(WR[0]);

        double[] fluxL = diffusiveFlux(WL, dWL, n, elem);
        double[] fluxR = diffusiveFlux(WR, dWR, n, elem);
        double[] flux = new double[nEqs];
        for (int j = 0; j < nEqs; j++) {
            flux[j] = (fluxL[j] + fluxR[j]) / 2;
        }

        if (TT < 0) {
            if (TT == Aerodynamics.BoundaryType.WALL) {
                flux[dim + 1] = .0;
            } else {
                Arrays.fill(flux, .0);
            }
        }
        return flux;
    }

    @Override
    public boolean isSourcePresent() {
        return true;
    }

    @Override
    public double[] sourceTerm(double[] W, double[] dW, ElementData elem) { // zdrojovy clen
        W[0] = limiteRho(W[0]);
        double rho = W[0];

        double lam = -2. / 3; // Stokesuv vztah

        double[] velocity = new double[dim];
        for (int d = 0; d < dim; ++d) {
            velocity[d] = W[d + 1] / rho;
        }

        double[] velocityJac = new double[dim * dim];
        for (int d = 0; d < dim; ++d) {
            for (int f = 0; f < dim; ++f) {
                velocityJac[dim * d + f] = (dW[f * nEqs + d + 1] - dW[f * nEqs] * velocity[d]) / rho;
            }
        }

        double k = W[dim + 2] / rho;
        double om = W[dim + 3] / rho;
        double expOm = Math.exp(om);
        double mut = max(0, W[dim + 2] / expOm);

        // turbulentStress tensor calculation
        double[] turbulentStress = new double[dim * dim];
        double trace = .0;
        for (int d = 0; d < dim; ++d) {
            trace += velocityJac[dim * d + d];
            for (int f = 0; f < dim; ++f) {
                turbulentStress[dim * d + f] = mut * (velocityJac[dim * d + f] + velocityJac[dim * f + d]);
            }
        }
        double kMax = 2.0 / 3 * max(0, W[dim + 2]);
        for (int d = 0; d < dim; ++d) {
            turbulentStress[dim * d + d] += mut * lam * trace - kMax;
        }

        double omDerSqr = 0;
        for (int d = 0; d < dim; ++d) {
            double omDer = (dW[d * nEqs + dim + 3] - dW[d * nEqs] * om) / rho;
            omDerSqr += omDer * omDer;
        }

        double Tv = 0;
        for (int d = 0; d < dim; d++) {
            for (int f = 0; f < dim; f++) {
                Tv += turbulentStress[dim * d + f] * velocityJac[dim * d + f];
            }
        }

        if (Tv > 10 * bk * rho * k * expOm) {
            Tv = 10 * bk * rho * k * expOm;
        }

        double Tvom = aom / max(1e-8, k) * Tv;
        if (Tvom > 100 * bom * rho * expOm) {
            Tvom = 100 * bom * rho * expOm;
        }

        double[] source = new double[nEqs];
        source[dim + 2] = max(Tv - bk * rho * k * expOm,0);
        source[dim + 3] = max(Tvom - bom * rho * expOm + 1 / Re * (1 + som * mut) * omDerSqr,0);
        return source;
    }

    @Override
    public double[] getResults(double[] W, double[] X, String name) {
        switch (name) {
            case "mach":
                double absVelocity = .0;
                for (int d = 0; d < dim; ++d) {
                    absVelocity += W[d + 1] * W[d + 1];
                }
                absVelocity = Math.sqrt(absVelocity) / W[0];

                double a = Math.sqrt(kapa * pressure(W) / W[0]);
                return new double[]{absVelocity / a};

            case "density":
                return new double[]{rhoRef * W[0]};

            case "xVelocity":
                return new double[]{velocityRef * W[1] / W[0]};

            case "yVelocity":
                if (dim > 1) {
                    return new double[]{velocityRef * W[2] / W[0]};
                } else {
                    throw new UnsupportedOperationException("undefined value" + name);
                }

            case "zVelocity":
                if (dim > 3) {
                    return new double[]{velocityRef * W[3] / W[0]};
                } else {
                    throw new UnsupportedOperationException("undefined value" + name);
                }

            case "velocity":
                double[] velocity = new double[dim];
                for (int i = 0; i < dim; i++) {
                    velocity[i] = velocityRef * W[i + 1] / W[0];
                }
                return velocity;

            case "temperature":
                throw new UnsupportedOperationException("undefined value" + name);

            case "energy":
                return new double[]{pRef * W[dim + 1]};

            case "pressure":
                return new double[]{pRef * pressure(W)};

            case "k":
                return new double[]{W[dim + 2] / W[0]};

            case "omega":
                return new double[]{W[dim + 3] / W[0]};

            case "turbulentViscosity":
                return new double[]{W[dim + 2] / Math.exp(W[dim + 3] / W[0]) / Re};

            default:
                throw new UnsupportedOperationException("undefined value " + name);
        }
    }

    public double max(double a, double b) {
        if (a > b) {
            return a;
        } else {
            return b;
        }
    }
}
