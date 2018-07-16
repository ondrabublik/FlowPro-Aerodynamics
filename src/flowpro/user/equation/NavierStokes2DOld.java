package flowpro.user.equation;

import flowpro.api.ElementData;
import flowpro.api.FlowProProperties;
import flowpro.api.Mat;
import static flowpro.user.equation.Aerodynamics.RHO_TOL;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

/**
 *
 * @author ales
 */
public class NavierStokes2DOld extends Aerodynamics {

    protected static final double P_TOL = 1e-1;

    @Override
    public void init(FlowProProperties props) throws IOException {
        int dimension = props.getInt("dimension");
        super.init(props, dimension, dimension + 2, props.getBoolean("isFlowViscous"));
    }

    @Override
    public double pressure(double[] W) {

        double momentum2 = .0;
        for (int d = 0; d < dim; ++d) {
            momentum2 += W[d + 1] * W[d + 1];
        }

        double p = (kapa - 1) * (W[dim + 1] - momentum2 / (2 * W[0]));
        return p;
    }

    @Override
    public void limitUnphysicalValues(double[] Ws, double[] W, int nBasis) { // limituje zaporne hodnoty
//        if (Ws[0] < RHO_TOL) {
//            for (int j = 0; j < nBasis; j++) {
//                W[j] = RHO_TOL;
//            }
//        }
//
//        double momentum2 = .0;
//        for (int d = 0; d < dim; ++d) {
//            momentum2 += Ws[d + 1] * Ws[d + 1];
//        }
//        double Ek = momentum2 / (2 * Ws[0]);
//        if (Ws[dim + 1] < Ek) {
//            for (int j = 0; j < nBasis; j++) {
//                W[(dim + 1) * nBasis + j] = Ek;
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

            return W;
        }
    }

    @Override
    public double[] numericalConvectiveFlux(double WL[], double WR[], double[] n, int TT, ElementData elem) {
        WL[0] = limiteRho(WL[0]);
        WR[0] = limiteRho(WR[0]);

        double[] fn = new double[nEqs];
        double[] fcL, fcR;
        double p, pL, pR;

        switch (TT) {
            case (-1): // stena
                p = pressure(WL);
                fn[0] = 0;
                fn[1] = p * n[0];
                fn[2] = p * n[1];
                fn[3] = 0;
                break;

            case (-2): // vstup
                fn = convectiveFlux(WR, n, elem);
                break;

            case (-3): // vystup
                fn = convectiveFlux(WR, n, elem);
                break;

            case (-4): // nevazka stena
                p = pressure(WL);
                fn[0] = 0;
                fn[1] = p * n[0];
                fn[2] = p * n[1];
                fn[3] = 0;
                break;

            default: // vnitrni stena
                pL = pressure(WL);
                pR = pressure(WR);
                double aL = Math.sqrt(kapa * pL / WL[0]);
                double aR = Math.sqrt(kapa * pR / WR[0]);
                double VnL = WL[1] / WL[0] * n[0] + WL[2] / WL[0] * n[1];
                double VnR = WR[1] / WR[0] * n[0] + WR[2] / WR[0] * n[1];
                fcL = convectiveFlux(WL, n, elem);
                fcR = convectiveFlux(WR, n, elem);
                double S = Math.max(Math.abs(VnL) + aL, Math.abs(VnR) + aR);
                for (int j = 0; j < nEqs; j++) {
                    fn[j] = (fcL[j] + fcR[j] - S * (WR[j] - WL[j])) / 2;
                }
                break;
        }

        return fn;
    }

    @Override
    public double[] convectiveFlux(double[] W, double[] n, ElementData elem) {
        double[] f = new double[nEqs];
        double V = W[1] / W[0] * n[0] + W[2] / W[0] * n[1];
        double p = pressure(W);
        f[0] = W[0] * V;
        f[1] = W[1] * V + p * n[0];
        f[2] = W[2] * V + p * n[1];
        f[3] = (W[3] + p) * V;

        return f;
    }

    @Override
    public double[] boundaryValue(double[] WL, double[] n, int TT, ElementData elem) {
        double[] WR = new double[nEqs];
        double p = pressure(WL);
        switch (TT) {
            case (-1): // stena
                if (Re != -1) {
                    WR[0] = WL[0];
                    WR[1] = 0;
                    WR[2] = 0;
                    WR[3] = WL[3];
                } else {
                    WR[0] = WL[0];
                    WR[1] = WL[1];
                    WR[2] = WL[2];
                    WR[3] = WL[3];
                }
                break;

            case (-4):  // nevazka stena
                WR = Arrays.copyOf(WL, nEqs);
                break;

            case (-2):
                // osetreni NaN
                if (pIn0 > p) {
                    double Minl = Math.sqrt((2 / (kapa - 1)) * (-1 + Math.pow(pIn0 / p, (kapa - 1) / kapa)));
                    double Rinl = rhoIn0*Math.pow((1+((kapa-1)/2)*Minl*Minl),1/(1-kapa));
                    double Vinl = Minl * Math.sqrt((kapa * p) / Rinl);
                    double uinl = Vinl * Math.cos(attackAngle[0]);
                    double vinl = Vinl * Math.sin(attackAngle[0]);
                    double Einl = p / (kapa - 1) + 0.5 * Rinl * Vinl * Vinl;
                    WR[0] = Rinl;
                    WR[1] = Rinl * uinl;
                    WR[2] = Rinl * vinl;
                    WR[3] = Einl;
                } else {
                    WR[0] = WL[0];
                    WR[1] = 0;
                    WR[2] = 0;
                    WR[3] = p / (kapa - 1);
                }
                break;

            case (-3):
                double ro = WL[0];
                double uo = WL[1] / WL[0];
                double vo = WL[2] / WL[0];
                double ao = Math.sqrt(kapa * p / ro);
                double Mo = Math.sqrt(uo * uo + vo * vo) / ao;

                if (Mo < 1) { // subsonicky vystup
                    double Eo, pout;
                    pout = pOut;
                    Eo = pout / (kapa - 1) + ro * (uo * uo + vo * vo) / 2;
                    WR[0] = ro;
                    WR[1] = ro * uo;
                    WR[2] = ro * vo;
                    WR[3] = Eo;
                } else { // supersonicky vystup
                    WR = Arrays.copyOf(WL, nEqs);
                }
                break;
        }
        return WR;
    }

    @Override
    public double[] diffusiveFlux(double[] W, double[] dW, double[] n, ElementData elem) {
        int nr = nEqs;
        double[] fvn = new double[nr];
        double r = W[0];
        double u = W[1] / r;
        double v = W[2] / r;
        double E = W[3];
        double rx = dW[0];
        double ry = dW[nr];
        double ux = 1 / r * (dW[1] - rx * u);
        double uy = 1 / r * (dW[nr + 1] - ry * u);
        double vx = 1 / r * (dW[2] - rx * v);
        double vy = 1 / r * (dW[nr + 2] - ry * v);
        double Ex = dW[3];
        double Ey = dW[nr + 3];
        double Tx = (Ex - E*rx/r - r*u*ux - r*v*vx)/r;
        double Ty = (Ey - E*ry/r - r*u*uy - r*v*vy)/r;

        double Sxx = ux - 1.0 / 3 * (ux + vy);
        double Sxy = 0.5 * (vx + uy);
        double Syy = vy - 1.0 / 3 * (ux + vy);

        double txx = 2 * Sxx;
        double txy = 2 * Sxy;
        double tyy = 2 * Syy;

        fvn[0] = 0;
        fvn[1] = 1 / Re * txx * n[0] + 1 / Re * txy * n[1];
        fvn[2] = 1 / Re * txy * n[0] + 1 / Re * tyy * n[1];
        fvn[3] = (1 / Re * (u * txx + v * txy + 1 / (kapa - 1) / Pr * Tx)) * n[0] + (1 / Re * (u * txy + v * tyy + 1 / (kapa - 1) / Pr * Ty)) * n[1];

        return fvn;
    }

    @Override
    public double[] numericalDiffusiveFlux(double WL[], double WR[], double dWL[], double dWR[],
            double[] n, int TT, ElementData elem) {

        double[] fvnL;
        double[] fvnR;

        fvnL = diffusiveFlux(WL, dWL, n, elem);
        fvnR = diffusiveFlux(WR, dWR, n, elem);

        double[] fvn = new double[nEqs];
        for (int j = 0; j < nEqs; j++) {
            fvn[j] = (fvnL[j] + fvnR[j]) / 2;
        }

        if (TT < 0) {
            if (TT == -1) {
                fvn[3] = 0;
            } else {
                for (int j = 0; j < nEqs; j++) {
                    fvn[j] = 0;
                }
            }
        }

        return fvn;
    }

    @Override
    public boolean isSourcePresent() {
        return false;
    }

    @Override
    public double[] sourceTerm(double[] W, double[] dW, ElementData elem) { // zdrojovy clen
        throw new UnsupportedOperationException("source is not present");
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
                if (dim > 2) {
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
                double velocity2 = 0;
                for (int i = 0; i < dim; i++) {
                    velocity2 += (W[i + 1] / W[0]) * (W[i + 1] / W[0]);
                }
                return new double[]{velocityRef * velocityRef * (W[dim + 1] / W[0] - velocity2 / 2) / cv};

            case "energy":
                return new double[]{pRef * W[dim + 1]};

            case "pressure":
                return new double[]{pRef * pressure(W)};

            default:
                throw new UnsupportedOperationException("undefined value " + name);
        }
    }

    @Override
    public boolean isEquationsJacobian() {
        return false;
    }

    @Override
    public double[] convectiveFluxJacobian(double[] W, double[] n, ElementData elemData) {
        return null;
    }

    @Override
    public double[] diffusiveFluxJacobian(double[] W, double[] dW, double n[], ElementData elemData) {
        return null;
    }
}
