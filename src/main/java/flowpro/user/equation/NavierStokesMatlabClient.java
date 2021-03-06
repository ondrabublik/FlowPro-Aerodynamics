package flowpro.user.equation;

import flowpro.api.ElementData;
import flowpro.api.FlowProProperties;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author ales
 */
public class NavierStokesMatlabClient extends Aerodynamics {

    protected static final double P_TOL = 1e-1;
    public MatlabClient mc;

    public void setState(double t, double dt) {
        try {
            mc.computeBoundaryLayer();
        } catch (IOException ex) {
            Logger.getLogger(NavierStokesMatlabClient.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ClassNotFoundException ex) {
            Logger.getLogger(NavierStokesMatlabClient.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    @Override
    public void init(FlowProProperties props) throws IOException {
        int dimension = props.getInt("dimension");
        super.init(props, dimension, dimension + 2, props.getBoolean("isFlowViscous"));

        // Launching Matlab
        try {
            mc = new MatlabClient();
            mc.init();
        } catch (Exception e) {
            System.out.println("Matlab init error " + e);
        }
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

        double[] f = new double[nEqs];

        switch (TT) {
            case (BoundaryType.WALL):
            case (BoundaryType.INVISCID_WALL):
                double p = pressure(WR);
                f[0] = 0;
                double V = .0;
                for (int d = 0; d < dim; ++d) {
                    f[d + 1] = p * n[d];
                    V += elem.meshVelocity[d] * n[d];
                }
                f[dim + 1] = p * V;

                // for ALE
                for (int j = 0; j < nEqs; j++) {
                    f[j] += V * WR[j];
                }
                break;

            case (BoundaryType.INLET):
            case (BoundaryType.OUTLET):
                f = convectiveFlux(WR, n, elem);
                break;

            default: // interior edge
                switch (numericalFluxType) {
                    case "vanLeer":
                        double[] fp = fplus(WL, pressure(WL), n);
                        double[] fm = fminus(WR, pressure(WR), n);
                        for (int j = 0; j < nEqs; j++) {
                            f[j] = fp[j] + fm[j];
                        }
                        break;
                    default:
                        double[] fL = convectiveFlux(WL, n, elem);
                        double[] fR = convectiveFlux(WR, n, elem);
                        double maxEigenValue = Math.max(maxEigenvalue(WL, elem), maxEigenvalue(WR, elem));
                        for (int j = 0; j < nEqs; j++) {
                            f[j] = (fL[j] + fR[j] - maxEigenValue * (WR[j] - WL[j])) / 2;
                        }
                        break;
                }
        }
        return f;
    }

    @Override
    public double[] convectiveFlux(double[] W, double[] n, ElementData elem) {
        W[0] = limiteRho(W[0]);

        double[] f = new double[nEqs];

        double V = .0;
        for (int d = 0; d < dim; ++d) {
            V += W[d + 1] * n[d];
        }
        V /= W[0];

        double p = pressure(W);
        f[0] = W[0] * V;
        for (int d = 0; d < dim; ++d) {
            f[d + 1] = W[d + 1] * V + p * n[d];
        }
        f[dim + 1] = (W[dim + 1] + p) * V;

        return f;
    }

    @Override
    public double[] boundaryValue(double[] WL, double[] n, int TT, ElementData elem) {
        WL[0] = limiteRho(WL[0]);
        double[] WR = new double[nEqs];
        double p = pressure(WL);

        switch (TT) {
            case (BoundaryType.WALL):
                if (isDiffusive) {
                    double[] u = elem.meshVelocity;
                    double absVelocity2 = .0;
                    for (int d = 0; d < dim; ++d) {
                        absVelocity2 += u[d] * u[d];
                    }

                    WR[0] = WL[0];
                    for (int d = 0; d < dim; ++d) {
                        WR[d + 1] = WR[0] * u[d];
                    }
                    WR[dim + 1] = p / (kapa - 1) + WR[0] * absVelocity2 / 2;

                } else {
                    WR = Arrays.copyOf(WL, nEqs);
                    double nu = 0;
                    for (int d = 0; d < dim; ++d) {
                        nu += WL[d + 1] * n[d];
                    }
                    for (int d = 0; d < dim; ++d) { //tangent to wall
                        WR[d + 1] = WL[d + 1] + n[d] * nu;
                    }
                }
                break;
            case (BoundaryType.INVISCID_WALL):
                WR = Arrays.copyOf(WL, nEqs);
                double nu = 0;
                for (int d = 0; d < dim; ++d) {
                    nu += WL[d + 1] * n[d];
                }
                for (int d = 0; d < dim; ++d) { //tangent to wall
                    WR[d + 1] = WL[d + 1] + n[d] * nu;
                }
                break;

            case (BoundaryType.INLET):
                if (!isInletSupersonic) { // subsonic inlet
                    double mach;
                    if (p > pIn0) {
                        mach = 0;
                    } else {
                        mach = Math.sqrt((2 / (kapa - 1)) * (-1 + Math.pow(pIn0 / p, (kapa - 1) / kapa)));
                    }
                    double rhoIn = rhoIn0 * Math.pow((1 + ((kapa - 1) / 2) * mach * mach), 1 / (1 - kapa));
                    double inletVelocity = mach * Math.sqrt((kapa * p) / rhoIn);
                    double E = p / (kapa - 1) + 0.5 * rhoIn * inletVelocity * inletVelocity;

                    WR[0] = rhoIn;
                    if (attackAngle == null) {
                        for (int d = 0; d < dim; ++d) {
                            WR[d + 1] = -rhoIn * inletVelocity * n[d];
                        }
                    } else {
                        double[] dir;
                        if (dim == 2) {
                            dir = new double[]{Math.cos(attackAngle[0]), Math.sin(attackAngle[0])};
                        } else {
                            dir = new double[]{Math.cos(attackAngle[0]) * Math.cos(attackAngle[1]), Math.sin(attackAngle[0]) * Math.cos(attackAngle[1]), Math.sin(attackAngle[1])};
                        }
                        for (int d = 0; d < dim; ++d) {
                            WR[d + 1] = rhoIn * inletVelocity * dir[d];
                        }
                    }
                    WR[dim + 1] = E;
                } else { // supersonic inlet
                    WR = Arrays.copyOf(WIn, nEqs);
                }
                break;

            case (BoundaryType.OUTLET):
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
            pOverRhoDer[d] = (rho * pDer - p * dW[d * nEqs]) / (rho * rho);
        }

        // stress tensor calculation
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

        double constant = kapa / (kapa - 1) / Pr;
        double[] flux = new double[nEqs];
        flux[0] = 0;
        for (int d = 0; d < dim; ++d) {
            double tmp = .0;
            for (int f = 0; f < dim; ++f) {
                flux[f + 1] += stress[dim * d + f] * n[d] / Re;
                tmp += velocity[f] * stress[dim * d + f];
            }
            flux[dim + 1] += (tmp + constant * pOverRhoDer[d]) * n[d] / Re;
        }
        return flux;
    }

    @Override
    public double[] numericalDiffusiveFlux(double Wc[], double dWc[], double[] n, int TT, ElementData elem) {
        Wc[0] = limiteRho(Wc[0]);

        double[] flux = diffusiveFlux(Wc, dWc, n, elem);
        if (TT < 0) {
            if (TT == BoundaryType.WALL) {
                flux[dim + 1] = .0;
            } else {
                Arrays.fill(flux, .0);
            }
        }
        return flux;
    }

    @Override
    public boolean isSourcePresent() {
        return false;
    }

    @Override
    public double[] sourceTerm(double[] W, double[] dW, ElementData elem) { // zdrojovy clen
        throw new UnsupportedOperationException("source is not present");
    }

    protected double[] fplus(double W[], double p, double[] n) {
        double[] f = new double[nEqs];
        double[] velocity = new double[dim];
        double Vn = 0;
        double q2 = 0;
        for (int d = 0; d < dim; ++d) {
            velocity[d] = W[d + 1] / W[0];
            Vn += velocity[d] * n[d];
            q2 += velocity[d] * velocity[d];
        }
        double a = Math.sqrt(kapa * p / W[0]); // rychlost zvuku
        double M = Vn / a;

        if (Math.abs(M) < 1) {
            double fm = W[0] * a * ((M + 1) * (M + 1)) / 4; // mass flux 
            f[0] = fm;
            for (int d = 0; d < dim; ++d) {
                f[d + 1] = fm * (velocity[d] + (-Vn + 2 * a) / kapa * n[d]);
            }
            f[dim + 1] = fm * ((q2 - Vn * Vn) / 2 + (((kapa - 1) * Vn + 2 * a)) * (((kapa - 1) * Vn + 2 * a)) / (2 * (kapa * kapa - 1)));
        } else if (M >= 1) {
            f[0] = W[0] * Vn;
            for (int d = 0; d < dim; ++d) {
                f[d + 1] = W[d + 1] * Vn + p * n[d];
            }
            f[dim + 1] = (W[dim + 1] + p) * Vn;
        }
        return f;
    }

    protected double[] fminus(double W[], double p, double[] n) {
        double[] f = new double[nEqs];
        double[] velocity = new double[dim];
        double Vn = 0;
        double q2 = 0;
        for (int d = 0; d < dim; ++d) {
            velocity[d] = W[d + 1] / W[0];
            Vn += velocity[d] * n[d];
            q2 += velocity[d] * velocity[d];
        }
        double a = Math.sqrt(kapa * p / W[0]); // rychlost zvuku
        double M = Vn / a;

        if (Math.abs(M) < 1) {
            double fm = -W[0] * a * ((M - 1) * (M - 1)) / 4; // mass flux 
            f[0] = fm;
            for (int d = 0; d < dim; ++d) {
                f[d + 1] = fm * (velocity[d] + (-Vn - 2 * a) / kapa * n[d]);
            }
            f[dim + 1] = fm * ((q2 - Vn * Vn) / 2 + (((kapa - 1) * Vn - 2 * a) * ((kapa - 1) * Vn - 2 * a)) / (2 * (kapa * kapa - 1)));
        } else if (M <= -1) {
            f[0] = W[0] * Vn;
            for (int d = 0; d < dim; ++d) {
                f[d + 1] = W[d + 1] * Vn + p * n[d];
            }
            f[dim + 1] = (W[dim + 1] + p) * Vn;
        }
        return f;
    }

    @Override
    public boolean isEquationsJacobian() {
        return false;
    }

    @Override
    public double[] convectiveFluxJacobian(double[] W, double[] n, ElementData elemData) {
        throw new UnsupportedOperationException("convectiveFluxJacobian is not present");
    }

    @Override
    public double[] diffusiveFluxJacobian(double[] W, double[] dW, double n[], ElementData elemData) {
        throw new UnsupportedOperationException("diffusiveFluxJacobian is not present");
    }

    @Override
    public double[] getResults(double[] W, double[] dW, double[] X, String name) {
        switch (name.toLowerCase()) {
            case "temperature":
                double velocity2 = 0;
                for (int i = 0; i < dim; i++) {
                    velocity2 += (W[i + 1] / W[0]) * (W[i + 1] / W[0]);
                }
                return new double[]{velocityRef * velocityRef * (W[dim + 1] / W[0] - velocity2 / 2) / cv};

            case "energy":
                return new double[]{pRef * W[dim + 1]};

            default:
                return super.getResults(W, dW, X, name);
        }
    }

    class MatlabClient {

        Socket socket = null;
        ObjectOutputStream out;
        ObjectInputStream in;

        MatlabClient() {
            try (ServerSocket listener = new ServerSocket(5767)) {
                socket = listener.accept();
                socket.setTcpNoDelay(true);
                socket.setKeepAlive(true);
                out = new ObjectOutputStream(new BufferedOutputStream(socket.getOutputStream()));
                out.writeObject("test");
                out.flush();
                in = new ObjectInputStream(new BufferedInputStream(socket.getInputStream()));
                in.readObject();
                System.out.println("Succesfully connect with Matlab ...");
            } catch (Exception e) {
                System.out.println(e);
            }
        }

        void init() throws IOException, ClassNotFoundException {
            out.writeObject("init");
            out.flush();
            //-----------------------------

            //-----------------------------
            in.readObject();
        }

        void computeBoundaryLayer() throws IOException, ClassNotFoundException {
            out.writeObject("def");
            out.flush();
            //-----------------------------

            out.writeDouble(1.0);
            out.flush();

            //-----------------------------
            in.readObject();
        }

        void nextTimeLevel() throws IOException, ClassNotFoundException {
            out.writeObject("next");
            out.flush();
            in.readObject();
        }

        void close() throws IOException {
            socket.close();
        }
    }
}
