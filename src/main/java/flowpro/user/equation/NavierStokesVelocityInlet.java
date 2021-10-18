package flowpro.user.equation;

import flowpro.api.ElementData;
import flowpro.api.FlowProProperties;
import java.io.IOException;
import java.util.Arrays;

public class NavierStokesVelocityInlet extends NavierStokes {
	
	double vInlet;
	
	@Override
	public void init(FlowProProperties props) throws IOException {
        super.init(props);
        
        // acceleration
        if(props.containsKey("vInlet")){
            vInlet = props.getDouble("vInlet");
        } else {
            throw new IOException("Velocity in inlet must be defined!");
        }
    }
	
	@Override
    public double[] boundaryValue(double[] WL, double[] n, int TT, ElementData elem) {
		
		switch (TT) {
			case BoundaryType.INLET:
				WL[0] = limiteRho(WL[0]);
				double[] WR = new double[nEqs];
				
				if (!isInletSupersonic) { // subsonic inlet
					double p = pressure(WL);				
					double inletVelocity = vInlet/velocityRef; 
					double rhoIn = 1;				               
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
				
				return WR;
			
			default:
				return super.boundaryValue(WL, n, TT, elem);
		}          
    }
	
	@Override
    public double[] constInitCondition() {		
        if (isInletSupersonic) {
            return WIn;
        } else {
			double p = pOut;
			double inletVelocity = vInlet/velocityRef;
			double rhoIn = 1;				               
			double E = p / (kapa - 1) + 0.5 * rhoIn * inletVelocity * inletVelocity;
			
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
                    vIn[0] = inletVelocity;
                    break;
                case 2:
                    vIn[0] = inletVelocity * Math.cos(alpha);
                    vIn[1] = inletVelocity * Math.sin(alpha);
                    break;
                case 3:
                    vIn[0] = inletVelocity * Math.cos(alpha) * Math.cos(beta);
                    vIn[1] = inletVelocity * Math.sin(alpha) * Math.cos(beta);
                    vIn[2] = inletVelocity * Math.sin(beta);
            }

            double[] W = new double[nEqs];
            W[0] = rhoIn;
            for (int d = 0; d < dim; ++d) {
                W[d + 1] = rhoIn * vIn[d];
            }
            W[dim + 1] = E;

            return W;
        }
    }
}
