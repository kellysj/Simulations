package LotkaVolteraRungeKutta;

import java.util.ArrayList;

public class LVRK {

	double x = 10;
	double y = 10;
	double h = 1;
	double t = 0;

	double alpha = .3;
	double beta = .5;
	double delta = .5;
	double gamma = .3;
	
	double upperError = .7;
	double lowerError = .0;
	double endT = 10000;
	
	ArrayList<Double> xStore;
	ArrayList<Double> yStore;
	ArrayList<Double> tStore;

	public static void main(String args[]) {
		LVRK run = new LVRK();
		run.RKcalc();
		//System.out.println(run.printToMatlabMatrix("x"));
	}

	private double delXcalc(double t, double y, double x) {
		double out = 0;
		out = alpha * x - beta * x * y;
		return out;
	}

	private double delYcalc(double t, double y, double x) {
		double out = 0;
		out = delta * x * y - gamma * y;
		return out;
	}

	public void RKcalc() {
		xStore = new ArrayList<Double>();
		yStore = new ArrayList<Double>();
		tStore = new ArrayList<Double>();
		long calcTime = System.nanoTime();
		while (t <= endT) {
			//System.out.println("t = " + t + " x = " + x + " y = " + y);
			double error = upperError * 2;
			xStore.add(x);
			yStore.add(y);
			tStore.add(t);
			while (!((error <= upperError) && (error >= lowerError))) {
				//K CALCULATIONS START HERE
				double st = t;
				double sy = y;
				double sx = x;
				double k1x = h * delXcalc(st, sy, sx);
				double k1y = h * delYcalc(st, sy, sx);
				st = t + 1 / 4 * h;
				sy = y + 1 / 4 * k1y;
				sx = x + 1 / 4 * k1x;
				double k2x = h * delXcalc(st, sy, sx);
				double k2y = h * delYcalc(st, sy, sx);
				st = t + (3 / 8) * h;
				sy = y 
						+ (3 / 32) * k1y 
						+ (9 / 32) * k2y;
				sx = x 
						+ (3 / 32) * k1x 
						+ (9 / 32) * k2x;
				double k3x = h * delXcalc(st, sy, sx);
				double k3y = h * delYcalc(st, sy, sx);
				st = t + (12 / 13) * h;
				sy = y 
						+ (1932 / 2197) * k1y 
						- (7200 / 2197) * k2y
						+ (7296 / 2197) * k3y;//coefficients cause large change in values of k...
				sx = x 
						+ (1932 / 2197) * k1x 
						- (7200 / 2197) * k2x
						+ (7296 / 2197) * k3x;
				double k4x = h * delXcalc(st, sy, sx);
				double k4y = h * delYcalc(st, sy, sx);
				st = t + h;
				sy = y 
						+ (439 / 216) * k1y 
						- 8 * k2y 
						+ (3680 / 513) * k3y
						- (845 / 4104) * k4y;
				sx = x 
						+ (439 / 216) * k1x 
						- 8 * k2x 
						+ (3680 / 513) * k3x
						- (845 / 4104) * k4x;
				double k5x = h * delXcalc(st, sy, sx);
				double k5y = h * delYcalc(st, sy, sx);
				st = t + (h / 2);
				sy = y 
						- (8 / 27) * k1y 
						+ 2 * k2y - (3544 / 2565) * k3y
						+ (1859 / 4104) * k4y 
						- (11 / 44) * k5y;
				sx = x 
						- (8 / 27) * k1x 
						+ 2 * k2x - (3544 / 2565) * k3x
						+ (1859 / 4104) * k4x 
						- (11 / 44) * k5x;
				double k6x = h * delXcalc(st, sy, sx);
				double k6y = h * delYcalc(st, sy, sx);
				//K CALCULATIONS END
				error = Math.abs(
								(
								(1.0 / 360.0) * (k1y + k1x) 
								- (128.0 / 4275.0) * (k3y + k3x) 
								- (2197.0 / 75240.0) * (k4y + k4x) 
								+ (1.0 / 50.0) * (k5y + k5x) 
								+ (2.0 / 55.0) * (k6y + k6x)
								) / 2 // averaging coupled k values
							)
						* Math.pow(h, 6);
				//ADAPTIVE STEP STUFF
				if (error > upperError) {
					h = h / 2;

				} else if (error < lowerError) {
					h = h * 2;

				}
				//END OF LOOP CONDITIONS
				else if (((error <= upperError) && (error >= lowerError))) {
					x = x + (
							(16.0 / 135.0) * k1x 
							+ (6656.0 / 12825.0) * k3x
							+ (28561.0 / 56430.0) * k4x 
							- (9.0 / 50.0)* k5x 
							+ (2.0 / 55.0) * k6x
							)* Math.pow(h, 6);
					y = y + (
							(16.0 / 135.0) * k1y 
							+ (6656.0 / 12825.0) * k3y
							+ (28561.0 / 56430.0) * k4y 
							- (9.0 / 50.0) * k5y 
							+ (2.0 / 55.0) * k6y
							)* Math.pow(h, 6);
					t = t + h;
				}
			}
		}
		calcTime = System.nanoTime() - calcTime;
		System.out.println("TIME TO RUN: " + calcTime + "ns");
	}
	public String printToMatlabMatrix(String var){
		if(var.equals("t")){
			return matLabFormat(tStore);
		}
		if(var.equals("x")){
			return matLabFormat(xStore);
		}
		if(var.equals("y")){
			return matLabFormat(yStore);
		}
		else return null;
	}
	private String matLabFormat(ArrayList<Double> in){
		String out;
		out = "[";
		for(int i = 0;i<in.size();i++){
			out = out + in.get(i);
			if(i==in.size()-1){
				out += "]";
			}
			else{
				out += " ";
			}
		}
		return out;
	}
}