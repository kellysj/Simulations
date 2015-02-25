package LotkaVolteraRungeKutta;

public class LVRK {
	
	double initialX;
	double initialY;
	double errorLevel;
	double x;
	double y;
	double delXT;
	double delYT;
	double endT=100;
	double h;
	double curT=0;
	String outputY;
	String outputX;
	double alpha =2;
	double beta = 4;
	double delta = 3;
	double gamma = 5;
	public static void main(String[] args[]){

		while(curT<=endT){
			double st = curT;
			double sx = x;
			double sy = y;
			
			while(%error here%){
				double k1x;
				double k1y;
				k1x = h*delXcalc(st,sy,sx);
				k1y = h*delYcalc(st,sy,sx);
				st = t + 1/4*h;
				sx = x + 1/4*k1x;
				sy = y + 1/4*k1y;
				double k2x;
				double k2y;
				k2x = h*delXcalc(st,sy,sx);
				k2y = h*delYcalc(st,sy,sx);
				st = t+(3/8)*h;
				sx = x+(3/32)*k1x+(9/32)*k2x;
				sy = y+(3/32)*k1y+(9/32)*k2y;
				double k3x;
				double k3y;
				k3x = h*delXcalc(st,sy,sx);
				k3y = h*delYcalc(st,sy,sx);
				st = t+(12/13)*h;
				sx = x+(1932/2197)*k1x-(7200/2197)*k2x+(7296/2197)*k3x;
				sy = y+(1932/2197)*k1y-(7200/2197)*k2y+(7296/2197)*k3y;
				double k4x;
				double k4y;
					k4x = h*delXcalc(st,sy,sx);
					k4y = h*delYcalc(st,sy,sx);
				st = t+h;
				sy = y+(439/216)*k1y-8*k2y+(3680/513)*k3y-(845/4104)*k4y;
				sx = x+(439/216)*k1x-8*k2x+(3680/513)*k3x-(845/4104)*k4x;
				double k5x;
				double k5y;
					k5x = h*delXcalc(st,sy,sx);

					k5y = h*delYcalc(st,sy,sx);
				st = t+(h/2);
				sx = x - (8/27)*k1x+2*k2x-(3544/2565)*k3x+(1859/4104)*k4x-(11/44)*k5x;
				sy = y - (8/27)*k1y+2*k2y-(3544/2565)*k3y+(1859/4104)*k4y-(11/44)*k5y;
				double k6x;
				double k6y;
				
					k6x = h*delXcalc(st,sy,sx);
					k6y = h*delYcalc(st,sy,sx);
				}
				
				
			}
		}
		
	}
	public double delXcalc(double t, double y, double x){
		double out = 0;
		out = alpha*x - beta*x*y;
		return out;
	}
	public double delYcalc(){
		double out = 0;
		out = delta*x*y-gamma*y;
		return out;
	}

}    