#include<iostream>
#include<iomanip>
using namespace std;

//step time k
int k = 0;
//step size delta
double tao, delta=0.5,  deltamax=1, nta=0.2;
const double epsi = pow(10, -8);

//analytic solution
double xx = 1.0;
double yy = 1.0;

double x, y, e, e2, df2, rho = 1.0;
double ex, ey, e2old, xold, yold, deltaold, gBg, px, py, rhok;

double f(double x, double y)
{
	return (1 - x) * (1 - x) + 100 * (y - x * x) * (y - x * x) ;
}
double dfx(double x, double y)
{
	return -2 * (1 - x) - 400 * (y - x * x) * x;
}
double dfy(double x, double y)
{
	return 200 * (y - x * x);
}
double dfxx(double x, double y)
{
	return 2 - 400 * y + 1200 * x * x;
}
double dfxy(double x, double y)
{
	return -400 * x;
}
double dfyy(double x, double y)
{
	return 200;
}

int main()
{
	cout << "Initial guess" << endl << "x0=";
	cin >> x;
	cout << "y0=";
	cin >> y;
	for (int k = 0; k <= 50; k++)
	{
		ex = x - xx;
		ey = y - yy;
		e2 = sqrt(ex * ex + ey * ey);
		df2 = sqrt(dfx(x, y) * dfx(x, y) + dfy(x, y) * dfy(x, y));
		cout << scientific << setprecision(3) << "x[" << k << "]= " << x << ' ' << y << endl;
		cout << "f[" << k << "]= " << f(x, y) << endl;
		cout << "df2[" << k << "]=" << df2 << endl;
		cout << "e2[" << k << "]=" << e2 << endl;
		if (k >= 1)
			cout << "rho[" << k << "]= " << e2 / e2old << endl;

		xold = x;
		yold = y;
		deltaold = delta;
		gBg = dfxx(xold, yold) * dfx(xold, yold) * dfx(xold, yold) + 2 * dfxy(xold, yold) * dfx(xold, yold) * dfy(xold, yold) + dfyy(xold, yold) * dfy(xold, yold) * dfy(xold, yold);

		if (gBg <= 0)
			tao = 1;
		else {
			if (df2 * df2 * df2 / deltaold / gBg < 1)
				tao = df2 * df2 * df2 / deltaold / gBg;
			else tao = 1;
		}

		px = -tao * delta * dfx(xold, yold) / df2;
		py = -tao * delta * dfy(xold, yold) / df2;
		
		rhok=(f(xold,yold)-f(xold+px,yold+py))/(-dfx(xold,yold)*px - dfy(xold, yold) * py-1/2*(dfxx(xold,yold)*px*px+2*dfxy(xold,yold)*px*py+dfyy(xold,yold)*py*py));
		if (rhok < 0.25)
			delta = 0.25 * deltaold;
		else {
			if ((rhok > 0.75) && (sqrt(px * px + py * py) == deltaold))
				if (2 * deltaold < deltamax)
					delta = 2 * deltaold;
				else
					delta = deltamax;
			else
				delta = deltaold;
		}
		if (rhok > nta)
		{
			x += px;
			y += py;
		}
		e2old = e2;
		if (df2 < epsi)
			break;
	}
	return 0;
}