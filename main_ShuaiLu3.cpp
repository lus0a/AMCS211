#include<iostream>
#include<iomanip>
using namespace std;

//step time k
int k = 0;
//step size s
double s;
const double sigma = 0.001, epsi = pow(10, -8);
//analytic solution
double xx = -17.0 / 8;
double yy = -15.0 / 8;

double x, y, e, e2, df2, rho = 1.0;
double ex, ey, e2old, xold, yold;

double f(double x, double y)
{
	return 5 * x * x + 5 * y * y - 6 * x * y + 10 * x + 6 * y + 5;
}
double dfx(double x, double y)
{
	return 10 * x - 6 * y + 10;
}
double dfy(double x, double y)
{
	return 10 * y - 6 * x + 6;
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

		//Coordinate Search method
		//x -= s * dfx(x, y);
		//y -= s * dfy(x, y);


		////Gradient method
		if (df2 <= epsi)
			break;
		xold = x;
		yold = y;
		for (int i = 0; i <= 8; i++)
		{
			s = pow(2, -i);
			x = xold - s * dfx(xold, yold);
			y = yold - s * dfy(xold, yold);
			if (f(x, y) <= f(xold, yold) - sigma * s * (dfx(xold, yold) * dfx(xold, yold) + dfy(xold, yold) * dfy(xold, yold)) || i == 8)
				break;
		}
		e2old = e2;
	}
	return 0;
}