#include<iostream>
using namespace std;

//step time k
int k = 0;
//step size alpha
double alpha = 0.1;
//analytic solution
double xx = -17.0 / 8;
double yy = -15.0 / 8;

double x, y, e, e2, rho = 1.0;
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
	for (int k = 0; k <= 20; k++)
	{
		ex = x - xx;
		ey = y - yy;
		e2 = sqrt(ex * ex + ey * ey);
		cout << "x[" << k << "]= " << x << ' ' << y << endl;
		cout << "e[" << k << "]= " << ex << ' ' << ey << endl;
		cout << "e2[" << k << "]=" << e2 << endl;
		if (k >= 1)
			cout << "rho[" << k << "]= " << e2 / e2old << endl;

		//Coordinate Search method
		//x -= s * dfx(x, y);
		//y -= s * dfy(x, y);


		////Gradient method

		xold = x;
		yold = y;
		alpha = (dfx(xold, yold) * dfx(xold, yold) + dfy(xold, yold) * dfy(xold, yold)) / (10 * dfx(xold, yold) * dfx(xold, yold) + 10 * dfy(xold, yold) * dfy(xold, yold) - 12 * dfx(xold, yold) * dfy(xold, yold));
		x -= alpha * dfx(xold, yold);
		y -= alpha * dfy(xold, yold);



		e2old = e2;
	}
	return 0;
}