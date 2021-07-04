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
double ex, ey, e2old, xold, yold, deltaold, gBg, px, py, pux, puy, rhok;

void GaussEli(double* A, double* b, int n)
{
	for (int k = 1; k <= n - 1; k++)
	{
		double a = A[(k - 1) * n + (k - 1)];
		for (int i = k + 1; i <= n; i++)
		{
			double aa = A[(i - 1) * n + (k - 1)] / a;
			for (int j = k; j <= n; j++)
			{
				A[(i - 1) * n + (j - 1)] -= A[(k - 1) * n + (j - 1)] * aa;
			}
			b[i - 1] -= b[k - 1] * aa;
		}
	}
}

void GaussEliPivot(double* A, double* b, int* piv, int n)
{
	for (int p = 1; p <= n; p++)
		piv[p - 1] = p - 1;

	for (int k = 1; k <= n - 1; k++)
	{

		int MaxID = k - 1;
		for (int p = k + 1; p <= n; p++)
		{
			if (abs(A[piv[MaxID] * n + (k - 1)]) < abs(A[piv[p - 1] * n + (k - 1)]))
				MaxID = p - 1;
		}
		if (MaxID != (k - 1))
		{
			int temp;
			temp = piv[k - 1];
			piv[k - 1] = piv[MaxID];
			piv[MaxID] = temp;
		}
		double a = A[piv[k - 1] * n + (k - 1)];
		for (int i = k + 1; i <= n; i++)
		{
			double aa = A[piv[i - 1] * n + (k - 1)] / a;
			for (int j = k; j <= n; j++)
			{
				A[piv[i - 1] * n + (j - 1)] -= A[piv[k - 1] * n + (j - 1)] * aa;
			}
			b[piv[i - 1]] -= b[piv[k - 1]] * aa;
		}
	}
}

void SolveUyPivot(double* U, double* y, double* x, int* piv, int n)
{
	for (int i = n; i >= 1; i--)
	{
		x[piv[i - 1]] = y[piv[i - 1]] / U[piv[i - 1] * n + (i - 1)];
		for (int j = 1; j <= i - 1; j++)
		{
			y[piv[j - 1]] -= U[piv[j - 1] * n + (i - 1)] * x[piv[i - 1]];
		}
	}
}

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

		pux = -df2 * df2 / gBg * dfx(xold, yold);
		puy = -df2 * df2 / gBg * dfy(xold, yold);
		if (sqrt(pux * pux + puy * puy) < delta)
		{
			int n = 2;
			double* A = new double[n * n];
			double* b = new double[n];
			double* df = new double[n];
			double* p = new double[n];
			int* piv = new int[n];

			A[0] = dfxx(xold, yold);
			A[1] = dfxy(xold, yold);
			A[2] = A[1];
			A[3] = dfyy(xold, yold);
			b[0] = -dfx(xold, yold);
			b[1] = -dfy(xold, yold);
			GaussEliPivot(A, b, piv, n);
			SolveUyPivot(A, b, df, piv, n);
			for (int i = 1; i <= n; i++)
			{
				p[i - 1] = df[piv[i - 1]];
			}
			if (sqrt(p[0] * p[0] + p[1] * p[1]) < delta)
			{
				px = p[0];
				py = p[1];
			}
			else 
			{
				double a = (p[0] - pux) * (p[0] - pux) + (p[1] - puy) * (p[1] - puy);
				double b = (p[0] - pux) * pux + (p[1] - puy) * puy;
				double c = pux * pux + puy * puy - delta * delta;
				px = pux + (-b + sqrt(b * b - 4 * a * c)) / (2 * a) * (p[0] - pux);
				py = puy + (-b + sqrt(b * b - 4 * a * c)) / (2 * a) * (p[1] - puy);
			}

		}
		else
		{
			px = delta / sqrt(pux * pux + puy * puy) * pux;
			py = delta / sqrt(pux * pux + puy * puy) * puy;
		}
		
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