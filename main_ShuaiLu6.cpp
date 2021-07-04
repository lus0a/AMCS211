#include<iostream>
#include<iomanip>
using namespace std;

//step time k
int k = 0;
//step size s
double s;
const double sigma = 0.001, epsi = pow(10, -8);

double x, y, df2, rho = 1.0;
double xold, yold;


double f(double x, double y)
{
	double r1 = 3 * x + 4 * y + 0.125 * x * y;
	double r2 = 2 * x + 2 * y + 0.1 * x * y;
	double r3 = 4 * x + 3 * y + x * x / 6;
	return 0.5 * ((r1 - 10) * (r1 - 10) + (r2 - 6) * (r2 - 6) + (r3 - 11) * (r3 - 11));
}
double dfx(double x, double y)
{
	return 0.5 * ((3 + 0.125 * y) * (3 * x + 4 * y + 0.125 * x * y - 10) + (2 + 0.1 * y) * (2 * x + 2 * y + 0.1 * x * y - 6) + (4 + x / 3) * (4 * x + 3 * y + x * x / 6 - 11));
}
double dfy(double x, double y)
{
	return 0.5 * ((4 + 0.125 * x) * (3 * x + 4 * y + 0.125 * x * y - 10) + (2 + 0.1 * x) * (2 * x + 2 * y + 0.1 * x * y - 6) + 3 * (4 * x + 3 * y + x * x / 6 - 11));
}
double dfxx(double x, double y)
{
	return 0.5 * ((3 + 0.125 * y) * (3 + 0.125 * y) + (2 + 0.1 * y) * (2 + 0.1 * y) + (4 + x / 3) * (4 + x / 3));
}
double dfxy(double x, double y)
{
	return 0.5 * ((3 + 0.125 * y) * (4 + 0.125 * x) + (2 + 0.1 * y) * (2 + 0.1 * x) + (4 + x / 3) * 3);
}
double dfyy(double x, double y)
{
	return 0.5 * ((4 + 0.125 * x) * (4 + 0.125 * x) + (2 + 0.1 * x) * (2 + 0.1 * x) + 3 * 3);
}


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

int main()
{
	cout << "Initial guess" << endl << "x0=";
	cin >> x;
	cout << "y0=";
	cin >> y;
	for (int k = 0; k <= 50; k++)
	{
		df2 = sqrt(dfx(x, y) * dfx(x, y) + dfy(x, y) * dfy(x, y));
		cout << scientific << setprecision(3) << "x[" << k << "]= " << x << ' ' << y << endl;
		cout << "f[" << k << "]= " << f(x, y) << endl;
		cout << "df2[" << k << "]=" << df2 << endl;

		////Gradient method
		if (df2 <= epsi)
			break;
		xold = x;
		yold = y;
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

		//A[0] = 1;
		//A[1] = 0;
		//A[2] = A[1];
		//A[3] = 1;

		b[0] = -dfx(xold, yold);
		b[1] = -dfy(xold, yold);
		//for (int i = 1; i <= 4; i++)
		//{
		//	cout << A[i - 1] << endl;
		//}
		//for (int i = 1; i <= 2; i++)
		//{
		//	cout << b[i - 1] << endl;
		//}
		GaussEliPivot(A, b, piv, n);
		SolveUyPivot(A, b, df, piv, n);
		for (int i = 1; i <= n; i++)
		{
			p[i - 1] = df[piv[i - 1]];
			//cout << p[i - 1]<<endl;
		}

		for (int i = 0; i <= 10; i++)
		{
			s = pow(2, -i);
			x = xold + s * p[0];
			y = yold + s * p[1];
			if (f(x, y) <= f(xold, yold) - sigma * s * (dfx(xold, yold) * p[0] + dfy(xold, yold) * p[1]) || i == 10)
				break;
		}
	}
	return 0;
}