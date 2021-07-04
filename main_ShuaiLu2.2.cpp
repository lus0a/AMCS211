#include<iostream>
#include<cmath>
#include <iomanip>

using namespace std;
#define PI acos(-1)

double x = PI / 3;
double h0 = 0.1;
double h, dplusf, d0f;

double f(double x)
{
	return sin(x);
}

double dfx(double x)
{
	return cos(x);
}

int main()
{
	cout << "Left difference error" << "   " << "Central difference error" << endl;
	for (int k = 0; k <= 6; k++)
	{
		h = h0 * pow(10, -k);
		if (k == 6)
			h = pow(10, -12);
		dplusf = (f(x + h) - f(x)) / h;
		d0f = (f(x + h) - f(x - h)) / (2 * h);
		cout << scientific << setprecision(8) << dplusf - dfx(x) << "         " << d0f - dfx(x) << endl;
	}

	return 0;
}