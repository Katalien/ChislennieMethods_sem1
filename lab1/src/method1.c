#include<stdio.h>
#include <math.h>

double f(double x) {
	return pow(x,2)+2*x-19;
}

double root(double a, double b, double eps) {


}


double bisection_method(double a, double b, double(*f)(double),  double eps) {
	int k = 0;
	 double x, xnow;
	 double c=0;
	 double dif = fabs(b - a);
	while (dif > 2* eps) {
		k++;
		c = (a + b) / 2;
		if ((*f)(c) == 0)
			return c;
		printf("a = %lf, b = %lf, f(a) = %lf, f(b) = %lf, f(c)=%lf \n", a, b, f(a), f(b), f(c));
		if ((*f)(a) * (*f)(c) < 0)
			b = c;
		else
			a = c;
		xnow = (a + b) / 2;

		//fall = fabs(x0 - xnow);
		//fallmass[k] = fall;
		dif = fabs(b - a);
	}
	return xnow;
}

double fixed_point_method(double a, double b) {
	double alpha = a + (b-a)/3;
	double beta = b - (b-a)/3;


}

 double wegstein_method (double prev02, double prev01, double itprev, double eps, double(*f)(double)) {
	 double itnow = 0;

}
 ///////////////////////////////////
 double f3(double x) { return x * x * log10(x) - 3.8; }

 //int main()
 //{
	// double x0, x1, f0, f1, eps;
	// int mt;
	// printf("Input x0, xl, epsilon, max itexation\n");
	// //scanf("%lf %lf %lf %d", &x0, &x1, &eps, &mt);
	// x0 = -6;
	// x1 = -4;
	// eps = 0.01;
	// f0 = f(x0);
	// f1 = f(x1);
	// int res;
	// while (fabs(x1 - x0) > eps)
	// {
	//	 double x = x1 - f1 * (x1 - x0) / (f1 - f0);
	//	 f0 = f1;
	//	 f1 = f(x);
	//	 x0 = x1;
	//	 x1 = x;

	//	
	// }
	// x0 = -6;
	// x1 = -4;
	// printf("Solution: x = %lf\n", x1);
	// res = bisection_method(x0, x1, f, eps);
	//printf("res = %f", res);
 //}


int main() {
	double a1, b1, a2, b2, c;
	 double eps = 0.01;
	 double res1, res2;
	int k=0;
	a1 = -6.0;
	b1 = -4.0;
	a2 = 2 ;
	b2 = 4;

	 double x0 = -5.472;
	 double fall ;
	 double fallmass[10];
	double dif = fabs(b1 - a1);

	res1 = bisection_method(a1, b1, f, eps);
	res2 = bisection_method(a2, b2, f, eps);
	

	printf("res2 = %f\n", res2);
	//printf("res2 = %f\n", res2);
}