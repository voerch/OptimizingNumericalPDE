#include "pde.h"
#include "fdm.h"

int main()
{

	double x_dom = 1.0;       // Spot goes from [0.0, 1.0]
	long J = 5;
	double t_dom = 0.075;         // Time period as for the option
	long N = 4;

	HeatEqn* heat_pde = new HeatEqn;

	ExplicitMethod euler(x_dom, J, t_dom, N, heat_pde);

	euler.stepMarch();

	delete heat_pde;

	system("pause");
	return 0;

}