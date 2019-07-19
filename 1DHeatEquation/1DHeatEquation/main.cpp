#include "PDE.h"
#include "FDM.h"
#include "Option.h"

int main()
{

	//double x_dom = 1.0;       // Spot goes from [0.0, 1.0]
	//long J = 20;
	//double t_dom = 0.075;         // Time period as for the option
	//long N = 20;

	//HeatEqn* heat_pde = new HeatEqn;

	//ExplicitMethod euler(x_dom, J, t_dom, N, heat_pde);
	//euler.stepMarch();

	//delete heat_pde;

	//double K = 0.5;  // Strike price
	//double r = 0.05;   // Risk-free rate (5%)
	//double v = 0.2;    // Volatility of the underlying (20%)
	//double T = 1.00;    // One year until expiry

	//double x_dom = 1.0;       // Spot goes from [0.0, 1.0]
	//unsigned long J = 20;
	//double t_dom = T;         // Time period as for the option
	//unsigned long N = 20;

	//VanillaOption* callOption = new EurCall(K, r, T, v);

	//BlackScholesPDE* bsPDE = new BlackScholesPDE(callOption);

	//ExplicitMethod euler(x_dom, J, t_dom, N, bsPDE);
	//euler.stepMarch();
	//delete bsPDE;




	system("pause");
	return 0;

}