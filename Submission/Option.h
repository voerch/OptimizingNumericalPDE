/*
Author: Mustafa Berke Erdis, August 2019

Header file for storing option data and analytical solution of the options. European call option is implemented.
*/

#pragma once
#include <algorithm>
const double PI = 4.0 * atan(1.0);

class VanillaOption 
{

public:
	double strike;
	double interestRate;
	double timeToExpiry;
	double sigma;

	virtual double Payoff(double sharePrice) = 0;
	virtual double PriceByBS(double S0) = 0;

	VanillaOption(double strike_, double interestRate_, double timeToExpiry_, double sigma_) :strike(strike_), interestRate(interestRate_), timeToExpiry(timeToExpiry_), sigma(sigma_) {};
};

class EurCall : public VanillaOption
{
private:
	double d_plus(double S0)
	{
		return (log(S0 / strike) + (interestRate + 0.5*pow(sigma, 2.0))*timeToExpiry) / (sigma*sqrt(timeToExpiry));
	}

	double d_minus(double S0)
	{
		return d_plus(S0) - sigma*sqrt(timeToExpiry);
	}

	double NormDist(double x)
	{
		return 0.5 * erfc(-x * sqrt(0.5));
	}

public:
	double Payoff(double sharePrice)
	{
		return std::max(sharePrice - strike, 0.0);
	}

	double PriceByBS(double S0)
	{
		return S0*NormDist(d_plus(S0)) - strike*exp(-interestRate*timeToExpiry)*NormDist(d_minus(S0));
	}

	EurCall(double strike_, double interestRate_, double timeToExpiry_, double sigma_) : VanillaOption(strike_, interestRate_, timeToExpiry_, sigma_) {}


};