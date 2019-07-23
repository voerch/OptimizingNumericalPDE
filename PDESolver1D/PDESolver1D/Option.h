#pragma once

#include <algorithm>


const double PI = 4.0 * atan(1.0);



class VanillaOption 
{

public:
	virtual double Payoff(double sharePrice) = 0;

	double strike;
	double interestRate;
	double timeToExpiry;
	double sigma;

	
	VanillaOption(double strike_, double interestRate_, double timeToExpiry_, double sigma_) :strike(strike_), interestRate(interestRate_), timeToExpiry(timeToExpiry_), sigma(sigma_) {};


	virtual double PriceByBS(double S0) = 0;
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
};

class EurCall : public VanillaOption
{
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

//class EurPut : public VanillaOption
//{
//public:
//	double Payoff(double sharePrice)
//	{
//		return std::max(strike - sharePrice , 0.0);
//	}
//
//	EurPut(double strike_, double interestRate_, double timeToExpiry_, double sigma_) : VanillaOption(strike_, interestRate_, timeToExpiry_, sigma_) {}
//
//};





