#pragma once

#include <algorithm>

class VanillaOption 
{

public:
	virtual double Payoff(double sharePrice) = 0;

	double strike;
	double interestRate;
	double timeToExpiry;
	double sigma;

	VanillaOption(double strike_, double interestRate_, double timeToExpiry_, double sigma_) :strike(strike_), interestRate(interestRate_), timeToExpiry(timeToExpiry_), sigma(sigma_) {};

};

class EurCall : public VanillaOption
{
public:
	double Payoff(double sharePrice)
	{
		return std::max(sharePrice - strike, 0.0);
	}

	EurCall(double strike_, double interestRate_, double timeToExpiry_, double sigma_) : VanillaOption(strike_, interestRate_, timeToExpiry_, sigma_) {}


};

class EurPut : public VanillaOption
{
public:
	double Payoff(double sharePrice)
	{
		return std::max(strike - sharePrice , 0.0);
	}

	EurPut(double strike_, double interestRate_, double timeToExpiry_, double sigma_) : VanillaOption(strike_, interestRate_, timeToExpiry_, sigma_) {}

};