#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <fstream>
#include <math.h>
#include "Exception.hpp"
#include "Quote.hpp"
#include "Stock.hpp"

void loadData(Stock&, std::string);
void updateMoney(double&, double, double, bool);

void updateMoney(double& m, double past, double future, bool buy) {
	if (buy) {
		m*=future/past;
		std::cout<<"Long. Money: "<<m<<"$. Change: "<<100*(future/past-1)<<"%.\n";
	} else {
		m*=past/future;
		std::cout<<"Short. Money: "<<m<<"$. Change: "<<100*(past/future-1)<<"%.\n";
	}
}

void loadData(Stock& S, std::string file) {
	/*

	Loads all historical data available in the given file into
	the object S, of class Stock.

	*/

	std::ifstream input0(file.c_str());
	assert(input0.is_open());
	int lines=0;
	std::string date;
	double price;
	while (!input0.eof()) {
		input0 >> date >> price;
		lines+=1;
	}
	input0.close();
	S=Stock("S&P 500", lines);

	std::ifstream input(file.c_str());

	assert(input.is_open());
	while (!input.eof()) {
		input >> date >> price;
		S.addData(Quote(S.getName(), price, date));
	}
	input.close();

}

void get_ret(double* ret, double* vec, int size){
	for(int i = 0; i < size-1; ++i){
		ret[i] = (vec[i+1] - vec[i]) / vec[i];
	}
}

double get_vec_mean(double* vec, int size){
	double vec_sum = 0;
	for(int i = 0; i < size; ++i){
		vec_sum += vec[i];
	}
	return (vec_sum / size);
}

double get_vec_std(double* vec, int size){
	
	double vec_pro = 0;
	double vec_mean = get_vec_mean(vec, size);
	for(int i = 0; i < size; ++i){
		vec_pro += (vec[i] - vec_mean) * (vec[i] - vec_mean);
	}
	return (sqrt(vec_pro / size));	
}

double get_avg_year_ret(double value, int size){
	double avg_ret;
	avg_ret = pow(value, (double)252 / size) - 1;
	return avg_ret;	
}

double get_sharpe_rat(double ret_mean, double* ret, int size){
	double ret_std = get_vec_std(ret, size);
	return((ret_mean - 0.02) / (sqrt(252) * ret_std));
}

void get_max_draw_down(double* result, double* value, double par, int size){
	double max_draw_down;
	int len = size;
	int low_index = 0;
	int peak_index = 0;
	int max_draw_down_period = 1;
	int li, pi;
	double low = value[0];
	double peak = value[0];
	for(int i = 1; i < len; ++i){
		if(value[i] > peak){
			peak = value[i];
			low = value[i];
			peak_index = i;
			low_index = i;
		}
		else if((value[i] < low) && ((peak - value[i]) > max_draw_down)) {
			low = value[i];
			low_index = i;
			max_draw_down = peak - low;
			max_draw_down_period = low_index - peak_index + 1;	
			li = i;
			pi = peak_index;
		}
	}
		result[0] = max_draw_down * par;
		result[1] = (double)(li - pi + 1);
		result[2] = (double)pi;
		result[3] = (double)li;
}
	
int main(int argc, char* argv[])
{	
	
	// We load the historical data for S&P 500 and NASDAQ
	Stock SP500, NASDAQ;
	loadData(SP500, "SP500.dat");
	loadData(NASDAQ, "NASDAQ100.dat");

	std::ofstream outfile1("C:\\Users\\Codyj\\OneDrive\\Documents\\CPP Scripts\\Final_Project\\pair-trading-master\\result.txt");
	outfile1 << "Date,Worth" << std::endl;

	const double INITIAL_MONEY = 1000.0; // Start with 1000$
	double money = INITIAL_MONEY; // Portfolio worth
	double portfolio[3] = {1.0, -1.0, 1}; // The portfolio {x, y, z} indicates weights in the portfolio of SP500, NASDAQ and cash, respectively.
	const int DAYS_CONSIDERED = 30; // Number of days we consider
	int size = SP500.getFilled() - DAYS_CONSIDERED + 1;
	double value[size];
	for(int i = DAYS_CONSIDERED; i < SP500.getFilled()-1; i++) {
		// Update portfolio worth
		money*=portfolio[2]+portfolio[0]*SP500.getData(i).getPrice()/SP500.getData(i-1).getPrice()+portfolio[1]*NASDAQ.getData(i).getPrice()/NASDAQ.getData(i-1).getPrice();
		value[i - DAYS_CONSIDERED] = money;
		// Show current state of portfolio.
		std::cout<<"> Portfolio weights: \n\t"<<100*portfolio[0]<<"% in SP500\n\t"<<100*portfolio[1]<<"% in NASDAQ\n\t"<<100*portfolio[2]<<"% in cash\n"<<"  Worth: "<<money<<"$\n  Return: "<<(money/INITIAL_MONEY-1)*100<<"%\n  Date: "<<SP500.getData(i).getDate();
		std::cout<<"\n\n";
		outfile1 << SP500.getData(i).getDate() << "," << money << std::endl;
  		
		double mean=0, sd=0;
		for (int j=i-DAYS_CONSIDERED; j<i; j++) {
			mean+=SP500.getData(j).getPrice()/NASDAQ.getData(j).getPrice();
			sd+=pow(SP500.getData(j).getPrice()/NASDAQ.getData(j).getPrice(),2);		
		}
		mean/=DAYS_CONSIDERED;
		sd=sqrt(sd/DAYS_CONSIDERED-pow(mean, 2));
		const double MARGIN=sd*2;
		if (SP500.getData(i).getPrice()/NASDAQ.getData(i).getPrice()>mean+MARGIN) {
			// SP500 is outperforming NASDAQ, as their ratio is big. Thus, short SP500 and buy NASDAQ.
			portfolio[0]=-1.0;
			portfolio[1]=1.0;
			portfolio[2]=1.0;
		} else if (SP500.getData(i).getPrice()/NASDAQ.getData(i).getPrice()<mean-MARGIN) {
			// SP500 is underperforming NASDAQ, as their ratio is small. Thus, buy SP500 and short NASDAQ.
			portfolio[0]=1.0;
			portfolio[1]=-1.0;
			portfolio[2]=1.0;
		} else {
			// The ratio of SP500 and NASDAQ is inside the "normal" values.
			portfolio[0]=0.0;
			portfolio[1]=0.0;
			portfolio[2]=1.0;
		}
	}
	
	outfile1.close();
	
	// calculate avg_return
	double ret[size-1]; 
	get_ret(ret,value,size);
	double avg_ret_port = get_avg_year_ret(money/INITIAL_MONEY,size-1);
	double sharpe_ret_port = get_sharpe_rat(get_vec_mean(ret,size-1),ret,size-1);
	double max_draw_down_port[4]; 
	get_max_draw_down(max_draw_down_port,value,1,size);
		
		
	std::cout<<"\n---------------------------------------------\n";
	std::cout<< "Results\n\n"<<"  Start date: "<<SP500.getData(0).getDate()<<"\n"<<"  End date: "<<SP500.getData(SP500.getFilled()-2).getDate()<<"\n";
	std::cout<<"  Return: "<<(money/INITIAL_MONEY-1)*100<<"%";
	std::cout<<"\n---------------------------------------------\n";
	
	
	std::cout<<"\n";
}
