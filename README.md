# Pair trading

Pair trading is a trading strategy that takes two stocks that are supposed to be highly correlated, and takes a short position when one of the stocks outperforms the other, and takes a long position in the other. As the stocks are highly correlated, the strategy makes the assumption that they will come back to their correlation, so that ideally the shorted stock would go down in price, and the other stock would go up. Moreover, as both stocks are highly correlated, this is a zero-beta strategy, as we are shorting a stock, while buying the other one.

# The code

The code has two main classes, Quote and Stock. In Quote, the information of the quote of a certain stock, in a certain date is stored. The class Stock is used to store the prices of a certain stock, for a period of time.

For the strategy, I have used S&P 500 and NASDAQ 100 indices, as they are highly correlated. The stretegy is relatively simple: for each day, the mean of the quotients of the price of S&P 500 and NASDAQ 100 for a specific number of previous days is calulated. If the cuotient for a certain day is greater than this mean plus twice its standard deviation, S&P 500 is bought and NASDAQ 100 is shorted. If it is less than the mean substracted twice the standard deviation, we short S&P 500, and buy NASDAQ 100.

The strategy is very simple, and its objective is to show how a simple version of the pair trading strategy could be implemented in C++. More sophisticated pair trading algorithms can be implemented, to try to obtain better results.

# Results

After running the code from 2006/12/18 to 2016/06/14, the return was 18.14%, that is, an average of 3.1% per year, more or less. The plot shows the evolution of the returns, with an initial capital of 1000$:

![alt text](https://raw.githubusercontent.com/imanolperez/pair-trading/master/Returns.PNG "Returns")


As we see, the returns are not spectacular, but the returns are quite stable and consistent. Moreover, as mentioned, the strategy is market-neutral as the beta must be almost zero, due to the correlation between both indices.
