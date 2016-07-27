/***********************************************************
	gamma.c -- ガンマ関数
	        -- ベータ関数
***********************************************************/
#include "stdafx.h"
#include <math.h>
#include "gamma.h"


#define PI      3.14159265358979324  /* $\pi$ */
#define LOG_2PI 1.83787706640934548  /* $\log 2\pi$ */
#define m_numData       8

#define B0  1                 /* 以下はBernoulli数 */
#define B1  (-1.0 / 2.0)
#define B2  ( 1.0 / 6.0)
#define B4  (-1.0 / 30.0)
#define B6  ( 1.0 / 42.0)
#define B8  (-1.0 / 30.0)
#define B10 ( 5.0 / 66.0)
#define B12 (-691.0 / 2730.0)
#define B14 ( 7.0 / 6.0)
#define B16 (-3617.0 / 510.0)

double loggamma(double m_data)  /* ガンマ関数の対数 */
{
	double v, w;

	v = 1;
	while (m_data < m_numData) {  v *= m_data;  m_data++;  }
	w = 1 / (m_data * m_data);
	return ((((((((B16 / (16 * 15))  * w + (B14 / (14 * 13))) * w
	            + (B12 / (12 * 11))) * w + (B10 / (10 *  9))) * w
	            + (B8  / ( 8 *  7))) * w + (B6  / ( 6 *  5))) * w
	            + (B4  / ( 4 *  3))) * w + (B2  / ( 2 *  1))) / m_data
	            + 0.5 * LOG_2PI - log(v) - m_data + (m_data - 0.5) * log(m_data);
}

double gamma(double m_data)  /* ガンマ関数 */
{
	if (m_data < 0)
		return PI / (sin(PI * m_data) * exp(loggamma(1 - m_data)));
	return exp(loggamma(m_data));
}

double beta(double m_data, double y)  /* ベータ関数 */
{
	return exp(loggamma(m_data) + loggamma(y) - loggamma(m_data + y));
}


double digamma(double xx)
{
	double v, w;

	v = 0;
	while(xx < m_numData){ v +=1/xx;	xx++;	}
	w = 1/(xx * xx);
	v += ((((((((B16 / 16) * w + (B14 /14)) * w + (B12 / 12)) * w + (B10 / 10)) * w + (B8 / 8)) * w + (B6 / 6)) * w + (B4 / 4))* w + (B2 / 2)) * w + 0.5 / xx;
	return log(xx) - v;
}

