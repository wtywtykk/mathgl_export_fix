/***************************************************************************
 * random.cpp is part of Math Graphic Library                              *
 * Copyright (C) 2020-??? Diego Sejas Viscarra <dsejas.math@pm.me>,        *
 *                        Alexey Balakin <mathgl.abalakin@gmail.ru>        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Library General Public License as       *
 *   published by the Free Software Foundation; either version 3 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include <math.h>
#include <mgl2/data.h>
#include <mgl2/parser.h>
//-----------------------------------------------------------------------------
//
//	Basic random functions
//
//-----------------------------------------------------------------------------
mreal MGL_EXPORT mgl_rnd_integer(long lo, long hi)
{	return round((hi - lo) * mgl_rnd() + lo);	}
double MGL_EXPORT mgl_rnd_integer_(int *lo, int *hi)
{	return mgl_rnd_integer(*lo, *hi);	}
//-----------------------------------------------------------------------------
mreal MGL_EXPORT mgl_rnd_uniform(mreal lo, mreal hi)
{	return (hi - lo) * mgl_rnd() + lo;	}
double MGL_EXPORT mgl_rnd_uniform_(double *lo, double *hi)
{	return mgl_rnd_uniform(*lo,*hi);	}
//-----------------------------------------------------------------------------
mreal MGL_EXPORT mgl_rnd_bernoulli(mreal p)
{	return (mgl_rnd() < p)? 1:0;	}
double MGL_EXPORT mgl_rnd_bernoulli_(double *p)
{	return mgl_rnd_bernoulli(*p);	}
//-----------------------------------------------------------------------------
long MGL_EXPORT mgl_rnd_binomial(long trials, mreal p)
{
	long heads=0;
	for(long i=0; i<trials; i++)
	{
		if(mgl_rnd() < p)  heads++;
	}
	return heads;
}
int MGL_EXPORT mgl_rnd_binomial_(int *trials, double *p)
{	return mgl_rnd_binomial(*trials, *p);	}
//-----------------------------------------------------------------------------
mreal MGL_EXPORT mgl_rnd_gaussian(mreal mu, mreal sigma)
{
	mreal x, y, r=0.0;
	while (r >= 1 || r == 0)
	{
		x = 2.0 * mgl_rnd() - 1.0;
		y = 2.0 * mgl_rnd() - 1.0;
		r = x*x + y*y;
	}
	return mu + sigma * x * sqrt(-2.0 * log(r) / r);
}
double MGL_EXPORT mgl_rnd_gaussian_(double *mu, double *sigma)
{	return mgl_rnd_gaussian(*mu, *sigma);	}
//-----------------------------------------------------------------------------
mreal MGL_EXPORT mgl_rnd_exponential(mreal lambda)
{	return -log(1.0 - mgl_rnd()) / lambda;	}
double MGL_EXPORT mgl_rnd_exponential_(double *lambda)
{	return mgl_rnd_exponential(*lambda);	}
//-----------------------------------------------------------------------------
long MGL_EXPORT mgl_rnd_discrete(HCDT A) // this assumes A to be 1d
{
	long i, n=A->GetNx();
	mreal amax=0.0, sum_prob=0.0;

	for(long i=0; i<n; i++) amax += A->v(i);

	mreal r=amax*mgl_rnd();
	for(i=0; i<n; i++)
	{
		sum_prob += A->v(i);
		if(sum_prob > r)   break;
	}
	return i;
}
double MGL_EXPORT mgl_rnd_discrete_(uintptr_t *d)
{	return mgl_rnd_discrete(_DT_);	}
//-----------------------------------------------------------------------------
//
//	Data filling randomly
//
//-----------------------------------------------------------------------------
void MGL_EXPORT mgl_data_rnd_integer(HMDT d, long lo, long hi)
{
	long n = d->GetNN();
	for (long i=0; i<n; i++)    d->a[i] = mgl_rnd_integer(lo, hi);
}
void MGL_EXPORT mgl_data_rnd_integer_(uintptr_t *d, int *lo, int *hi)
{	mgl_data_rnd_integer(_DT_,*lo,*hi);	}
//-----------------------------------------------------------------------------
void MGL_EXPORT mgl_data_rnd_uniform(HMDT d, mreal lo, mreal hi)
{
	long n = d->GetNN();
	for (long i=0; i<n; i++)    d->a[i] = mgl_rnd_uniform(lo, hi);
}
void MGL_EXPORT mgl_data_rnd_uniform_(uintptr_t *d, double *lo, double *hi)
{	mgl_data_rnd_uniform(_DT_,*lo,*hi);	}
//-----------------------------------------------------------------------------
void MGL_EXPORT mgl_data_rnd_bernoulli(HMDT d, mreal p=0.5)
{
	long n = d->GetNN();
	for (long i=0; i<n; i++)    d->a[i] = mgl_rnd_bernoulli(p);
}
void MGL_EXPORT mgl_data_rnd_bernoulli_(uintptr_t *d, double *p)
{	mgl_data_rnd_bernoulli(_DT_,*p);	}
//-----------------------------------------------------------------------------
void MGL_EXPORT mgl_data_rnd_binomial(HMDT d, long trials, mreal p=0.5)
{
	long n = d->GetNN();
	for (long i=0; i<n; i++)    d->a[i] = mgl_rnd_binomial(trials, p);
}
void MGL_EXPORT mgl_data_rnd_binomial_(uintptr_t *d, double *p)
{	mgl_data_rnd_binomial(_DT_,*p);	}
//-----------------------------------------------------------------------------
void MGL_EXPORT mgl_data_rnd_gaussian(HMDT d, mreal mu=0.0, mreal sigma=1.0)
{
	long n = d->GetNN();
	for (long i=0; i<n; i++)    d->a[i] = mgl_rnd_gaussian(mu, sigma);
}
void MGL_EXPORT mgl_data_rnd_gaussian_(uintptr_t *d, double *mu, double *s)
{	mgl_data_rnd_gaussian(_DT_,*mu,*s);	}
//-----------------------------------------------------------------------------
void MGL_EXPORT mgl_data_rnd_exponential(HMDT d, mreal lambda)
{
	long n = d->GetNN();
	for (long i=0; i<n; i++)    d->a[i] = mgl_rnd_exponential(lambda);
}
void MGL_EXPORT mgl_data_rnd_exponential_(uintptr_t *d, double *l)
{	mgl_data_rnd_exponential(_DT_,*l);	}
//-----------------------------------------------------------------------------
void MGL_EXPORT mgl_data_rnd_discrete(HMDT d, HCDT a)
{
	if (!d || !a)    return;
	long m = d->GetNN(), n=a->GetNx();
	mreal amax=0;
	for(long i=0; i<n; i++) amax += a->v(i);
#pragma omp parallel for
	for(long j=0;j<m;j++)
	{
		mreal r=amax*mgl_rnd(), sum_prob=0;
		long i;
		for(i=0; i<n; i++)
		{
			sum_prob += a->v(i);
			if(sum_prob > r)   break;
		}
		d->a[j] = i;
	}
}
void MGL_EXPORT mgl_data_rnd_discrete_(uintptr_t *d, uintptr_t *a)
{	mgl_data_rnd_discrete(_DT_,_DM_(a));	}
//-----------------------------------------------------------------------------
void MGL_EXPORT mgl_shuffle(HMDT d, char dir)
{
	if(dir=='x')
	{
		long n = d->GetNx(), m = d->ny*d->nz;
		for(long i=0;i<n-1;i++)
		{
			long j = long((n-i)*mgl_rnd() + i);
			for(long k=0;k<m;k++)
			{	long ii = i+k*n, jj = j+k*n;
				mreal temp = d->a[ii];	d->a[ii] = d->a[jj];	d->a[jj] = temp;	}
		}
	}
	if(dir=='y')
	{
		long n = d->GetNy(), m = d->nx, l = d->nz;
		for(long i=0;i<n-1;i++)
		{
			long j = long((n-i)*mgl_rnd() + i);
			for(long q=0;q<l;q++)	for(long k=0;k<m;k++)
			{	long ii = k+m*(i+q*n), jj = k+m*(j+q*n);
				mreal temp = d->a[ii];	d->a[ii] = d->a[jj];	d->a[jj] = temp;	}
		}
	}
	if(dir=='z')
	{
		long n = d->GetNz(), m = d->ny*d->nx;
		for(long i=0;i<n-1;i++)
		{
			long j = long((n-i)*mgl_rnd() + i);
			for(long k=0;k<m;k++)
			{	long ii = m*i+k, jj = m*j+k;
				mreal temp = d->a[ii];	d->a[ii] = d->a[jj];	d->a[jj] = temp;	}
		}
	}
	if(dir=='a')
	{
		long n = d->GetNN();
		for(long i=0;i<n-1;i++)
		{
			long j = long((n-i)*mgl_rnd() + i);
			mreal temp = d->a[i];	d->a[i] = d->a[j];	d->a[j] = temp;
		}
	}
}
void MGL_EXPORT mgl_shuffle_(uintptr_t *d, char *dir, int)
{	mgl_shuffle(_DT_,*dir);	}
//-----------------------------------------------------------------------------
void MGL_NO_EXPORT mgl_fill_brownian(HMDT d, long n1, long n2, mreal sigma, mreal alpha)
{
	if (n1+1<n2)
	{
		long n = d->nx, nn = d->ny*d->nz, m = (n1+n2)/2;
		for(long i=0;i<nn;i++)
		{
			mreal delta = mgl_rnd_gaussian(0, sigma);
			d->a[m+n*i] = (d->a[n1+n*i]+d->a[n2+n*i])/2 + delta;
		}
		mgl_fill_brownian(d, n1, m, sigma/alpha, alpha);
		mgl_fill_brownian(d, m, n2, sigma/alpha, alpha);
	}
}
void MGL_EXPORT mgl_data_brownian(HMDT d, mreal y1, mreal y2, mreal sigma, mreal alpha)
{
	long n = d->nx, nn = d->ny*d->nz;
	for(long i=0;i<nn;i++)	{	d->a[n*i] = y1;   d->a[n*(i+1)-1] = y2;	}
	mgl_fill_brownian(d, 0, n-1, sigma, alpha);
}
void MGL_EXPORT mgl_data_brownian_(uintptr_t *d, double *y1, double *y2, double *sigma, double *alpha)
{	mgl_data_brownian(_DT_,*y1,*y2,*sigma,*alpha);	}
//-----------------------------------------------------------------------------
//
//	MGL commands
//
//-----------------------------------------------------------------------------
int MGL_NO_EXPORT mgls_data_rnd_integer(mglGraph *, long, mglArg *a, const char *k, const char *)
{
	int res=0;
	if(k[0]=='d' && a[0].d->temp)	return 5;
	mglData *d = dynamic_cast<mglData*>(a[0].d);
	if (d && !strcmp(k, "dnn"))	d->RndInteger(mgl_int(a[1].v), mgl_int(a[2].v));
	else    res = 1;
	return res;
}
//-----------------------------------------------------------------------------
int MGL_NO_EXPORT mgls_data_rnd_uniform(mglGraph *, long, mglArg *a, const char *k, const char *)
{
	int res=0;
	if(k[0]=='d' && a[0].d->temp)	return 5;
	mglData *d = dynamic_cast<mglData*>(a[0].d);
	if (d && !strcmp(k, "dnn"))	d->RndUniform(a[1].v,a[2].v);
	else    res = 1;
	return res;
}
//-----------------------------------------------------------------------------
int MGL_NO_EXPORT mgls_data_rnd_bernoulli(mglGraph *, long, mglArg *a, const char *k, const char *)
{
	int res=0;
	if(k[0]=='d' && a[0].d->temp)	return 5;
	mglData *d = dynamic_cast<mglData*>(a[0].d);
	if (d && !strcmp(k, "dn"))	d->RndBernoulli(a[1].v);
	else if (d && !strcmp(k, "d"))	d->RndBernoulli();
	else    res = 1;
	return res;
}
//-----------------------------------------------------------------------------
int MGL_NO_EXPORT mgls_data_rnd_binomial(mglGraph *, long, mglArg *a, const char *k, const char *)
{
	int res=0;
	if(k[0]=='d' && a[0].d->temp)	return 5;
	mglData *d = dynamic_cast<mglData*>(a[0].d);
	if (d && !strcmp(k, "dnn"))	d->RndBinomial(mgl_int(a[1].v),a[2].v);
	else if (d && !strcmp(k, "dn"))	d->RndBinomial(mgl_int(a[1].v));
	else    res = 1;
	return res;
}
//-----------------------------------------------------------------------------
int MGL_NO_EXPORT mgls_data_rnd_gaussian(mglGraph *, long, mglArg *a, const char *k, const char *)
{
	int res=0;
	if(k[0]=='d' && a[0].d->temp)	return 5;
	mglData *d = dynamic_cast<mglData*>(a[0].d);
	if (d && !strcmp(k, "dnn"))	d->RndGaussian(a[1].v,a[2].v);
	else if (d && !strcmp(k, "d"))	d->RndGaussian();
	else    res = 1;
	return res;
}
//-----------------------------------------------------------------------------
int MGL_NO_EXPORT mgls_data_rnd_exponential(mglGraph *, long, mglArg *a, const char *k, const char *)
{
	int res=0;
	if(k[0]=='d' && a[0].d->temp)	return 5;
	mglData *d = dynamic_cast<mglData*>(a[0].d);
	if (d && !strcmp(k, "dn"))	d->RndExponential(a[1].v);
	else    res = 1;
	return res;
}
//-----------------------------------------------------------------------------
int MGL_NO_EXPORT mgls_data_rnd_discrete(mglGraph *, long, mglArg *a, const char *k, const char *)
{
	int res=0;
	if(k[0]=='d' && a[0].d->temp)	return 5;
	mglData *d = dynamic_cast<mglData*>(a[0].d);
	if (d && !strcmp(k, "dd"))	mgl_data_rnd_discrete(d, a[1].d);
	else    res = 1;
	return res;
}
//-----------------------------------------------------------------------------
int MGL_NO_EXPORT mgls_shuffle(mglGraph *, long, mglArg *a, const char *k, const char *)
{
	int res=0;
	if(k[0]=='d' && a[0].d->temp)	return 5;
	mglData *d = dynamic_cast<mglData*>(a[0].d);
	if (d && !strcmp(k, "ds"))	d->Shuffle(a[1].s[0]);
	else if (d && !strcmp(k, "d"))	d->Shuffle();
	else    res = 1;
	return res;
}
//-----------------------------------------------------------------------------
int MGL_NO_EXPORT mgls_brownian(mglGraph *, long, mglArg *a, const char *k, const char *)
{
	int res=0;
	if(k[0]=='d' && a[0].d->temp)	return 5;
	mglData *d = dynamic_cast<mglData*>(a[0].d);
	if (d && !strcmp(k, "dnnnn"))	d->Brownian(a[1].v,a[2].v,a[3].v,a[4].v);
	else    res = 1;
	return res;
}
//-----------------------------------------------------------------------------
mglCommand mgls_rnd_cmd[] = {
	{"bernoulli", "Fills by random numbers according to Bernoulli distribution with probability p", "bernoulli A [p]", mgls_data_rnd_bernoulli, 3},
	{"binomial", "Fills by random numbers according to binomial distribution in n coin flips with probability p", "binomial A n [p]", mgls_data_rnd_binomial, 3},
	{"brownian", "Fills by fractional brownian motion", "brownian A y1 y2 sigma h", mgls_brownian, 3},
	{"discrete", "Fills by random numbers according to discrete distribution", "discrete A D", mgls_data_rnd_discrete, 3},
	{"exponential", "Fills by random numbers according to exponential distribution with probability p", "exponential A p", mgls_data_rnd_exponential, 3},
	{"gaussian", "Fills by random numbers according to Gaussian distribution", "gaussian A [mu sigma]", mgls_data_rnd_gaussian, 3},
	{"shuffle", "Shuffle data cells (for dir='a') or slices (for dir='xyz')", "shuffle A ['dir']", mgls_shuffle, 3},
	{"uniform", "Fills by random numbers uniformly chosen in [lo,hi)", "uniform A lo hi", mgls_data_rnd_uniform, 3},
	{"uniformint", "Fills by random integers uniformly chosen in [lo, hi)", "uniformint A lo hi", mgls_data_rnd_integer, 3},
{"","","",NULL,0}};
//-----------------------------------------------------------------------------
