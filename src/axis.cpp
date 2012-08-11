/***************************************************************************
 * axis.cpp is part of Math Graphic Library
 * Copyright (C) 2007 Alexey Balakin <balakin@appl.sci-nnov.ru>            *
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
#include <time.h>
#include <ctype.h>
#include "mgl2/data.h"
#include "mgl2/canvas.h"
#include "mgl2/prim.h"
//-----------------------------------------------------------------------------
#define islog(a, b) (((a)>0 && (b)>10*(a)) || ((b)<0 && (a)<10*(b)))
#define sign(a)	((a)<0 ? -1:1)
//-----------------------------------------------------------------------------
inline struct tm* mgl_localtime (const time_t *clock, struct tm *result, bool use_utc)
{	if (!clock || !result) return NULL;
memcpy(result,use_utc?gmtime(clock):localtime(clock),sizeof(*result));
	return result;	}
//-----------------------------------------------------------------------------
long mgl_have_color(const char *stl)
{
	register long i,j=0;
	if(stl)	for(i=0;stl[i];i++)	if(strchr(MGL_COLORS,stl[i]))	j++;
	return j;
}
//-----------------------------------------------------------------------------
wchar_t *mgl_wcsdup(const wchar_t *s)
{
	wchar_t *r = (wchar_t *)malloc((wcslen(s)+1)*sizeof(wchar_t));
	memcpy(r,s,(wcslen(s)+1)*sizeof(wchar_t));
	return r;
}
//-----------------------------------------------------------------------------
void mgl_wcstrim(wchar_t *str)
{
	wchar_t *c = mgl_wcsdup(str);
	size_t n=wcslen(str), k;
	for(k=0;k<n;k++)	if(str[k]>' ')	break;
	wcscpy(c,str+k);	n = wcslen(c);
	for(k=n;k>0;k--)	if(c[k-1]>' ')	break;
	c[k] = 0;
	wcscpy(str,c);	free(c);
}
//-----------------------------------------------------------------------------
//		Ticks setup
//-----------------------------------------------------------------------------
void mglCanvas::SetAxisStl(const char *stl, const char *tck, const char *sub)
{
	if(!stl || !(*stl))		strcpy(AxisStl,"k");
	else if(strlen(stl)<32)	strcpy(AxisStl,stl);
	if(!tck || !(*tck))		strcpy(TickStl,AxisStl);
	else if(strlen(tck)<32)	strcpy(TickStl,tck);
	if(!sub || !(*sub))		strcpy(SubTStl,AxisStl);
	else if(strlen(sub)<32)	strcpy(SubTStl,sub);
}
//-----------------------------------------------------------------------------
void mglCanvas::SetTickLen(mreal tlen, mreal stt)
{	TickLen = tlen?tlen:0.02;	st_t=stt>0?stt:1;	}
//-----------------------------------------------------------------------------
void mglCanvas::SetTicks(char dir, mreal d, int ns, mreal org)
{
	if(!strchr("xyzca",dir))	return;
	mglAxis &aa = (dir=='x' ? ax : (dir=='y' ? ay : (dir=='z' ? az : ac)));

	if(aa.f==1)	aa.t[0]=0;
	aa.d=d;	aa.f=0;	aa.ns=ns;	aa.o=org;
	aa.txt.clear();
}
//-----------------------------------------------------------------------------
void mglCanvas::SetTicksVal(char dir, HCDT v, const wchar_t *lbl, bool add)
{
	if(!strchr("xyzca",dir))	return;
	mglAxis &aa = (dir=='x' ? ax : (dir=='y' ? ay : (dir=='z' ? az : ac)));
	bool ff = (dir=='x' ? fx : (dir=='y' ? fy : (dir=='z' ? fz : fa)));

	aa.txt.clear();
	if(add)	{	UpdateAxis();	AdjustTicks(aa,ff);	}
	if(!v || !lbl || !lbl[0])	{	aa.f = 0;	return;	}
	aa.f = 2;	aa.ns=0;	aa.ds=0;
	register long i,j,l=0,n=v->GetNx();
	for(i=j=l=0;i<n && lbl[j];j++)
	{
		if(lbl[j]=='\n')
		{
			aa.AddLabel(std::wstring(lbl+l,j-l),v->v(i));
			i++;	l=j+1;
		}
		if(j>1 && lbl[j]=='n' && lbl[j-1]=='\\')
		{
			aa.AddLabel(std::wstring(lbl+l,j-l-1),v->v(i));
			i++;	l=j+1;
		}
	}
	if(i<n && lbl[l])	aa.AddLabel(lbl+l,v->v(i));
}
//-----------------------------------------------------------------------------
void mglCanvas::SetTicksVal(char dir, HCDT v, const char *lbl, bool add)
{
	long ll=strlen(lbl);
	wchar_t *wcs = new wchar_t[ll+1];
	mbstowcs(wcs,lbl,ll);	wcs[ll]=0;
	SetTicksVal(dir,v,wcs,add);
	delete []wcs;
}
//-----------------------------------------------------------------------------
void mglCanvas::SetTicksVal(char dir, const wchar_t *lbl, bool add)
{
	register long i,j,len=wcslen(lbl);
	for(i=0,j=1;j<len;j++)
		if(lbl[j]=='\n' || (lbl[j]=='n' && lbl[j-1]=='\\'))	i++;
	if(i>63)	i=63;
	mglData val(i+1);	val.Fill(Min.x,Max.x);
	SetTicksVal(dir, &val, lbl, add);
}
//-----------------------------------------------------------------------------
void mglCanvas::SetTicksVal(char dir, const char *lbl, bool add)
{
	register long i,j,len=strlen(lbl);
	for(i=0,j=1;j<len;j++)
		if(lbl[j]=='\n' || (lbl[j]=='n' && lbl[j-1]=='\\'))	i++;
	if(i>63)	i=63;
	mglData val(i+1);	val.Fill(Min.x,Max.x);
	SetTicksVal(dir, &val, lbl, add);
}
//-----------------------------------------------------------------------------
void mglCanvas::SetTicksVal(char dir, HCDT v, const wchar_t **lbl, bool add)
{
	if(!strchr("xyzca",dir))	return;
	mglAxis &aa = (dir=='x' ? ax : (dir=='y' ? ay : (dir=='z' ? az : ac)));
	bool ff = (dir=='x' ? fx : (dir=='y' ? fy : (dir=='z' ? fz : fa)));

	aa.txt.clear();
	if(add)	{	UpdateAxis();	AdjustTicks(aa,ff);	}
	if(!v || !lbl)	{	aa.f = 0;	return;	}
	aa.f = 2;	aa.ns=0;	aa.ds=0;
	register long i,n=v->GetNx();
	for(i=0;i<n;i++)	aa.AddLabel(lbl[i],v->v(i));
}
//-----------------------------------------------------------------------------
void mglCanvas::SetTicksVal(char dir, HCDT v, const char **lbl, bool add)
{
	if(!strchr("xyzca",dir))	return;
	mglAxis &aa = (dir=='x' ? ax : (dir=='y' ? ay : (dir=='z' ? az : ac)));
	bool ff = (dir=='x' ? fx : (dir=='y' ? fy : (dir=='z' ? fz : fa)));

	aa.txt.clear();
	if(add)	{	UpdateAxis();	AdjustTicks(aa,ff);	}
	if(!v || !lbl)	{	aa.f = 0;	return;	}
	aa.f = 2;	aa.ns=0;	aa.ds=0;
	register unsigned long i,n=v->GetNx(),l=0;
	for(i=0;i<n;i++)	if(strlen(lbl[i])>l)	l=strlen(lbl[i]);
	wchar_t *str=new wchar_t[l+1];
	for(i=0;i<n;i++)
	{
		mbstowcs(str,lbl[i],strlen(lbl[i])+1);
		aa.AddLabel(str,v->v(i));
	}
	delete []str;
}
//-----------------------------------------------------------------------------
void mglCanvas::SetTickTempl(char dir, const wchar_t *t)
{
	if(!strchr("xyzca",dir))	return;
	mglAxis &aa = (dir=='x' ? ax : (dir=='y' ? ay : (dir=='z' ? az : ac)));

	if(aa.f==1)	aa.f = 0;	// remove time ticks
	if(!t || !t[0])	aa.t[0]=0;
	else if(wcslen(t)<255) wcscpy(aa.t,t);
}
//-----------------------------------------------------------------------------
void mglCanvas::SetTickTempl(char dir, const char *t)
{
	if(!strchr("xyzca",dir))	return;
	mglAxis &aa = (dir=='x' ? ax : (dir=='y' ? ay : (dir=='z' ? az : ac)));

	if(aa.f==1)	aa.f = 0;	// remove time ticks
	if(!t || !t[0])	aa.t[0]=0;
	else if(strlen(t)<255) mbstowcs(aa.t,t,strlen(t)+1);
}
//-----------------------------------------------------------------------------
double mgl_adj_val(double v,mreal *ds=0)
{
	double n = floor(log10(v)), s;
	v = floor(v*pow(10.,-n));	n = pow(10.,n);

	if(v==1)	{	v = n/5;	s=n/10;	}
	else if(v<4){	v = n/2;	s=n/10;	}
	else if(v<7){	v = n;		s=n/5;	}
	else		{	v = 2*n;	s=n/2;	}
	if(ds)	*ds=s;
	return v;
}
//-----------------------------------------------------------------------------
void mglCanvas::SetTickTime(char dir, mreal d, const char *t)
{
	if(!strchr("xyzca",dir))	return;
	mglAxis &aa = (dir=='x' ? ax : (dir=='y' ? ay : (dir=='z' ? az : ac)));
	UpdateAxis();

	time_t tt;	tm t1,t2;
	tt=aa.v1;	mgl_localtime(&tt, &t1, get(MGL_USE_GMTIME));
	tt=aa.v2;	mgl_localtime(&tt, &t2, get(MGL_USE_GMTIME));
	if(aa.v1<aa.v2)	// adjust periodic values
	{
		if(abs(t1.tm_year-t2.tm_year)==1)	t2.tm_yday += 365;
		if(abs(t1.tm_yday-t2.tm_yday)==1)	t2.tm_hour += 24;
		if(abs(t1.tm_hour-t2.tm_hour)==1)	t2.tm_min += 60;
		if(abs(t1.tm_min-t2.tm_min)==1)		t2.tm_sec += 60;
	}
	else
	{
		if(abs(t1.tm_year-t2.tm_year)==1)	t1.tm_yday += 365;
		if(abs(t1.tm_yday-t2.tm_yday)==1)	t1.tm_hour += 24;
		if(abs(t1.tm_hour-t2.tm_hour)==1)	t1.tm_min += 60;
		if(abs(t1.tm_min-t2.tm_min)==1)		t1.tm_sec += 60;
	}
	if(!t || !t[0])		// adjust template
	{
		t = abs(t1.tm_yday-t2.tm_yday)>1 ? "%x" : "%X";
		if(abs(t1.tm_year-t2.tm_year)>3)	t = "%Y";
	}
	if(d==0)	// try to select opimal step
	{

		if(abs(t1.tm_year-t2.tm_year)>1)
			d = 365.25*24*3600*mgl_adj_val(abs(t1.tm_year-t2.tm_year));	// number of second in year NOTE: improve it
		// NOTE here should be months ... but it is too unregular ... so omit it now
// 		else if(t1.tm_mon!=t2.tm_mon)	d = 30*24*3600;	// number of second in month
		else if(abs(t1.tm_yday-t2.tm_yday)>1)	// localtime("%x") cannot print time < 1 day
		{	d = 24*3600.*mgl_adj_val(abs(t1.tm_yday-t2.tm_yday));	d = d>24*3600?d:24*3600;	}
		else if(abs(t1.tm_hour-t2.tm_hour)>1)
			d = 3600.*mgl_adj_val(abs(t1.tm_hour-t2.tm_hour));
		else if(abs(t1.tm_min-t2.tm_min)>1)
			d = 60*mgl_adj_val(abs(t1.tm_min-t2.tm_min));
		else if(abs(t1.tm_sec-t2.tm_sec)>1)	// localtime("%X") cannot print time < 1 sec
		{	d = mgl_adj_val(abs(t1.tm_sec-t2.tm_sec));	d = d>1?d:1;	}
		else	// adjust msec. NOTE: this is not supported by localtime() !!!
			d = mgl_adj_val(fabs(aa.v2-aa.v1));
		aa.ds = 0;
	}

	aa.dv = d;	aa.f = 1;	aa.txt.clear();
	if(strlen(t)<255) mbstowcs(aa.t,t,strlen(t)+1);

	if(strchr("xyztuvw",aa.ch))
		aa.org = mglPoint(GetOrgX(aa.ch), GetOrgY(aa.ch), GetOrgZ(aa.ch));
	if(aa.ch=='x')	aa.v0 = aa.org.x;
	if(aa.ch=='y')	aa.v0 = aa.org.y;
	if(aa.ch=='z')	aa.v0 = aa.org.z;

	wchar_t buf[64];
	mreal v, v0 = mgl_isnan(aa.o) ? aa.v0 : aa.o, v1;
	if(aa.v2>aa.v1)
	{	v1 = aa.v2;		v0 = v0 - aa.dv*floor((v0-aa.v1)/aa.dv+1e-3);	}
	else
	{	v1 = aa.v1;		v0 = v0 - aa.dv*floor((v0-aa.v2)/aa.dv+1e-3);	}
	if(v0+aa.dv!=v0 && v1+aa.dv!=v1)	for(v=v0;v<=v1;v+=aa.dv)
	{
		tt = v;	tm tp;		mgl_localtime(&tt, &tp, get(MGL_USE_GMTIME));
		wcsftime(buf,64,aa.t,&tp);	aa.AddLabel(buf,v);
	}
}
//-----------------------------------------------------------------------------
void mglCanvas::AdjustTicks(const char *dir, bool force)
{
	if(force)	SetTuneTicks(-1);
	UpdateAxis();
	if(strchr(dir,'x'))	// NOTE dir have to be non-NULL here !!!
	{	if(force)	ax.d=0;	AdjustTicks(ax,fx);	}
	if(strchr(dir,'y'))
	{	if(force)	ay.d=0;	AdjustTicks(ay,fy);	}
	if(strchr(dir,'z'))
	{	if(force)	az.d=0;	AdjustTicks(az,fz);	}
	if(strchr(dir,'a') || strchr(dir,'c'))
	{	if(force)	ac.d=0;	AdjustTicks(ac,fa);	}
}
//-----------------------------------------------------------------------------
void mglCanvas::AdjustTicks(mglAxis &aa, bool ff)
{
	double d = fabs(aa.v2-aa.v1), n;
	if(aa.f>0)	return;
	if(ff && islog(aa.v1,aa.v2))
	{	aa.dv = 0;	aa.ds=0;	}
	else if(aa.d>0)
	{	aa.dv = aa.d;	aa.ds = aa.d/(abs(aa.ns)+1);	}
	else if(aa.d>-1.5)	// like =0 or =-1
	{	aa.dv = mgl_adj_val(d,&aa.ds);	aa.o=0;	}
	else
	{
		d /= -aa.d;		n = floor(log10(d));
		int k = int(d*pow(10.,-n)+0.5);
		aa.dv = pow(10.,n)*k;
		aa.o=0;	aa.ds = pow(10.,n);
	}
	LabelTicks(aa);
}
//-----------------------------------------------------------------------------
int mgl_tick_ext(mreal a, mreal b, wchar_t s[32], mreal &v)
{
	int kind = 0;
	if(fabs(a-b)<=0.001*fabs(a))
	{
		kind = 1;
		v = fabs(a-b);
		if(v>1000.f)
		{
			int k=int(log10(v)-0.01);
			kind=3;		v=mgl_ipow(10.,k);
			mglprintf(s, 32, L" (@{\\times{}10^{%d}})", k);
		}
		if(v<0.02f)
		{
			int k=int(log10(v)-0.01)-1;
			kind=3;		v=mgl_ipow(10.,k);
			mglprintf(s, 32, L" (@{\\times{}10^{%d}})", k);
		}
	}
	else
	{
		v = fabs(b)>fabs(a)?fabs(b):fabs(a);
		if(v>=1000.f)
		{
			kind = 2;
			int k=int(log10(v)-0.01);
			v=mgl_ipow(10.,k);
			mglprintf(s, 32, L" \\times 10^{%d}", k);
		}
		if(v<=1e-3f)
		{
			kind = 2;
			int k=int(log10(v)-0.01)-1;
			v=mgl_ipow(10.,k);
			mglprintf(s, 32, L" \\times 10^{%d}", k);
		}
	}
	return kind;
}
//-----------------------------------------------------------------------------
void mgl_tick_text(mreal z, mreal z0, mreal d, mreal v, int kind, wchar_t str[64], bool tune)
{
	mreal u = fabs(z)<d ? 0:z;
	if((kind&1) && z>z0)	u = fabs(z-z0)<d ? 0:(z-z0);
	if(kind==2 || (kind==3 && z>z0))	u /= v;
	if((kind&1) && z>z0)
	{
		size_t n1,n2;
		mglprintf(str, 64, L"@{(+%.2g)}",u);
//		mglprintf(str, 64, fabs(u)<1 ? L"@{(+%.2g)}" : L"@{(+%.3g)}",u);
		n1=wcslen(str);	mglprintf(str, 64, L"@{(+%g)}",u);	n2=wcslen(str);
//		if(n1<n2)	mglprintf(str, 64, fabs(u)<1 ? L"@{(+%.2g)}" : L"@{(+%.3g)}",u);
		if(n1<n2)	mglprintf(str, 64, L"@{(+%.2g)}",u);
	}
	else
	{
		size_t n1,n2;
		mglprintf(str, 64, fabs(u)<1 ? L"%.2g" :  L"%.3g",u);
		n1=wcslen(str);	mglprintf(str, 64, L"%g",u);	n2=wcslen(str);
		if(n1<n2 && tune)	mglprintf(str, 64, fabs(u)<1 ? L"%.2g" :  L"%.3g",u);
	}
}
//-----------------------------------------------------------------------------
void mglCanvas::LabelTicks(mglAxis &aa)
{
	if(strchr("xyztuvw",aa.ch))
		aa.org = mglPoint(GetOrgX(aa.ch), GetOrgY(aa.ch), GetOrgZ(aa.ch));
	if(aa.ch=='x')	aa.v0 = aa.org.x;
	if(aa.ch=='y')	aa.v0 = aa.org.y;
	if(aa.ch=='z')	aa.v0 = aa.org.z;

	wchar_t buf[64]=L"";
	mreal v,v0,v1,w;
	int d,ds;
	if(aa.f)	return;
	aa.txt.clear();
	if(aa.dv==0 && aa.v1>0)	// positive log-scale
	{
		v0 = exp(M_LN10*floor(0.1+log10(aa.v1)));
		ds = int(floor(0.1+log10(aa.v2/v0))/7)+1;
		for(v=v0;v<=aa.v2*MGL_EPSILON;v*=10)	if(v*MGL_EPSILON>=aa.v1)
		{
			d = int(floor(0.1+log10(v)));
			if(d==0)	wcscpy(buf,L"1");
			else if(d==1)	wcscpy(buf,L"10");
			else mglprintf(buf,64,L"10^{%d}",d);
			if(d%ds!=0)	*buf=0;	//	remove too often log ticks
			aa.AddLabel(buf,v);
		}
	}
	else if(aa.dv==0 && aa.v2<0)	// negative log-scale
	{
		v0 = -exp(M_LN10*floor(0.1+log10(-aa.v2)));
		ds = int(floor(0.1+log10(aa.v1/v0))/7)+1;
		for(v=v0;v>=aa.v1*MGL_EPSILON;v*=10)	if(v*MGL_EPSILON<=aa.v2)
		{
			d = int(floor(0.1+log10(-v)));
			if(d==0)	wcscpy(buf,L"-1");
			else if(d==1)	wcscpy(buf,L"-10");
			else mglprintf(buf,64,L"-10^{%d}",d);
			if(d%ds!=0)	*buf=0;	//	remove too often log ticks
			aa.AddLabel(buf,v);
		}
	}
	else if(aa.dv)	// ticks drawing
	{
		int kind=0;
		wchar_t s[32]=L"";
		if(aa.t[0]==0 && TuneTicks) kind = mgl_tick_ext(aa.v2, aa.v1, s, w);
		if((TuneTicks&1)==0 && kind==2)	kind=0;
		if((TuneTicks&2)==0 && kind!=2)	kind=0;

		v0 = mgl_isnan(aa.o) ? aa.v0 : aa.o;
		if(aa.v2>aa.v1)
		{	v1 = aa.v2;		v0 = v0 - aa.dv*floor((v0-aa.v1)/aa.dv+1e-3);	}
		else
		{	v1 = aa.v1;		v0 = v0 - aa.dv*floor((v0-aa.v2)/aa.dv+1e-3);	}

		if(v0+aa.dv!=v0 && v1+aa.dv!=v1)	for(v=v0;v<=v1;v+=aa.dv)
		{
			if(aa.t[0])
				mglprintf(buf, 64, aa.t, fabs(v)<aa.dv/100 ? 0 : v);
			else
				mgl_tick_text(v,v0,aa.dv/100,w,kind,buf,TuneTicks);
			mgl_wcstrim(buf);	aa.AddLabel(buf,v);
		}
		if(kind&2)	aa.AddLabel(s,FactorPos*(aa.v2-aa.v1)+aa.v1);
	}
}
//-----------------------------------------------------------------------------
void mglCanvas::Axis(const char *dir, const char *stl)
{
	if(!dir || !dir[0])	dir="xyz";
	bool text = !strchr(dir,'_');
	const char *ar = "AKDTVISO";
	char arr=0;
	for(unsigned i=0;i<strlen(ar);i++)
		if(strchr(dir,ar[i]))	{	arr=ar[i];	break;	}
	bool adjust = stl && strchr(stl,'a');

	bool ret = get(MGL_ENABLE_RTEXT);
	if(dir && strchr(dir,'U'))	clr(MGL_ENABLE_RTEXT);
	AdjustTicks(dir,adjust);
	// TODO: Ternary axis labeling ...
	if(strchr(dir,'x'))	DrawAxis(ax, text, arr, stl);
	if(strchr(dir,'z'))	DrawAxis(az, text, arr, stl);
	if((TernAxis&3))
	{
		mglAxis ty(ay);					ty.ch='T';
		ty.dir = mglPoint(-1,1);		ty.org = mglPoint(1,0,ay.org.z);
		DrawAxis(ty, text, arr, stl);	ty.ch='t';
		ty.dir = mglPoint(0,-1);		ty.org = mglPoint(0,1,ay.org.z);
		DrawAxis(ty, text, arr, stl);
	}
	else if(strchr(dir,'y'))	DrawAxis(ay, text, arr, stl);
	set(ret, MGL_ENABLE_RTEXT);
}
//-----------------------------------------------------------------------------
void mglCanvas::DrawAxis(mglAxis &aa, bool text, char arr,const char *stl)
{
	if(strchr("xyz",aa.ch))
		aa.org = mglPoint(GetOrgX(aa.ch), GetOrgY(aa.ch), GetOrgZ(aa.ch));
	if(aa.ch=='x')	aa.v0 = aa.org.x;
	if(aa.ch=='y')	aa.v0 = aa.org.y;
	if(aa.ch=='z')	aa.v0 = aa.org.z;

	mglPoint d = aa.dir, o = aa.org, q(NAN);	// "transverse" org
	if(strchr("xyz",aa.ch))	o -= d*(o*d);
	mglPoint av=(Min+Max)/2, dv,da,db, p;
	dv = mglPoint(sign(av.x-o.x), sign(av.y-o.y), sign(av.z-o.z));
	da = aa.a*(dv*aa.a);	db = aa.b*(dv*aa.b);

	register long i,j,k1,k2;
	SetPenPal(mgl_have_color(stl) ? stl:AxisStl);
	static int cgid=1;	StartGroup("Axis",cgid++);
	
	p = o + d*aa.v1;	k1 = AddPnt(p,CDef,q,-1,3);
	for(i=1;i<31;i++)	// axis itself
	{
		p = o + d*(aa.v1+(aa.v2-aa.v1)*i/30.);
		k2 = k1;	k1 = AddPnt(p,CDef,q,-1,3);
		line_plot(k2,k1);
	}
	if(arr)
	{
		p = o + d*(aa.v1+(aa.v2-aa.v1)*1.05);
		k2 = k1;	k1 = AddPnt(p,CDef,q,-1,3);
		arrow_plot(k1,k2,arr);
	}

	k2 = aa.txt.size();
	mreal v, u, v0 = mgl_isnan(aa.o) ? aa.v0 : aa.o;
	if(k2>0)	for(i=0;i<k2;i++)
	{
		v = aa.txt[i].val;	u = fabs(v);
		if(v>=aa.v1 && v<=aa.v2)	tick_draw(o+d*v, da, db, 0, stl);
//		else	tick_draw(o+d*v, da, db, 0, " ");
		if(aa.dv==0 && fabs(u-exp(M_LN10*floor(0.1+log10(u))))<0.01*u)
			for(j=2;j<10 && v*j<aa.v2;j++)	tick_draw(o+d*(v*j),da,db,1,stl);
	}
	if(aa.ds>0 && !get(MGL_NOSUBTICKS))
	{
		if(aa.v2>aa.v1)	v0 = v0 - aa.ds*floor((v0-aa.v1)/aa.ds+1e-3);
		else			v0 = v0 - aa.ds*floor((v0-aa.v2)/aa.ds+1e-3);
		if(v0+aa.ds!=v0 && aa.v2+aa.ds!=aa.v2)
			for(v=v0;v<aa.v2;v+=aa.ds)	tick_draw(o+d*v,da,db,1,stl);
	}
	if(text)	DrawLabels(aa);
	EndGroup();
}
//-----------------------------------------------------------------------------
void mglCanvas::DrawLabels(mglAxis &aa)
{
	if(strchr("xyz",aa.ch))
		aa.org = mglPoint(GetOrgX(aa.ch), GetOrgY(aa.ch), GetOrgZ(aa.ch));
	mglPoint d = aa.dir, o = aa.org;	// "transverse" org
	if(strchr("xyz",aa.ch))	o -= d*(o*d);
	mglPoint p,q, s=(Min+Max)/2, nn;
	s = s - d*(s*d);

	register long i,n = aa.txt.size();
	char pos[4]="t:C";
//	if(get(MGL_DISABLE_SCALE) && ((aa.dir.x==0 && aa.org.x<0) || (aa.dir.y==0 && aa.org.y>0)))	pos[0]='T';
	if(aa.ch=='c')	pos[0]=(aa.ns==0 || aa.ns==3)?'t':'T';
	if(aa.ch=='T')	pos[0]='T';
	mreal *w=new mreal[n], h = TextHeight(FontDef,-1)/3, c=NAN, l=NAN, tet=0, v, vv;	// find sizes
	long *kk=new long[n];
	for(i=0;i<n;i++)
	{
		w[i] = TextWidth(aa.txt[i].text.c_str(),FontDef,-1)/2;
		kk[i] = AddPnt(o+d*aa.txt[i].val,-1,d,0,7);
	}

	for(l=0,c=1e7,i=0;i<n-1;i++)
	{
		// exclude factors
		if(aa.ch!='c' && (aa.txt[i].val<aa.v1 || aa.txt[i+1].val<aa.v1 || aa.txt[i].val>aa.v2 || aa.txt[i+1].val>aa.v2))
			continue;
		v = (GetPntP(kk[i+1])-GetPntP(kk[i])).norm();	// distance between ticks
		vv = (w[i]+w[i+1])/2;	// length of labels
		if(v>0 && l < vv/v)	l = vv/v;
		if(c>v)	c = v;
	}
	if(get(MGL_ENABLE_RTEXT) && get(MGL_TICKS_ROTATE) && l>1 && c>0)	// try rotate first
	{	tet = c>1.1*h ? asin(1.1*h/c) : M_PI/2;	pos[2]=aa.ch=='c'?'R':'L';
		l=0.99*h/sin(tet)/c;	for(i=0;i<n;i++)	w[i]=l*c;	}
	// TODO: do clever points exclusion (i.e. longest and so on)
	long k = get(MGL_TICKS_SKIP) ? 1+l : 1;
	if(n>0)	for(i=0;i<n;i++)
	{
		if(kk[i]<0)	continue;	// should be never here?!
		c = aa.txt[i].val;
		if(c>aa.v1 && c<aa.v2 && i%k!=0)	continue;
		p = o+d*c;	nn = s-o;	ScalePoint(p,nn);
		mglPnt &qq = Pnt[kk[i]];
		mreal ux=qq.u*cos(tet) + qq.v*sin(tet), uy=qq.v*cos(tet) - qq.u*sin(tet);
		if(qq.u*nn.x+qq.v*nn.y < ux*nn.x+uy*nn.y)
		{	ux=qq.u*cos(tet) - qq.v*sin(tet);	uy=qq.v*cos(tet) + qq.u*sin(tet);	}
		qq.u = ux;	qq.v = uy;

		if(!get(MGL_ENABLE_RTEXT) && nn.x!=0)	pos[2] = nn.x<0 ? 'L':'R';
		if(aa.ch=='c' && aa.txt[i].text[0]==' ')	qq.u = qq.v = NAN;
		if(!get(MGL_DISABLE_SCALE))	pos[0]=(qq.u*nn.y-qq.v*nn.x>0) ? 'T':'t';
		if(aa.ch=='T' && pos[0]=='T')	pos[0]='t';
		if(aa.ch=='T' && pos[0]=='t')	pos[0]='T';
		text_plot(kk[i], aa.txt[i].text.c_str(), pos, -1, aa.sh+0.07,CDef,tet?false:true);
	}
	delete []w;	delete []kk;
}
//-----------------------------------------------------------------------------
void mglCanvas::tick_draw(mglPoint o, mglPoint d1, mglPoint d2, int f, const char *stl)
{
	if(TickLen==0)	return;
	// try to exclude ticks out of axis range
	if(f && ((o.x-Min.x)*(o.x-Max.x)>0 || (o.y-Min.y)*(o.y-Max.y)>0 || (o.z-Min.z)*(o.z-Max.z)>0))
		return;
	mreal v = font_factor*TickLen/sqrt(1.f+f*st_t);
	mglPoint p=o;
	long k1,k2,k3=mgl_have_color(stl);

	if(*TickStl && !f)	SetPenPal(k3 ? stl:TickStl);
	if(*SubTStl && f)	SetPenPal(k3 ? stl:SubTStl);

	ScalePoint(o, d1, false);	d1.Normalize();
	ScalePoint(p, d2, false);	d2.Normalize();
	k2 = AddPnt(p, CDef, mglPoint(NAN), 0, 0);
	p += d1*v;	k1 = AddPnt(p, CDef, mglPoint(NAN), 0, 0);
	p = o+d2*v;	k3 = AddPnt(p, CDef, mglPoint(NAN), 0, 0);
	line_plot(k1,k2);	line_plot(k2,k3);
}
//-----------------------------------------------------------------------------
void mglCanvas::Grid(const char *dir, const char *pen)
{
	if(!dir || !dir[0])	dir="xyz";
	AdjustTicks(dir,false);
	SetPenPal(pen);
	// TODO: Ternary axis labeling ...
	static int cgid=1;	StartGroup("AxisGrid",cgid++);
	if(strchr(dir,'x'))	DrawGrid(ax);
	if(strchr(dir,'y'))	DrawGrid(ay);
	if(strchr(dir,'z'))	DrawGrid(az);
	EndGroup();
}
//-----------------------------------------------------------------------------
void mglCanvas::DrawGrid(mglAxis &aa)
{
	aa.org = mglPoint(GetOrgX(aa.ch), GetOrgY(aa.ch), GetOrgZ(aa.ch));
	mglPoint d=aa.dir, da1,da2,db1,db2,oa,ob, p,q;
	da1 = aa.a*(aa.a*Min);	da2 = aa.a*(aa.a*Max);
	db1 = aa.b*(aa.b*Min);	db2 = aa.b*(aa.b*Max);
	oa  = aa.b*(aa.b*aa.org);	ob  = aa.a*(aa.a*aa.org);

	register long i,j,n=aa.txt.size(),k1,k2;
	mreal v;

	Reserve(62*n);
	if(n>0)	for(i=0;i<n;i++)
	{
		q = oa+d*aa.txt[i].val;	p = q+da1;	// lines along 'a'
		k1 = AddPnt(p,CDef);
		for(j=1;j<31;j++)
		{
			v = j/30.;
			p = q+da1*(1-v)+da2*v;
			k2 = k1;	k1 = AddPnt(p,CDef);
			line_plot(k2,k1);
		}
		q = ob+d*aa.txt[i].val;	p = q+db1;	// lines along 'b'
		k1 = AddPnt(p,CDef);
		for(j=1;j<31;j++)
		{
			v = j/30.;
			p = q+db1*(1-v)+db2*v;
			k2 = k1;	k1 = AddPnt(p,CDef);
			line_plot(k2,k1);
		}
	}
}
//-----------------------------------------------------------------------------
void mglCanvas::Label(char dir, const char *str, mreal pos, mreal shift)
{
	unsigned s = strlen(str)+1;
	wchar_t *wcs = new wchar_t[s];
	mbstowcs(wcs,str,s);
	Labelw(dir, wcs, pos, shift);
	delete []wcs;
}
//-----------------------------------------------------------------------------
void mglCanvas::Labelw(char dir, const wchar_t *text, mreal pos, mreal shift)
{
	mreal t, x0, y0, z0;
	x0 = GetOrgX(dir);	y0 = GetOrgY(dir);	z0 = GetOrgZ(dir);
	mglPoint p,ss=(Min+Max)/2,q,qq,nn;

	if(dir=='x')	//	TODO: Tern axis & colorbar labels!!!
	{
		AdjustTicks(ax,fx);
		if(ax.dv)	t = (Min.x+Max.x+pos*(Max.x-Min.x))/2;
		else	t = Min.x*pow(Max.x/Min.x, (pos+1)/2);
		nn = mglPoint(0,ss.y-y0,ss.z-z0);
		p = mglPoint(t,y0,z0);	q = mglPoint(1,0,0);	shift += ax.sh;
	}
	if(dir=='y')
	{
		AdjustTicks(ay,fy);
		if(ay.dv)	t = (Min.y+Max.y+pos*(Max.y-Min.y))/2;
		else	t = Min.y*pow(Max.y/Min.y, (pos+1)/2);
		nn = mglPoint(ss.x-x0,0,ss.z-z0);
		p = mglPoint(x0,t,z0);	q = mglPoint(0,1,0);	shift += ay.sh;
		if(TernAxis&3)
		{
			q = mglPoint(-1,1,0);	pos=-pos;
			nn = mglPoint((ss.x-x0)/2,(ss.y-y0)/2,ss.z-z0);
		}
	}
	if(dir=='t' && (TernAxis&3))
	{
		AdjustTicks(ay,fy);	pos = -pos;
		if(ay.dv)	t = (Min.y+Max.y+pos*(Max.y-Min.y))/2;
		else	t = Min.y*pow(Max.y/Min.y, (pos+1)/2);
		nn = mglPoint(ss.x-x0,0,ss.z-z0);
		p = mglPoint(x0,t,z0);	q = mglPoint(0,1,0);	shift += ay.sh;
	}
	if(dir=='z')
	{
		AdjustTicks(az,fz);
		if(az.dv)	t = (Min.z+Max.z+pos*(Max.z-Min.z))/2;
		else	t = Min.z*pow(Max.z/Min.z, (pos+1)/2);
		nn = mglPoint(ss.x-x0,ss.y-y0,0);
		p = mglPoint(x0,y0,t);	q = mglPoint(0,0,1);	shift += az.sh;
	}
	ss = p;	ScalePoint(ss,nn,false);
	char font[33],ff[3]=":C";
	if(pos<-0.2)	ff[1]='L';	if(pos>0.2)	ff[1]='R';
	strcpy(font,FontDef);	strcat(font,ff);
	strcat(font,nn.y>1e-5 || nn.x<0 ? "T":"t");
	text_plot(AddPnt(p,-1,q,0,7),text,font,-1.4,0.35+shift);
}
//-----------------------------------------------------------------------------
void mglCanvas::Box(const char *col, bool ticks)
{
	mglPoint o = Org;
	mreal tl=TickLen;
	if(!ticks)	TickLen=0;
	set(MGL_NOSUBTICKS);	Org = Min;
	static int cgid=1;	StartGroup("Box",cgid++);
	Axis("xyz_",col);
	if(TernAxis&1)
	{
		Org.z=Max.z;	Org.x=Max.x;	Axis("xz_",col);
		Org.x=Min.x;	Org.y=Max.y;	Axis("z_",col);

		mglAxis ty(ay);				ty.ch='T';
		ty.dir = mglPoint(-1,1);	ty.org = mglPoint(1,0,Max.z);
		DrawAxis(ty, false, 0,col);	ty.ch='t';
		ty.dir = mglPoint(0,-1);	ty.org = mglPoint(0,1,Max.z);
		DrawAxis(ty, false, 0,col);
	}
	else if(TernAxis&2)
	{
		mglAxis ty(az);
		ty.ch='T';	ty.a=mglPoint(1,0);	ty.b=mglPoint(-1,1);
		ty.dir = mglPoint(-1,0,1);	ty.org = mglPoint(1,0,0);
		DrawAxis(ty, false, 0,col);
		ty.ch='t';	ty.a=mglPoint(0,1);	ty.b=mglPoint(-1,1);
		ty.dir = mglPoint(0,-1,1);	ty.org = mglPoint(0,1,0);
		DrawAxis(ty, false, 0,col);
	}
	else
	{
		Org.z=Max.z;	Axis("xy_",col);
		Org = Max;		Axis("xyz_",col);
		Org.z=Min.z;	Axis("xy_",col);
		Org.x=Min.x;	DrawAxis(az,0,0,col);
		Org.x=Max.x;	Org.y=Min.y;	DrawAxis(az,0,0,col);
		if(col && strchr(col,'@'))
		{
			// edge points
			mglPoint p[8]={Min,Min,Min,Min,Max,Max,Max,Max},nan=mglPoint(NAN),oo[8];
			p[1].x=Max.x;	p[2].y=Max.y;	p[3].z=Max.z;
			p[4].x=Min.x;	p[5].y=Min.y;	p[6].z=Min.z;
			mreal zm=1e5;	int im=0;
			memcpy(oo,p,8*sizeof(mglPoint));
			for(int i=0;i<8;i++)	// find deepest point
			{
				ScalePoint(p[i],nan,false);
				if(p[i].z<zm)	{	zm=p[i].z;	im=i;	}
			}
			// now draw faces
			char color[5]="{y9}";
			register int i;	// first color used for faces, last one for edges
			for(i=0;col[i];i++)	if(strchr(MGL_COLORS,col[i]))
			{
				if(i>1 && col[i-1]=='{')	{	color[1]=col[i];	color[2]=col[i+1];	break;	}
				else	{	color[0]=col[i];	color[1]=0;	break;	}
			}
			mgl_facex(this, oo[im].x, Min.y, Min.z, Max.y-Min.y, Max.z-Min.z, color,0,0);
			mgl_facey(this, Min.x, oo[im].y, Min.z, Max.x-Min.x, Max.z-Min.z, color,0,0);
			mgl_facez(this, Min.x, Min.y, oo[im].z, Max.x-Min.x, Max.y-Min.y, color,0,0);
		}
	}
	EndGroup();
	clr(MGL_NOSUBTICKS);	Org=o;	TickLen=tl;
}
//-----------------------------------------------------------------------------
void mglCanvas::Colorbar(const char *sch)
{
	bool in = sch && strchr(sch,'I');
	mreal s=1/B.pf, x=1, y=0;
	if(sch && strchr(sch,'>'))	{	x=in?(1+s)/2:1;	y=0;	}
	if(sch && strchr(sch,'<'))	{	x=in?(1-s)/2:0;	y=0;	}
	if(sch && strchr(sch,'^'))	{	x=0;	y=in?(1+s)/2:1;	}
	if(sch && strchr(sch,'_'))	{	x=0;	y=in?(1-s)/2:0;	}
	Colorbar(sch, x, y, 1, 1);
}
//-----------------------------------------------------------------------------
void mglCanvas::Colorbar(const char *sch, mreal x, mreal y, mreal w, mreal h)
{
	bool in = sch && strchr(sch,'I');
	int where = 0;		// ‘0’ - right, ‘1’ - left, ‘2’ - above, ‘3’ - under
	if(sch && strchr(sch,'>'))	where = in?1:0;
	if(sch && strchr(sch,'<'))	where = in?0:1;
	if(sch && strchr(sch,'^'))	where = in?3:2;
	if(sch && strchr(sch,'_'))	where = in?2:3;
	if(sch && strchr(sch,'A'))	{	Push();	Identity();	}

	long n=256, s = AddTexture(sch);
	mglData v(n);
	if(ac.d || Min.c*Max.c<=0)	v.Fill(Min.c,Max.c);
	else if(Max.c>Min.c && Min.c>0)
	{	v.Fill(log(Min.c), log(Max.c));		v.Modify("exp(u)");		}
	else if(Min.c<Max.c && Max.c<0)
	{	v.Fill(log(-Min.c), log(-Max.c));	v.Modify("-exp(u)");	}
	mreal *c=new mreal[n];
	for(long i=0;i<n;i++)	c[i] = GetC(s,v.a[i]);
	colorbar(&v, c, where, x, y, w, h);
	delete []c;
	if(sch && strchr(sch,'A'))	Pop();
}
//-----------------------------------------------------------------------------
void mglCanvas::Colorbar(HCDT v, const char *sch)
{
	bool in = sch && strchr(sch,'I');
	mreal s=1/B.pf, x=1, y=0;
	if(sch && strchr(sch,'>'))	{	x=in?(1+s)/2:1;	y=0;	}
	if(sch && strchr(sch,'<'))	{	x=in?(1-s)/2:0;	y=0;	}
	if(sch && strchr(sch,'^'))	{	x=0;	y=in?(1+s)/2:1;	}
	if(sch && strchr(sch,'_'))	{	x=0;	y=in?(1-s)/2:0;	}
	Colorbar(v, sch, x, y, 1, 1);
}
//-----------------------------------------------------------------------------
void mglCanvas::Colorbar(HCDT v, const char *sch, mreal x, mreal y, mreal w, mreal h)
{
	bool in = sch && strchr(sch,'I');
	int where = 0;
	if(sch && strchr(sch,'>'))	where = in?1:0;
	if(sch && strchr(sch,'<'))	where = in?0:1;
	if(sch && strchr(sch,'^'))	where = in?3:2;
	if(sch && strchr(sch,'_'))	where = in?2:3;
	if(sch && strchr(sch,'A'))	{	Push();	Identity();	}

	mreal *c=new mreal[v->GetNx()];
	if(!mgl_have_color(sch))	sch = MGL_DEF_PAL;
	long s = AddTexture(sch);
	int nc = GetNumPal(s*256);
	mreal dc = nc>1 ? 1/(MGL_EPSILON*(nc-1)):0;
	for(long i=0;i<v->GetNx();i++)	c[i] = s+i*dc;
	colorbar(v, c, where, x, y, w, h);
	delete []c;
	if(sch && strchr(sch,'A'))	Pop();
}
//-----------------------------------------------------------------------------
void mglCanvas::colorbar(HCDT vv, const mreal *c, int where, mreal x, mreal y, mreal w, mreal h)
{
	static int cgid=1;	StartGroup("Colorbar",cgid++);
	register unsigned long i,n=vv->GetNx();
	long n1,n2,n3,n4;
	mreal d,s3=B.pf,ss=1/s3;		// NOTE: colorbar was wider ss=0.9;
	mglPoint p1,p2;

	Push();	set(MGL_DISABLE_SCALE);	B=B1;	B.pf=s3;
	x = s3*(2*x-1);	y = s3*(2*y-1);	w *= s3;	h *= s3;
	for(i=0;i<n-1;i++)
	{
		d = GetA(vv->v(i))*2-1;
		p1 = p2 = mglPoint((ss*d+1)*w+x, (ss*d+1)*h+y, s3);
		switch(where)
		{
			case 1:	p1.x = x;	p2.x = x+0.1*w;	break;
			case 2:	p1.y = y-0.1*h;	p2.y = y;	break;
			case 3:	p1.y = y;	p2.y = y+0.1*h;	break;
			default:p1.x = x-0.1*w;	p2.x = x;	break;
		}
		n1 = AddPnt(p1,c[i]);	n2 = AddPnt(p2,c[i]);
		d = GetA(vv->v(i+1))*2-1;
		p1 = p2 = mglPoint((ss*d+1)*w+x, (ss*d+1)*h+y, s3);
		switch(where)
		{
			case 1:	p1.x = x;	p2.x = x+0.1*w;	break;
			case 2:	p1.y = y-0.1*h;	p2.y = y;	break;
			case 3:	p1.y = y;	p2.y = y+0.1*h;	break;
			default:p1.x = x-0.1*w;	p2.x = x;	break;
		}
		n3 = AddPnt(p1,c[i]);	n4 = AddPnt(p2,c[i]);
		quad_plot(n1,n2,n3,n4);
	}
	if(n<64)
	{
		wchar_t buf[64];
		for(i=0;i<n;i++)
		{
			d = vv->v(i);
			mglprintf(buf,64,ac.t[0]?ac.t:(fabs(d)<1 ? L"%.2g" :  L"%.3g"),d);
			ac.AddLabel(buf,d);
		}
	}
	else	{	UpdateAxis();	AdjustTicks(ac,fa);	}
	// hint for using standard label drawing function
	SetPenPal(*TickStl ? TickStl:"k-1");
	for(i=0;i<ac.txt.size();i++)
	{
		d = ac.txt[i].val = GetA(ac.txt[i].val)*2-1;
		p1 = p2 = mglPoint((ss*d+1)*w+x, (ss*d+1)*h+y, s3);
		switch(where)
		{
			case 1:	p1.x = x;	p2.x = x+0.1*w;	break;
			case 2:	p1.y = y-0.1*h;	p2.y = y;	break;
			case 3:	p1.y = y;	p2.y = y+0.1*h;	break;
			default:p1.x = x-0.1*w;	p2.x = x;	break;
		}
		n1 = AddPnt(p1);	n2 = AddPnt(p2);
		line_plot(n1,n2);
	}
	ac.dir = mglPoint(ss*w,ss*h,0);
	ac.org = mglPoint(w+x,h+y,s3+1);
	switch(where)
	{
		case 1:	ac.dir.x = 0;	ac.org.x = x+0.1*w;	break;
		case 2:	ac.dir.y = 0;	ac.org.y = y-0.1*h;	break;
		case 3:	ac.dir.y = 0;	ac.org.y = y+0.1*h;	break;
		default:ac.dir.x = 0;	ac.org.x = x-0.1*w;	break;
	}
	ac.ns = where;	DrawLabels(ac);	// NOTE ns isn't used for colorbar
	Pop();	clr(MGL_DISABLE_SCALE);	EndGroup();
}
//-----------------------------------------------------------------------------
