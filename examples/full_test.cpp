/***************************************************************************
 * full_test.cpp is part of Math Graphic Library
 * Copyright (C) 2007-2012 Alexey Balakin <mathgl.abalakin@gmail.ru>       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include <time.h>
#include <locale.h>
#if !defined(_MSC_VER) && !defined(__BORLANDC__)
#include <getopt.h>
#endif
#include <vector>
#include "mgl2/mgl.h"
#include "mgl2/eval.h"
//-----------------------------------------------------------------------------
void mgl_create_cpp_font(HMGL gr, const wchar_t *how);
//-----------------------------------------------------------------------------
struct mglSample	/// Structure for list of samples
{
	const char *name;
	void (*func)(mglGraph*);
	const char *mgl;
};
extern mglSample samp[];
extern const char *mmgl_dat_prepare;
//-----------------------------------------------------------------------------
int mgl_cmd_smp(const void *a, const void *b)
{
	const mglSample *aa = (const mglSample *)a;
	const mglSample *bb = (const mglSample *)b;
	return strcmp(aa->name, bb->name);
}
//-----------------------------------------------------------------------------
int type = 0;
int dotest  = 0;
int width  = 800;
int height = 600;
int mini = 0;
int big  = 0;
int srnd = 0;
int use_mgl = 0;
int verbose = 0;
//-----------------------------------------------------------------------------
void mgls_prepare1d(mglData *y, mglData *y1=0, mglData *y2=0, mglData *x1=0, mglData *x2=0);
void mgls_prepare2d(mglData *a, mglData *b=0, mglData *v=0);
void mgls_prepare3d(mglData *a, mglData *b=0);
void mgls_prepare2v(mglData *a, mglData *b);
void mgls_prepare3v(mglData *ex, mglData *ey, mglData *ez);
//-----------------------------------------------------------------------------
void save(mglGraph *gr,const char *name,const char *suf);
void test(mglGraph *gr)
{
	
	mglParse par;
	par.AllowSetSize(true);
	setlocale(LC_CTYPE, "");
	par.Execute(gr,"text 0 0 'test':text 0 1 'test='1+1:text 0 -1 'test='[1,3]");
//	FILE *fp=fopen("/home/balakin/progr/mathgl-code/mathgl-2x/build/test.mgl","r");
//	par.Execute(gr,fp,true);
//	fclose(fp);
}
//-----------------------------------------------------------------------------
#if !defined(_MSC_VER) && !defined(__BORLANDC__)
static struct option longopts[] =
{
	{ "big",	no_argument,	&big,		1 },
	{ "bps",	no_argument,	&type,		8 },
	{ "help",	no_argument,	NULL,		'?' },
	{ "height",	required_argument,	NULL,	'h' },
	{ "png",	no_argument,	&type,		0 },
	{ "eps",	no_argument,	&type,		1 },
	{ "gif",	no_argument,	&type,		6 },
	{ "jpeg",	no_argument,	&type,		4 },
	{ "kind",	required_argument,	NULL,	'k' },
	{ "list",	no_argument,	NULL,		'l' },
	{ "mgl",	no_argument,	&use_mgl,	1 },
	{ "mini",	no_argument,	&mini,		1 },
	{ "none",	no_argument,	&type,		7 },
	{ "obj",	no_argument,	&type,		11 },
	{ "obj_old",no_argument,	&type,		10 },
	{ "off",	no_argument,	&type,		12 },
	{ "prc",	no_argument,	&type,		5 },
	{ "pdf",	no_argument,	&type,		9 },
	{ "solid",	no_argument,	&type,		3 },
	{ "srnd",	no_argument,	&srnd,		1 },
	{ "svg",	no_argument,	&type,		2 },
	{ "stl",	no_argument,	&type,		13 },
	{ "tex",	no_argument,	&type,		14 },
	{ "json",	no_argument,	&type,		15 },
	{ "test",	no_argument,	&dotest,	1 },
	{ "font",	no_argument,	&dotest,	2 },
	{ "thread",	required_argument,	NULL,	't' },
	{ "verbose",no_argument,	&verbose,	1 },
	{ "width",	required_argument,	NULL,	'w' },
	{ NULL,		0,				NULL,		0 }
};
//-----------------------------------------------------------------------------
void usage()
{
	puts (
//		"--png		- output png\n"
		"--width=num	- png picture width\n"
		"--height=num	- png picture height\n"
		"--mini		- png picture is 200x150\n"
		"--big		- png picture is 1920x1440\n"
		"--prc		- output prc\n"
		"--pdf		- output pdf\n"
		"--eps		- output EPS\n"
		"--eps		- output LaTeX\n"
		"--jpeg		- output JPEG\n"
		"--json		- output JSON\n"
		"--solid	- output solid PNG\n"
		"--svg		- output SVG\n"
		"--obj		- output obj/mtl\n"
		"--obj_old	- output obj/mtl in old way\n"
		"--off		- output off\n"
		"--stl		- output stl\n"
		"--none		- none output\n"
		"--srnd		- use the same random numbers in any run\n"
		"--list		- print list of sample names\n"
		"--kind=name	- produce only this sample\n"
		"--thread=num	- number of threads used\n"
		"--mgl		- use MGL scripts for samples\n"
		"--test		- perform test\n"
		"--font		- write current font as C++ file\n"
	);
}
#endif
//-----------------------------------------------------------------------------
void save(mglGraph *gr,const char *name,const char *suf="")
{
	//	return;
	char buf[128];
	printf("%s ",name);	fflush(stdout);
	switch(type)
	{
		case 1:	// EPS
			sprintf(buf,"%s%s.eps",name,suf);
			gr->WriteEPS(buf);
			break;
		case 2:	// SVG
			sprintf(buf,"%s%s.svg",name,suf);
			gr->WriteSVG(buf);	break;
		case 3:	// PNG
			sprintf(buf,"%s%s.png",name,suf);
			gr->WritePNG(buf,0,true);	break;
		case 4:	// JPEG
			sprintf(buf,"%s%s.jpg",name,suf);
			gr->WriteJPEG(buf);	break;
		case 5:	// PRC
			sprintf(buf,"%s%s.prc",name,suf);
			gr->WritePRC(buf,"",false);	break;
		case 6:	// GIF
			sprintf(buf,"%s%s.gif",name,suf);
			gr->WriteGIF(buf);	break;
		case 7:	// none
			break;
		case 8:	// EPS to PNG
			sprintf(buf,"%s%s.png",name,suf);
			gr->WritePNG(buf,0,false);
			break;
 		case 9:	// PDF
			sprintf(buf,"%s%s.prc",name,suf);
			gr->WritePRC(buf);	remove(buf);	break;
		case 10:	// old OBJ
			sprintf(buf,"%s%s.obj",name,suf);
			gr->WriteOBJold(buf);	break;
		case 11:	// OBJ
			sprintf(buf,"%s%s.obj",name,suf);
			gr->WriteOBJ(buf);	break;
		case 12:	// OFF
			sprintf(buf,"%s%s.off",name,suf);
			gr->WriteOFF(buf);	break;
		case 13:	// STL
			sprintf(buf,"%s%s.stl",name,suf);
			gr->WriteSTL(buf);	break;
		case 14:	// TeX
			sprintf(buf,"%s%s.tex",name,suf);
			gr->WriteTEX(buf);	break;
		case 15:	// JSON
			sprintf(buf,"%s%s.json",name,suf);
			gr->WriteJSON(buf);	break;
		default:// PNG (no alpha)
			sprintf(buf,"%s%s.png",name,suf);
			gr->WritePNG(buf,0,false);	break;
	}
}
//-----------------------------------------------------------------------------
int main(int argc,char **argv)
{
	const char *suf = "";
	char name[256]="", *tmp;
	int ch;
	mglGraph *gr = NULL;
	mglSample *s=samp;
#if !defined(_MSC_VER) && !defined(__BORLANDC__)
	while(( ch = getopt_long_only(argc, argv, "", longopts, NULL)) != -1)
		switch(ch)
		{
			case 0:		break;
			case 'w':	width =atoi(optarg);	break;
			case 'h':	height=atoi(optarg);	break;
			case 'k':	strcpy(name, optarg);
						tmp=strchr(name,'.');	if(tmp)	*tmp=0;
						tmp=strchr(name,'-');	if(tmp)	*tmp=0;
						break;
			case 't':	mglNumThr=atoi(optarg);	break;
			case 'l':
				while(s->name[0])	{	printf("%s ",s->name);	s++;	}
				printf("\n");	return 0;
			case '?':
			default:	usage();	return 0;
		}
#endif

	if(dotest==1)	printf("Global (before):%s\n",mglGlobalMess.c_str());
	gr = new mglGraph;	//gr->SetQuality(0);
	if(	type==11|| type==12|| type==5 || type==9)	width=height;
	if(mini)		{	gr->SetSize(190,145);	suf = "-sm";	}
	else if(big)
	{	gr->SetSize(1920,1440);	suf = "-lg";	}
	else	gr->SetSize(width,height);

	if(dotest==1)
	{
		mgl_set_test_mode(true);	test(gr);
		gr->WritePNG("test.png","",false);
		gr->WriteEPS("test.eps");
		printf("Messages:%s\n",gr->Message());
		printf("Global:%s\n",mglGlobalMess.c_str());
		delete gr;	return 0;
	}
	else if(dotest==2)
	{	mgl_create_cpp_font(gr->Self(), L"!-~,¡-ÿ,̀-̏,Α-ω,ϑ,ϕ,ϖ,ϰ,ϱ,ϵ,А-я,ℏ,ℑ,ℓ,ℜ,←-↙,∀-∯,≠-≯,⟂");
		delete gr;	return 0;	}

	if(type==15)	mini=1;	// save mini version for json
	
	if(srnd)	mgl_srnd(1);
	gr->VertexColor(false);	gr->Compression(false);
	if(name[0]==0)	while(s->name[0])	// all samples
	{
		gr->DefaultPlotParam();	gr->Clf();
		if(use_mgl)
		{
			mglParse par;
			par.AllowSetSize(true);
			setlocale(LC_CTYPE, "");
			char *buf = new char[strlen(s->mgl)+strlen(mmgl_dat_prepare)+1];
			strcpy(buf,s->mgl);		strcat(buf,mmgl_dat_prepare);
			printf("\n-------\n%s\n-------\n",verbose?buf:s->mgl);
			par.Execute(gr,buf);	delete []buf;
			const char *mess = gr->Message();
			if(*mess)	printf("Warnings: %s\n-------\n",mess);
		}
		else	s->func(gr);
		save(gr, s->name, suf);
		fflush(stdout);	s++;
	}
	else	// manual sample
	{
		mglSample tst;	tst.name=name;
		int i=0;
		for(i=0;samp[i].name[0];i++);	// determine the number of samples
		s = (mglSample *) bsearch(&tst, samp, i, sizeof(mglSample), mgl_cmd_smp);
		if(s)
		{
			gr->DefaultPlotParam();	gr->Clf();
			if(use_mgl)
			{
				mglParse par;
				par.AllowSetSize(true);
				setlocale(LC_CTYPE, "");
				char *buf = new char[strlen(s->mgl)+strlen(mmgl_dat_prepare)+1];
				strcpy(buf,s->mgl);		strcat(buf,mmgl_dat_prepare);
				printf("\n-------\n%s\n-------\n",verbose?buf:s->mgl);
				par.Execute(gr,buf);	delete []buf;
				const char *mess = gr->Message();
				if(*mess)	printf("Warnings: %s\n-------\n",mess);
			}
			else	s->func(gr);
			save(gr, s->name, suf);
		}
		else	printf("no sample %s\n",name);
	}
	printf("\n");
	delete gr;	return 0;
}
//-----------------------------------------------------------------------------
