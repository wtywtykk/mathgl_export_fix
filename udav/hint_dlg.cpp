/***************************************************************************
 *   Copyright (C) 2008 by Alexey Balakin                                  *
 *   mathgl.abalakin@gmail.com                                             *
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
#include <QSettings>
#include <QCheckBox>
#include <QLayout>
#include <QPushButton>
#include <QTextEdit>
#include "hint_dlg.h"
//-----------------------------------------------------------------------------
//
//	Hint dialog
//
//-----------------------------------------------------------------------------
#define qtr		HintDialog::tr
int numHints=-1;
QString hints[] = {
	qtr("You can rotate plot by mouse. Just press 'Rotate' toolbutton, click image and hold a mouse button: left button for rotation, right button for zooming/perspective, middle button for shifting."),
	qtr("You may quickly draw the data from file. Just use: udav 'filename.dat' in command line."),
	qtr("You can copy the current image to clipboard by pressing Ctrl-Shift-C. Later you can paste it directly into yours document or presentation."),
	qtr("You can export image into a set of format (EPS, SVG, PNG, JPEG) by pressing right mouse button inside image and selecting 'Export as ...'."),
	qtr("You can setup colors for script highlighting in Property dialog. Just select menu item 'Edit/Properties'."),
	qtr("You can save the parameter of animation inside MGL script by using comment started from '##a '."),
	qtr("New drawing never clears things drawn already. For example, you can make a surface with contour lines by calling commands 'surf' and 'cont' one after another (in any order). "),
	qtr("You can put several plots in the same image by help of commands 'subplot' or 'inplot'."),
	qtr("All indexes (of data arrays, subplots and so on) are always start from 0."),
	qtr("You can edit MGL file in any text editor. Also you can run it in console by help of commands: mgl2png, mgl2eps, mgl2svg and so on."),
	qtr("You can use command 'once on|off' for marking the block which should be executed only once. For example, this can be the block of large data reading/creating/handling. Press F9 (or menu item 'Graphics/Reload') to re-execute this block."),
	qtr("You can use command 'stop' for terminating script parsing. It is useful if you don't want to execute a part of script."),
	qtr("You can type arbitrary expression as input argument for data or number. In last case, the first value of data array is used."),
	qtr("There is powerful calculator with a lot of special functions. You can use buttons or keyboard to type the expression. Also you can use existed variables in the expression."),
	qtr("The calculator can help you to put complex expression in the script. Just type the expression (which may depend on coordinates x,y,z and so on) and put it into the script."),
	qtr("You can easely insert file or folder names, last fitted formula or numerical value of selection by using menu Edit|Insert."),
	qtr("The special dialog (Edit|Insert|New Command) help you select the command, fill its arguments and put it into the script."),
	qtr("")};

//-----------------------------------------------------------------------------
extern "C"{double mgl_rnd();}
HintDialog::HintDialog(QWidget *parent) : QDialog(parent)
{
	if(numHints<0)
		for(numHints=0;!hints[numHints].isEmpty();numHints++){};
	cur = int(mgl_rnd()*numHints);
	setWindowTitle(tr("UDAV - Hint"));
	QHBoxLayout *a;
	QPushButton *b;
	QVBoxLayout *o = new QVBoxLayout(this);
	text = new QTextEdit(this);	o->addWidget(text);
	text->setReadOnly(true);	text->setText(hints[cur]);

	start = new QCheckBox(tr("&Show at startup"), this);	o->addWidget(start);
	start->setChecked(true);

	a = new QHBoxLayout;	o->addLayout(a);
	b = new QPushButton(tr("&Prev"), this);		a->addWidget(b);
	connect(b, SIGNAL(clicked()), this, SLOT(prevClicked()));
	b = new QPushButton(tr("&Next"), this);		a->addWidget(b);
	connect(b, SIGNAL(clicked()), this, SLOT(nextClicked()));
	b = new QPushButton(tr("&Close"), this);	a->addWidget(b);
	connect(b, SIGNAL(clicked()),this, SLOT(close()));
}
//-----------------------------------------------------------------------------
HintDialog::~HintDialog()	{}
//-----------------------------------------------------------------------------
void HintDialog::nextClicked()
{	cur = (cur+1)%numHints;	text->setText(hints[cur]);	}
//-----------------------------------------------------------------------------
void HintDialog::prevClicked()
{	cur--;	if(cur<0)	cur=numHints-1;	text->setText(hints[cur]);	}
//-----------------------------------------------------------------------------
void HintDialog::closeEvent(QCloseEvent *)
{
	QSettings settings("udav","UDAV");
	settings.setPath(QSettings::IniFormat, QSettings::UserScope, "UDAV");
	settings.beginGroup("/UDAV");
	settings.setValue("/showHint", start->isChecked());
	settings.endGroup();
}
//-----------------------------------------------------------------------------
