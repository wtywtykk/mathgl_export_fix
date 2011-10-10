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
#ifndef INFODIALOG_H
#define INFODIALOG_H
//-----------------------------------------------------------------------------
#include <QDialog>
#include <QMainWindow>
//-----------------------------------------------------------------------------
class QMathGL;
class mglDrawScript;
class mglVar;
class QComboBox;
class QTextEdit;
//-----------------------------------------------------------------------------
/// Display picture and short information about variable
class InfoDialog : public QDialog
{
Q_OBJECT
public:
	InfoDialog(QWidget *parent = 0);
	~InfoDialog();
	void setVar(mglVar *v);
	void setSlice(int k)	{	kz=k;	refresh();	};
	bool allowRefresh;
//	void setName(QString &name);
	void showEvent(QShowEvent *ev);

public slots:
	void refresh(bool force=false);

private:
	int kz;
	mglVar *var;
	QMathGL *mgl;
	QTextEdit *info;
	QComboBox *kind;
	mglDrawScript *draw;
};
//-----------------------------------------------------------------------------
#endif
//-----------------------------------------------------------------------------
