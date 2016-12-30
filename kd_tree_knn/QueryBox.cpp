/***************************************************************************
 *   Copyright (C) 2010 by Alok Watve   *
 *   alokkw@gmail.com   *
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

/*
 * QueryBox.cpp
 *
 *  Created on: Oct 21, 2010
 *      Author: alok
 */

#include "QueryBox.h"
#include <iostream>
using namespace std;
long QueryBox::dim = 2;

QueryBox::QueryBox()
{
    cntBox = new float[2*dim];
}

QueryBox::~QueryBox()
{
    if (cntBox != NULL)
        delete[] cntBox;
}

/**
 * This function is added for debugging only
 */
void QueryBox::print()
{
    for(long i=0;i<dim;i++)
    {
        cout<<cntBox[2*i]<<":"<<cntBox[2*i+1]<<",";
    }
    cout<<endl;
}

float *QueryBox::getCntCenter()
{
    float * center = new float[dim];
    for (long i=0;i < dim;i++)
    {
        center[i] = (cntBox[2*i] + cntBox[2*i+1])/2.;
    }
    return center;
}

QueryBox::QueryBox(const QueryBox &qb)
{
    cntBox = new float[dim];
    for(long i=0;i<dim;i++)
    {
        cntBox[i] = qb.cntBox[i];
    }
}

QueryBox & QueryBox::operator=(const QueryBox &qb)
{
    if (cntBox == NULL)
    {
        cntBox = new float[dim];
    }
    for(long i=0;i<dim;i++)
    {
        cntBox[i] = qb.cntBox[i];
    }
    return *this;
}

bool QueryBox::overlaps(QueryBox *qbox)
{
}

