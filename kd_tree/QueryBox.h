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
 * QueryBox.h
 *
 *  Created on: Oct 21, 2010
 *      Author: alok
 */

#ifndef QUERYBOX_H_
#define QUERYBOX_H_

#include <iostream>

using namespace std;

class QueryBox
{
public:
    static long dim;
    float *cntBox;

    // Default constructor - Allocates space.
    QueryBox();

    // Destructor - Frees the space.
    ~QueryBox();

    /**
     * This function is added for debugging only
     */
    void print();

    /**
     * Returns the center of the query box
     */
    float * getCntCenter();

    /**
     * Copy constructor
     */
    QueryBox(const QueryBox &qb);

    /**
     * Assignment operator
     */
    QueryBox &operator=(const QueryBox &qb);

    /**
     * Returns true of this box overlaps with the other
     * box
     * @param qbox The other box
     */
    bool overlaps(QueryBox *qbox);
};

#endif /* QUERYBOX_H_ */
