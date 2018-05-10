/*  Function to count line breaks in file
    Copyright (C) 2018  Benedict Kalus
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef COUNTLINES
#define COUNTLINES

#include <stdio.h>

double countlines(FILE* File){
	// save current position of carriage
	long pos=ftell(File);
	// set carriage to beginning of file
	fseek(File, 0L, SEEK_SET);
	int numlines=0;
	// loop through file and count line brakes found
	while(!feof(File)) {
		double ch = fgetc(File);
		if(ch == '\n') {
			++numlines;
		}
	}
	// set carriage to initial position
	fseek(File, pos, SEEK_SET);
	return numlines;
};

#endif
