%module tab_interp

%{
#define SWIG_FILE_WITH_INIT

#include <stdlib.h>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <math.h>

#include "dust_temp_interp.h"
%}

%include "dust_temp_interp.h"
