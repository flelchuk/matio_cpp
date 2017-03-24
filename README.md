# matio_cpp - a C++ 11 interface to MAT I/O

## What this is

If you need to exchange data between MATLAB and your C++ application
and using MEX is not an option, you might consider the structured
.MAT "workspace" files as a data exchange format.
Multiple variables of various types can be packaged together.
Compact binary representation and compression are properties of the
.mat file format which render it superior to text-based exchange
formats, such as CSV files.

The [matio](https://sourceforge.net/projects/matio/) library can do
the reading and writing of the .MAT format for you. The library
however, offers only a lower-level C interface.
Plenty of code must be written to get variables out to MATLAB and
the other way around.

matio_cpp is a convenience C++/11 interface to the matio library,
which can significantly reduce the effort needed to read or write
.MAT.

## What it *can* help you with

* shoving *numeric* data between MATLAB/C++ via .MAT files
* automatic selection of .MAT variable type from C++ datatype
* painless subscripting (abstraction from row-major / column-major addressing)
* do more with less code - *see the example*

## What it *can't* help you with (yet)

* cell arrays, structures, ...
* char arrays, std::string <-> char array conversion

## Usage

Before you can make use of matio_cpp, you will first need to
integrate the matio library into your project and get that
to compile and link successfully.

Once matio works, making use of matio_cpp is very simple.
Just include the matio_cpp.h file:

```C++
#include "matio_cpp.h"
```

It features the source code of all class templates, so you need
not link to additional libraries.
All features - such as classes and data types - end up in the
*MatioCPP* namespace:

```C++
using namespace MatioCPP;
```

## Example code

```C++
#include <cstdlib>
#include <cmath>
#include "matio_cpp.h"  // also does #include <matio.h>

using namespace MatioCPP;

#define NUM_ROWS 7
#define NUM_COLS 9

int main()
{
    // create .MAT file, as usual with matio
    mat_t * mat_file = Mat_CreateVer(matName, NULL, MAT_FT_MAT73);
    
    if (mat_file)
    {
        // create a two-dimensional array of 'single' (matrix)
        MatVar<float, 2> mat_of_roots("matrix_of_roots", NUM_ROWS, NUM_COLS);
        
        for (int i = 0; i < NUM_ROWS; i++)
        {
            for (int j = 0; j < NUM_COLS; j++)
            {
                float x = static_cast<float>( i + 1 );
                float y = static_cast<float>( j + 1 );
                
                // Noteworthy in the following line of code:
                // (1) use of '->' to access the MultiArray behind the variable
                // (2) at(i, j, ...) takes as many arguments as there are dimensions
                // (3) at(...) returns a lvalue (float&) which can be assigned to
                mat_of_roots->at(i, j) = std::sqrt(x * y);
            }
        }
        
        // write the variable to the MAT file
        mat_of_roots.write(mat_file);
        
        // close the .MAT file
        Mat_Close(mat_file);
        
        // mat_of_roots gets out of scope at the end of the block. It's destructor
        // also frees the buffer behind mat_of_roots.
    }
    
    return 0;
}
```

## License

See LICENSE.txt.

## TODO...

* Tutorial
* (Online?) docs
