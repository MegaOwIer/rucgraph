# rucgraph

This is a library of cpp header files mainly used for graph computing tasks.


## How to use rucgraph

You can just download the whole repository into a rucgraph folder, and add the path of this folder when compiling cpp codes.

### An example of using rucgraph on a Linux server

1. Download the whole rucgraph folder, and save it on the server, e.g., the path of the saved folder is "root/rucgraph".

2. Given a cpp file named as "try.cpp", the contents in which are
```
using namespace std;

#include <data_structures/PairingHeapYS_with_offset.h>

int main()
{
	example_PairingHeapYS_with_offset();
}
```
, where "PairingHeapYS_with_offset.h" is a header file in rucgraph. In the terminal, compile and run the above codes using the following commands:
```
g++ -std=c++17 -I/root/rucgraph try.cpp a.out
./a.out
```
, where "-I/root/rucgraph" is to add the path of the rucgraph folder when compiling. Then, you successfully run an example application of an augmented pairing heap, detailed in "PairingHeapYS_with_offset.h".
![MarineGEO circle logo](/assets/images/202212171254231.png "MarineGEO logo")

## Notes for updating files.

1. It is preferable to add an example in the end of each h file, for describing how to use the codes in the file.



## Catalog of rucgraph

### Folder: build_in_progress

This folder contains informal codes.

### Folder: data_structures

This folder contains some special data structures. 

#### PairingHeapYS_with_offset.h

This file contains an augmented pairing heap. In this heap, there is an offset value for every inside node. Using these values, we can change the key values of all inside nodes in O(1) time!