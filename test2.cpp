
// #include<iostream>

// #include<fstream>

// #include <math.h>

// using namespace std;



// int main () {
// char * buffer;
// long size;
// long halfSize;
// cout<<"k"<<endl;

// ifstream infile ("frag_1.fastq",ifstream::binary);
// ofstream outfile ("new.fastq",ofstream::binary);
// ofstream outfile2 ("new2.fastq",ofstream::binary);
// cout<<"k"<<endl;

// // get size of file
// infile.seekg(0,ifstream::end);
// size=infile.tellg();
// infile.seekg(0);
// cout<<"k"<<endl;
// halfSize = static_cast<int>(floor(size/2));
// // allocate memory for file content
// cout<<"k"<<endl;

// char* buffer1 = new char[halfSize];
// cout<<"k"<<endl;

// char* buffer2 = new char[size-halfSize];
// cout<<"k"<<endl;

// // read content of infile
// infile.read (buffer1,halfSize);
// infile.read (buffer2,size-halfSize);

// // write to outfile
// outfile.write (buffer1,halfSize);
// outfile2.write (buffer2,size-halfSize);

// // release dynamically-allocated memory
// delete[] buffer;
// delete[] buffer2;

// outfile.close();
// infile.close();
// outfile2.close();
// return 0;
// }


#include <fstream>
using namespace std;

int main () {

char * buffer;
long size;

ifstream infile ("frag_1.fastq");
ofstream outfile ("new.fastq");
ofstream outfile2 ("new2.fastq");

// get size of file
cout<<"k"<<endl;
infile.seekg(0,ifstream::end);
size=infile.tellg();
infile.seekg(0);

// allocate memory for file content
cout<<"k";
buffer = new char [size/2];

// read content of infile
cout<<"k";

infile.read (buffer,size/2);

// write to outfile
cout<<"k";

outfile.write (buffer,size/2);
//outfile2.write (buffer+size/2,size);

// release dynamically-allocated memory
delete[] buffer;

outfile.close();
infile.close();
outfile2.close();
return 0;
}