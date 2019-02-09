#include<iostream>

#include<fstream>

using namespace std;

 

int main()

{

    int count = 0;

    string line;

 

    /* Creating input filestream */ 

    ifstream file("frag_1.fastq");
    ofstream outfile ("new.fastq");
	ofstream outfile2 ("new2.fastq");

    while (getline(file, line))
	{
        count++;
		//cout<<line<<endl;
 	}
    cout << "Numbers of lines in the file : " << count << endl;

    int first=0,second=0;
    int half = count/2;
   
   	ifstream file1 ("frag_1.fastq");
    while(first!=half){
    	getline(file1,line);
    	outfile << line<<endl;
    	first++;
    }

    while(getline(file1,line)){
    	outfile2<< line<<endl;
    }


    cout << "Numbers of lines in the file : " << count << endl;
outfile.close();
file.close();
outfile2.close();
    return 0;

}