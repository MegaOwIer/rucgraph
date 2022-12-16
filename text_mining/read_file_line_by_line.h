#pragma once
#include <fstream>
#include <string>

void read_file_line_by_line(string file_name, int line_num_between_getchar) {

	string line_content;
	ifstream myfile(file_name); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		int count = 0;
		while (getline(myfile, line_content)) // read file line by line
		{
			count++;
			std::cout << line_content << endl;
			if (count%line_num_between_getchar == 0) {
				getchar();
			}
		}
		myfile.close(); //close the file

		cout << "Total Line Num: " << count << endl;
	}
	else
	{
		std::cout << "Unable to open file " << file_name << endl << "Please check the file location or file name." << endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}

}