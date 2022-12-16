#pragma once

#include <queue>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <text mining/list_all_files_in_a_directory.h>
#pragma warning(disable : 4996) // https://stackoverflow.com/questions/13550864/error-c4996-ctime-this-function-or-variable-may-be-unsafe 


void ysgraph_statistics() {

    queue<string> unchecked_directories;
    unchecked_directories.push("C:/Users/Yahui/Google Drive/ysgraph"); // path of ysgraph

    /*get all_header_files_paths inside the above path*/
    vector<string> all_header_files_paths;
    while (unchecked_directories.size() > 0) {
        string path1 = unchecked_directories.front();
        unchecked_directories.pop();
        vector<pair<string, bool>> files = list_all_files_in_a_directory(path1);
        for (int i = 0; i < files.size() && i < 3e2; i++) {
            if (files[i].second) { // sub-directory
                unchecked_directories.push(path1 + "/" + files[i].first); // an inside directory
            }
            else {
                if (files[i].first.substr(files[i].first.size() - 2, 2) == ".h") { // a header file
                    all_header_files_paths.push_back(path1 + "/" + files[i].first);
                }
            }
        }
    }

    /*header_files statistics*/
    int header_files_total_lines = 0, header_files_total_non_empty_lines = 0, header_files_total_string_length = 0;
    double per_header_file_avg_lines = 0, per_header_file_avg_non_empty_lines = 0, per_header_file_avg_string_length = 0;
    for (int i = 0; i < all_header_files_paths.size(); i++) {

        /*read a header file*/
        string line_content;
        ifstream myfile(all_header_files_paths[i]); // open the file
        int total_lines = 0, total_non_empty_lines = 0, total_string_length = 0;
        while (getline(myfile, line_content)) // read file line by line
        {
            total_lines++;
            total_string_length = total_string_length + line_content.length();
            if (line_content.length() != 0) {
                total_non_empty_lines++;
            }
        }
        myfile.close(); //close the file

        header_files_total_lines = header_files_total_lines + total_lines;
        header_files_total_non_empty_lines = header_files_total_non_empty_lines + total_non_empty_lines;
        header_files_total_string_length = header_files_total_string_length + total_string_length;
        per_header_file_avg_lines = per_header_file_avg_lines + (double)total_lines / all_header_files_paths.size();
        per_header_file_avg_non_empty_lines = per_header_file_avg_non_empty_lines + (double)total_non_empty_lines / all_header_files_paths.size();
        per_header_file_avg_string_length = per_header_file_avg_string_length + (double)total_string_length / all_header_files_paths.size();
    }


    /*print statistics*/
    time_t timetoday;
    time(&timetoday);
    cout << asctime(localtime(&timetoday));
    cout << "Total number of header files: " << all_header_files_paths.size() << endl;
    cout << "header_files_total_lines: " << header_files_total_lines << endl;
    cout << "header_files_total_non_empty_lines: " << header_files_total_non_empty_lines << endl;
    cout << "header_files_total_string_length: " << header_files_total_string_length << endl;
    cout << "per_header_file_avg_lines: " << per_header_file_avg_lines << endl;
    cout << "per_header_file_avg_non_empty_lines: " << per_header_file_avg_non_empty_lines << endl;
    cout << "per_header_file_avg_string_length: " << per_header_file_avg_string_length << endl;
}