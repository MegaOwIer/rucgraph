#pragma once


template<typename T>
void binary_save_vector_of_vectors(std::string path, const vector<vector<T> >& myVector)
{
    std::ofstream FILE(path, std::ios::out | std::ofstream::binary);

    // Store size of the outer vector
    int s1 = myVector.size();
    FILE.write(reinterpret_cast<const char*>(&s1), sizeof(s1));

    // Now write each vector one by one
    for (auto& v : myVector) {
        // Store its size
        int size = v.size();
        FILE.write(reinterpret_cast<const char*>(&size), sizeof(size));

        // Store its contents
        FILE.write(reinterpret_cast<const char*>(&v[0]), v.size() * sizeof(T));
    }
    FILE.close();
}

template<typename T>
void binary_read_vector_of_vectors(std::string path, vector<vector<T>>& myVector)
{
    vector<vector<T>>().swap(myVector);

    ifstream FILE(path, std::ios::in | std::ifstream::binary);

    int size = 0;
    FILE.read(reinterpret_cast<char*>(&size), sizeof(size));
    if (!FILE)
    {
        std::cout << "Unable to open file " << path << endl << "Please check the file location or file name." << endl; // throw an error message
        exit(1); // end the program
    }
    myVector.resize(size);
    for (int n = 0; n < size; ++n) {
        int size2 = 0;
        FILE.read(reinterpret_cast<char*>(&size2), sizeof(size2));
        T f;
        for (int k = 0; k < size2; ++k) {
            FILE.read(reinterpret_cast<char*>(&f), sizeof(f));
            myVector[n].push_back(f);
        }
        vector<T>(myVector[n]).swap(myVector[n]);
    }
    vector<vector<T>>(myVector).swap(myVector);
}

