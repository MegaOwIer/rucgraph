#pragma once


template<typename T>
void binary_save_vector(std::string path, const vector<T>& myVector)
{
    std::ofstream FILE(path, std::ios::out | std::ofstream::binary);

    // Store its size
    int size = myVector.size();
    FILE.write(reinterpret_cast<const char*>(&size), sizeof(size));

    // Store its contents
    FILE.write(reinterpret_cast<const char*>(&myVector[0]), myVector.size() * sizeof(T));
    FILE.close();
}

template<typename T>
void binary_read_vector(std::string path, vector<T>& myVector)
{
    vector<T>().swap(myVector);

    ifstream FILE(path, std::ios::in | std::ifstream::binary);

    int size = 0;
    FILE.read(reinterpret_cast<char*>(&size), sizeof(size));
    if (!FILE)
    {
        std::cout << "Unable to open file " << path << endl << "Please check the file location or file name." << endl; // throw an error message
        exit(1); // end the program
    }
    myVector.resize(size);
    T f;
    for (int k = 0; k < size; ++k) {
        FILE.read(reinterpret_cast<char*>(&f), sizeof(f));
        myVector[k] = f;
    }
    vector<T>(myVector).swap(myVector);
}

