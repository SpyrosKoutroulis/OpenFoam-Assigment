#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <iomanip>
using namespace std;

bool startsWith(const string& line, const string& prefix) {
    return line.compare(0, prefix.size(), prefix) == 0;
}

//Function to read/initialize the OpenFoam values.
void InitializeValues(const string& filename, int& nCells, int& nFaces, double*& diagPtr,
                      double*& psiPtr, int*& uPtr, int*& lPtr, double*& upperPtr, double*& lowerPtr, double*& ApsiPtr) {
    // Open the file
    ifstream file(filename);

    string line;
    int cellCount = 0, faceCount = 0;
    double* currentArray = nullptr; // Pointer to the currently parsed array
    int currentArraySize = 0;       // Size of the current array being parsed

    while (getline(file, line)) {
        istringstream iss(line);
        if (line == "----------"){
            break;
        }

        if (startsWith(line, "nCells:")) {
            //We use ignore to skip the name "nCells:"
            iss.ignore(7);
            iss >> nCells;

            //Initialize arrays for nCells
            diagPtr = new double[nCells];
            psiPtr = new double[nCells];
            ApsiPtr = new double[nCells];
        }
        else if(startsWith(line, "nFaces:")){
            iss.ignore(7);
            iss >> nFaces;

            //Initialize arrays for nFaces
            uPtr = new int[nFaces];
            lPtr = new int[nFaces];
            upperPtr = new double[nFaces];
            lowerPtr = new double[nFaces];
        }
        //Initialize each array and create a flag to know which array we will be filling.
        else if (startsWith(line, "diagPtr:")){
            currentArray = diagPtr;
            currentArraySize = nCells;
            cellCount = 0;
            iss.ignore(9);
        }
        else if(startsWith(line, "psiPtr:")) {
            currentArray = psiPtr;
            currentArraySize = nCells;
            cellCount = 0;
            iss.ignore(8);
        }
        else if(startsWith(line, "uPtr:")){
            currentArray = reinterpret_cast<double*>(uPtr);
            currentArraySize = nFaces;
            faceCount = 0;
            iss.ignore(6);
        }
        else if (startsWith(line, "lPtr:")){
            currentArray = reinterpret_cast<double*>(lPtr);
            currentArraySize = nFaces;
            faceCount = 0;
            iss.ignore(6);
        }
        else if (startsWith(line, "upperPtr:")){
            currentArray = upperPtr;
            currentArraySize = nFaces;
            faceCount = 0;
            iss.ignore(10);
        }
        else if (startsWith(line, "lowerPtr:")){
            currentArray = lowerPtr;
            currentArraySize = nFaces;
            faceCount = 0;
            iss.ignore(10);
        }
        else if (startsWith(line, "ApsiPtr:")){
            currentArray = ApsiPtr;
            currentArraySize = nCells;
            cellCount = 0;
            iss.ignore(9);
        }

        //Fill the arrays
        if (currentArray != nullptr){
            double value;
            while (iss >> value) {
                if (cellCount < currentArraySize || faceCount < currentArraySize) {
                    if (currentArraySize == nCells) {
                        currentArray[cellCount++] = value;
                    }
                    else if (currentArraySize == nFaces){
                        //Some arrays are double and some are int.
                        if (startsWith(line, "uPtr:") || startsWith(line, "lPtr:")) {
                            reinterpret_cast<int*>(currentArray)[faceCount++] = static_cast<int>(value);
                        }
                        else{
                            currentArray[faceCount++] = value;
                        }
                    }
                }
            }
        }
    }

    //Close the file
    file.close();
}

bool areClose(double a, double b, double tolerance = 1e-2) {
    return fabs(a - b) <= tolerance;
}

int main() {

    const string filename = "AmulLog.txt";

    //Number of cells and faces
    int nCells = 0;
    int nFaces = 0;

    //Pointers for dynamically allocated arrays
    double *diagPtr = nullptr, *psiPtr = nullptr, *upperPtr = nullptr, *lowerPtr = nullptr, *ApsiPtr = nullptr;
    int *uPtr = nullptr, *lPtr = nullptr;

    //Call the InitializeValues function
    InitializeValues(filename, nCells, nFaces, diagPtr, psiPtr, uPtr, lPtr, upperPtr, lowerPtr, ApsiPtr);

    double *myApsiPtr = nullptr;
    myApsiPtr = new double[nCells];

    cout << fixed << setprecision(17);

    for (int cell = 0; cell < nCells; cell++)
    {
        myApsiPtr[cell] = diagPtr[cell] * psiPtr[cell];
        //if (cell == 0 ) {cout << myApsiPtr[cell] << " " << diagPtr[cell] << " " << psiPtr[cell] << endl;}
    }


    for (int face = 0; face < nFaces; face++)
    {
        //if (face == 0 ) {cout << myApsiPtr[uPtr[face]] << endl;}
        myApsiPtr[uPtr[face]] += lowerPtr[face]* psiPtr[lPtr[face]];
        //if (face == 0 ) {cout << myApsiPtr[uPtr[face]] << " " << lowerPtr[face] << " " << psiPtr[lPtr[face]] << endl;}
        //if (face == 0 ) {cout << myApsiPtr[lPtr[face]] << endl;}
        myApsiPtr[lPtr[face]] += upperPtr[face] * psiPtr[uPtr[face]];
        //if (face == 0 ) {cout << myApsiPtr[lPtr[face]] << " " << upperPtr[face] << " " << psiPtr[uPtr[face]] << endl;}
    }

    //cout << myApsiPtr[0] << " " << ApsiPtr[0] << endl;


    bool flag = true;
    int k = 0;
    for (int i = 0; i < nCells; i++) {
        if (myApsiPtr[i] != ApsiPtr[i]){
            //cout << i << " " << myApsiPtr[i] << " " << ApsiPtr[i] << endl;
            if( k < 100 ){
                cout << i << " " << myApsiPtr[i] << " " << ApsiPtr[i] << endl;
                k++;
            }
            flag = false;
        }
    }
    cout << flag << endl;






    //Display the results
    cout << "nCells: " << nCells << ", nFaces: " << nFaces << endl;
    cout << "diagPtr[0]: " << diagPtr[0] << ", diagPtr[nCells-1]: " << diagPtr[nCells - 1] << endl;
    cout << "psiPtr[0]: " << psiPtr[0] << ", psiPtr[nCells-1]: " << psiPtr[nCells - 1] << endl;
    cout << "uPtr[0]: " << uPtr[0] << ", uPtr[nFaces-1]: " << uPtr[nFaces - 1] << endl;
    cout << "lPtr[0]: " << lPtr[0] << ", lPtr[nFaces-1]: " << lPtr[nFaces - 1] << endl;
    cout << "upperPtr[0]: " << upperPtr[0] << ", upperPtr[nFaces-1]: " << upperPtr[nFaces - 1] << endl;
    cout << "lowerPtr[0]: " << lowerPtr[0] << ", lowerPtr[nFaces-1]: " << lowerPtr[nFaces - 1] << endl;
    cout << "ApsiPtr[0]: " << ApsiPtr[0] << ", ApsiPtr[nCells-1]: " << ApsiPtr[nCells - 1] << endl;

    //Delete the arrays
    delete[] diagPtr;
    delete[] psiPtr;
    delete[] uPtr;
    delete[] lPtr;
    delete[] upperPtr;
    delete[] lowerPtr;
    delete[] ApsiPtr;

    return 0;
}
