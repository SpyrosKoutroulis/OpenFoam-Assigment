void Foam::lduMatrix::Amul
(
    solveScalarField& Apsi,
    const tmp<solveScalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{

    //Edited code
    //Open file
    static bool isFirstTime = true;

    std::ofstream logFile;
    logFile.open("AmulLog.txt", std::ios::app);

    if(!logFile.is_open())
    {
        Info << "Error: Could not open file" << endl;
        return;
    }
    //End of code

    solveScalar* __restrict__ ApsiPtr = Apsi.begin();
    
    const solveScalarField& psi = tpsi();
    const solveScalar* const __restrict__ psiPtr = psi.begin();

    const scalar* const __restrict__ diagPtr = diag().begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = upper().begin();
    const scalar* const __restrict__ lowerPtr = lower().begin();

    const label startRequest = Pstream::nRequests();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    //Edited Code
    //Log Inputs
    if (isFirstTime) {
    logFile << "nCells: " << diag().size() << "\n";
    logFile << "nFaces: " << upper().size() << "\n";
    logFile << "diagPtr: ";
    for (label cell=0; cell<diag().size(); cell++)
    {
        logFile << diagPtr[cell] << " ";
    }
    logFile << "\npsiPtr: ";
    for (label cell=0; cell<diag().size(); cell++)
    {
        logFile << psiPtr[cell] << " ";
    }
    logFile << "\nuPtr: ";
    for (label face=0; face<upper().size(); face++)
    {
        logFile << uPtr[face] << " ";
    }
    logFile << "\nlPtr: ";
    for (label face=0; face<upper().size(); face++)
    {
        logFile << lPtr[face] << " ";
    }
    logFile << "\nupperPtr: ";
    for (label face=0; face<upper().size(); face++)
    {
        logFile << upperPtr[face] << " ";
    }
    logFile << "\nlowerPtr: ";
    for (label face=0; face<upper().size(); face++)
    {
        logFile << lowerPtr[face] << " ";
    }
    logFile << "\n";
    }
    //end code

    const label nCells = diag().size();
    for (label cell=0; cell<nCells; cell++)
    {
        ApsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
    }


    const label nFaces = upper().size();

    for (label face=0; face<nFaces; face++)
    {
        ApsiPtr[uPtr[face]] += lowerPtr[face]*psiPtr[lPtr[face]];
        ApsiPtr[lPtr[face]] += upperPtr[face]*psiPtr[uPtr[face]];
    }


    //Edited Code
    //Log Inputs
    if (isFirstTime) {
    logFile << "ApsiPtr: ";
    for (label cell=0; cell<Apsi.size(); cell++)
    {
        logFile << ApsiPtr[cell] << " ";
    }
    logFile << "\n ---------- \n";

    //Clode the file
    logFile.close();
    isFirstTime = false;
    }
    //end code


    // Update interface interfaces
    updateMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt,
        startRequest
    );

    tpsi.clear();
}
