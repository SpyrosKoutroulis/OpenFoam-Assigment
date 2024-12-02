template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::Amul
(
    Field<Type>& Apsi,
    const tmp<Field<Type>>& tpsi
) const
{
    Type* __restrict__ ApsiPtr = Apsi.begin();

    const Field<Type>& psi = tpsi();
    const Type* const __restrict__ psiPtr = psi.begin();

    const DType* const __restrict__ diagPtr = diag().begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const LUType* const __restrict__ upperPtr = upper().begin();
    const LUType* const __restrict__ lowerPtr = lower().begin();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        true,
        interfacesUpper_,
        psi,
        Apsi
    );

    const label nCells = diag().size();
    for (label cell = 0; cell < nCells; cell++)
    {
        ApsiPtr[cell] = dot(diagPtr[cell], psiPtr[cell]);
    }


    const label nFaces = upper().size();
    for (label face = 0; face < nFaces; face++)
    {
        ApsiPtr[uPtr[face]] += dot(lowerPtr[face], psiPtr[lPtr[face]]);
        ApsiPtr[lPtr[face]] += dot(upperPtr[face], psiPtr[uPtr[face]]);
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        true,
        interfacesUpper_,
        psi,
        Apsi
    );

    tpsi.clear();
}