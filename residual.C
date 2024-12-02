template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::residual
(
    Field<Type>& rA,
    const Field<Type>& psi
) const
{
    Type* __restrict__ rAPtr = rA.begin();

    const Type* const __restrict__ psiPtr = psi.begin();
    const DType* const __restrict__ diagPtr = diag().begin();
    const Type* const __restrict__ sourcePtr = source().begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const LUType* const __restrict__ upperPtr = upper().begin();
    const LUType* const __restrict__ lowerPtr = lower().begin();

    // Parallel boundary initialisation.
    // Note: there is a change of sign in the coupled
    // interface update to add the contibution to the r.h.s.

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        false,          // negate interface contributions
        interfacesUpper_,
        psi,
        rA
    );

    const label nCells = diag().size();
    for (label cell = 0; cell < nCells; cell++)
    {
        rAPtr[cell] = sourcePtr[cell] - dot(diagPtr[cell], psiPtr[cell]);
    }


    const label nFaces = upper().size();
    for (label face = 0; face < nFaces; face++)
    {
        rAPtr[uPtr[face]] -= dot(lowerPtr[face], psiPtr[lPtr[face]]);
        rAPtr[lPtr[face]] -= dot(upperPtr[face], psiPtr[uPtr[face]]);
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        false,
        interfacesUpper_,
        psi,
        rA
    );
}