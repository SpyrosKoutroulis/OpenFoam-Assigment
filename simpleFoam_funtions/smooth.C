template<class Type, class DType, class LUType>
void Foam::TGaussSeidelSmoother<Type, DType, LUType>::smooth
(
    const word& fieldName_,
    Field<Type>& psi,
    const LduMatrix<Type, DType, LUType>& matrix_,
    const Field<DType>& rD_,
    const label nSweeps
)
{
    Type* __restrict__ psiPtr = psi.begin();

    const label nCells = psi.size();

    Field<Type> bPrime(nCells);
    Type* __restrict__ bPrimePtr = bPrime.begin();

    const DType* const __restrict__ rDPtr = rD_.begin();

    const LUType* const __restrict__ upperPtr =
        matrix_.upper().begin();

    const LUType* const __restrict__ lowerPtr =
        matrix_.lower().begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();

    const label* const __restrict__ ownStartPtr =
        matrix_.lduAddr().ownerStartAddr().begin();


    // Parallel boundary initialisation.  The parallel boundary is treated
    // as an effective jacobi interface in the boundary.
    // Note: there is a change of sign in the coupled
    // interface update to add the contibution to the r.h.s.

    for (label sweep = 0; sweep < nSweeps; sweep++)
    {
        bPrime = matrix_.source();

        matrix_.initMatrixInterfaces
        (
            false,
            matrix_.interfacesUpper(),
            psi,
            bPrime
        );

        matrix_.updateMatrixInterfaces
        (
            false,
            matrix_.interfacesUpper(),
            psi,
            bPrime
        );

        Type curPsi;
        label fStart;
        label fEnd = ownStartPtr[0];

        for (label celli = 0; celli < nCells; celli++)
        {
            // Start and end of this row
            fStart = fEnd;
            fEnd = ownStartPtr[celli + 1];

            // Get the accumulated neighbour side
            curPsi = bPrimePtr[celli];

            // Accumulate the owner product side
            for (label curFace = fStart; curFace < fEnd; curFace++)
            {
                curPsi -= dot(upperPtr[curFace], psiPtr[uPtr[curFace]]);
            }

            // Finish current psi
            curPsi = dot(rDPtr[celli], curPsi);

            // Distribute the neighbour side using current psi
            for (label curFace = fStart; curFace < fEnd; curFace++)
            {
                bPrimePtr[uPtr[curFace]] -= dot(lowerPtr[curFace], curPsi);
            }

            psiPtr[celli] = curPsi;
        }
    }
}
