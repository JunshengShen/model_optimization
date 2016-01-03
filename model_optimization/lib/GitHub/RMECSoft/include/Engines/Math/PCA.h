template <class T>
class PCA
{
public:
    void InitFromDensePoints(const suVector<const T*> &Points, UINT Dimension);
    //void InitFromPointMatrix(const SparseMatrix<T> &B);
    void InitFromPointMatrix(DenseMatrix<T> &B);
    void InitFromMATLAB(UINT Dimension, const String &EigensuVectorFilename, const String &EigenvalueFilename, const String &MeanFilename);
	
    UINT ReducedDimension(double EnergyPercent);
	
    void Transform(suVector<T> &Result, const suVector<T> &Input, UINT ReducedDimension);
	void InverseTransform(suVector<T> &Result, const suVector<T> &Input);
	
    void Transform(T *Result, const T *Input, UINT ReducedDimension);
	void InverseTransform(T *Result, const T *Input, UINT ReducedDimension);

private:
    void InitFromCorrelationMatrix(const DenseMatrix<T> &M);
    void FinalizeFromEigenSystem();

	suVector<T> _Means;
	suVector<T> _Eigenvalues;
	DenseMatrix<T> _EigensuVectors;
};

#include "PCA.cpp"
