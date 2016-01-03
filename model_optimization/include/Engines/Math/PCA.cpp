#pragma once

template<class T>
void PCA<T>::InitFromPointMatrix(DenseMatrix<T> &B)
{
    const UINT PointCount = B.ColCount();
    const UINT Dimension = B.RowCount();

    Console::WriteLine(String("Initializing PCA, ") + String(B.ColCount()) + String(" points, ") + String(B.RowCount()) + String(" dimensions"));
	
    _Means.Allocate(Dimension);
	_Means.Clear(0.0);
	for(UINT PointIndex = 0; PointIndex < PointCount; PointIndex++)
	{
		for(UINT DimensionIndex = 0; DimensionIndex < Dimension; DimensionIndex++)
		{
            _Means[DimensionIndex] += B.Cell(DimensionIndex, PointIndex);
		}
	}
	for(UINT DimensionIndex = 0; DimensionIndex < Dimension; DimensionIndex++)
	{
		_Means[DimensionIndex] /= PointCount;
	}
	for(UINT PointIndex = 0; PointIndex < PointCount; PointIndex++)
	{
        for(UINT DimensionIndex = 0; DimensionIndex < Dimension; DimensionIndex++)
		{
		    B.Cell(DimensionIndex, PointIndex) -= _Means[DimensionIndex];
        }
	}

    DenseMatrix<T> C;
	Console::WriteLine("Building correlation matrix...");
	DenseMatrix<T>::MultiplyMMTranspose(C, B);
	DenseMatrix<T>::MultiplyInPlace(C, T(1.0) / T(PointCount));
    InitFromCorrelationMatrix(C);
}

template<class T>
void PCA<T>::InitFromDensePoints(const suVector<const T*> &Points, UINT Dimension)
{
	Console::WriteLine(String("Initializing PCA, ") + String(Points.Length()) + String(" points, ") + String(Dimension) + String(" dimensions"));
	
    const UINT PointCount = Points.Length();
    DenseMatrix<T> B(Dimension, PointCount);
	
    _Means.Allocate(Dimension);
	_Means.Clear(0.0);
	for(UINT PointIndex = 0; PointIndex < PointCount; PointIndex++)
	{
		const T *CurPoint = Points[PointIndex];
		for(UINT DimensionIndex = 0; DimensionIndex < Dimension; DimensionIndex++)
		{
			_Means[DimensionIndex] += CurPoint[DimensionIndex];
		}
	}
	for(UINT DimensionIndex = 0; DimensionIndex < Dimension; DimensionIndex++)
	{
		_Means[DimensionIndex] /= PointCount;
	}
	for(UINT PointIndex = 0; PointIndex < PointCount; PointIndex++)
	{
		const T *CurPoint = Points[PointIndex];
		for(UINT DimensionIndex = 0; DimensionIndex < Dimension; DimensionIndex++)
		{
			B[DimensionIndex][PointIndex] = CurPoint[DimensionIndex] - _Means[DimensionIndex];
		}
	}

	DenseMatrix<T> C;
	Console::WriteLine("Building correlation matrix...");
	DenseMatrix<T>::MultiplyMMTranspose(C, B);
	DenseMatrix<T>::MultiplyInPlace(C, T(1.0) / T(PointCount));
    InitFromCorrelationMatrix(C);
}

template<class T>
void PCA<T>::InitFromCorrelationMatrix(const DenseMatrix<T> &C)
{
    const UINT Dimension = C.RowCount();
	const bool OutputMatlabAndExit = false;
    if(OutputMatlabAndExit)
    {
        Console::WriteLine("Outputting MATLAB matrix and exiting...");
        C.SaveToMATLAB("C:\\MATLAB7\\work\\PCA.txt");
        _Means.SaveToASCIIFile("C:\\JReaderData\\RawMeansuVector.txt");
        exit(0);
    }
    else
    {
        Console::WriteLine("Computing eigensystem...");
        C.EigenSystem(_Eigenvalues, _EigensuVectors);
        Console::WriteLine(C.EigenTest(_Eigenvalues, _EigensuVectors));

        //_EigensuVectors.SaveToMATLAB("C:\\MATLAB7\\work\\MyEigensuVectors.txt");
        //_Eigenvalues.SaveToASCIIFile("C:\\MATLAB7\\work\\MyEigenvalues.txt");
    }

    FinalizeFromEigenSystem();
}

template<class T>
void PCA<T>::InitFromMATLAB(UINT Dimension, const String &EigensuVectorFilename, const String &EigenvalueFilename, const String &MeanFilename)
{
    Console::WriteLine(String("Initializing PCA from MATLAB eigensystem files..."));
    //M = load('PCA.txt');
    //[EigensuVectors, Eigenvalues] = eig(M);
    //dlmwrite('EigensuVectors.txt', EigensuVectors, ' ');
    //dlmwrite('Eigenvalues.txt', diag(Eigenvalues), ' ');

    _Means.LoadFromASCIIFile("C:\\JReaderData\\RawMeansuVector.txt");
    _EigensuVectors.Allocate(Dimension, Dimension);
    _Eigenvalues.Allocate(Dimension);

    ifstream EigensuVectorFile(EigensuVectorFilename.CString());
    for(UINT Row = 0; Row < Dimension; Row++)
    {
        for(UINT Col = 0; Col < Dimension; Col++)
        {
            EigensuVectorFile >> _EigensuVectors.Cell(Row, Dimension - 1 - Col);
        }
    }

    ifstream EigenvalueFile(EigenvalueFilename.CString());
    for(UINT Index = 0; Index < Dimension; Index++)
    {
        EigenvalueFile >> _Eigenvalues[Dimension - 1 - Index];
    }

    //_EigensuVectors.SaveToMATLAB("C:\\MATLAB7\\work\\RepeatedEigensuVectors.txt");
    //_Eigenvalues.SaveToASCIIFile("C:\\MATLAB7\\work\\RepeatedEigenvalues.txt");

    FinalizeFromEigenSystem();
}

template<class T>
void PCA<T>::FinalizeFromEigenSystem()
{
    const UINT Dimension = _Eigenvalues.Length();
    double EigenvalueSum = 0.0;
	for(UINT DimensionIndex = 0; DimensionIndex < Dimension; DimensionIndex++)
	{
        EigenvalueSum += Math::Max(T(0.0), _Eigenvalues[DimensionIndex]);
	}
	double CumulativeEnergy = 0.0;
	for(UINT DimensionIndex = 0; DimensionIndex < Dimension; DimensionIndex++)
	{
		CumulativeEnergy += _Eigenvalues[DimensionIndex];
		Console::WriteLine(String("Energy at ") + String(DimensionIndex + 1) + String(" terms: ") + String(CumulativeEnergy / EigenvalueSum * 100.0f) + String("%"));
	}
}

template<class T>
UINT PCA<T>::ReducedDimension(double EnergyPercent)
{
	double EigenvalueSum = 0.0;
	for(UINT DimensionIndex = 0; DimensionIndex < _Eigenvalues.Length(); DimensionIndex++)
	{
		EigenvalueSum += _Eigenvalues[DimensionIndex];
	}
	double CumulativeEnergy = 0.0;
	for(UINT DimensionIndex = 0; DimensionIndex < _Eigenvalues.Length(); DimensionIndex++)
	{
		CumulativeEnergy += _Eigenvalues[DimensionIndex];
		if(CumulativeEnergy / EigenvalueSum >= EnergyPercent)
		{
			return DimensionIndex + 1;
		}
	}
	return _Eigenvalues.Length();
}

template<class T>
void PCA<T>::Transform(suVector<T> &Result, const suVector<T> &Input, UINT ReducedDimension)
{
    if(Result.Length() != ReducedDimension)
    {
	    Result.Allocate(ReducedDimension);
    }
	Transform(Result.CArray(), Input.CArray(), ReducedDimension);
}

template<class T>
void PCA<T>::InverseTransform(suVector<T> &Result, const suVector<T> &Input)
{
	const UINT Dimension = _Means.Length();
	Result.Allocate(Dimension);
	InverseTransform(Result.CArray(), Input.CArray(), Input.Length());
}

template<class T>
void PCA<T>::Transform(T *Result, const T *Input, UINT ReducedDimension)
{
	const UINT Dimension = _Means.Length();
	for(UINT Row = 0; Row < ReducedDimension; Row++)
    {
		T Total = 0.0;
		for(UINT Index = 0; Index < Dimension; Index++)
		{
			Total += _EigensuVectors[Row][Index] * (Input[Index] - _Means[Index]);
		}
		Result[Row] = Total;
    }
}

template<class T>
void PCA<T>::InverseTransform(T *Result, const T *Input, UINT ReducedDimension)
{
	UINT Dimension = _Means.Length();
	for(UINT Col = 0; Col < Dimension; Col++)
    {
		T Total = 0.0;
		for(UINT Index = 0; Index < ReducedDimension; Index++)
		{
			Total += _EigensuVectors[Index][Col] * Input[Index];
		}
		Result[Col] = Total + _Means[Col];
    }
}
